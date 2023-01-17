/*
  Copyright 2010 Statoil ASA.

  This file is part of The Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <config.h>

#include <opm/upscaling/RelPermUtils.hpp>

#include <opm/input/eclipse/Parser/ParserKeywords/P.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/S.hpp>

#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#include <opm/input/eclipse/Units/Units.hpp>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace {
    double to_double(const std::string& x)
    {
        return std::strtod(x.c_str(), nullptr);
    }

    int to_int(const std::string& x)
    {
        return static_cast<int>(std::strtol(x.c_str(), nullptr, 10));
    }

    std::vector<double>::size_type
    sub2ind(std::array<std::vector<double>::size_type, 2>&& n,
            std::array<std::vector<double>::size_type, 3>&& ijk)
    {
        return ijk[0] + n[0]*(ijk[1] + n[1]*ijk[2]);
    }

    std::vector<std::vector<double>::size_type>
    zcorn_offset(const std::vector<double>::size_type nx,
                 const std::vector<double>::size_type ny)
    {
        using sz_t = std::vector<double>::size_type;
        auto  off  = std::vector<sz_t>{};
        off.reserve(2 * 2 * 2);

        const auto n1 = 2 * nx;
        const auto n2 = 2 * ny;

        for (sz_t k = 0; k < 2; ++k) {
        for (sz_t j = 0; j < 2; ++j) {
        for (sz_t i = 0; i < 2; ++i) {
            off.push_back(sub2ind({n1, n2}, { i, j, k }));
        }}}

        return off;
    }

    std::vector<double>::size_type
    zcorn_start(const std::vector<double>::size_type nx,
                const std::vector<double>::size_type ny,
                const std::array<int,3>&             ijk)
    {
        using sz_t = std::vector<double>::size_type;

        auto x2 = [](const sz_t x) { return 2 * x; };

        return sub2ind({ x2(nx)    , x2(ny)   },
                       { x2(ijk[0]), x2(ijk[1]), x2(ijk[2]) });
    }
} // Anonymous

class ArithmeticAverage
{
public:
    void add(const double x);

    double value() const;

private:
    double x_{ 0.0 };
    double n_{ 0.0 };
};

void ArithmeticAverage::add(const double x)
{
    x_ += x;
    n_ += 1.0;
}

double ArithmeticAverage::value() const
{
    return x_ / n_;
}

class ModelThickness
{
public:
    ModelThickness(const std::array<int, 3>&  cdim,
                   const std::vector<double>& zcorn);

    void activateCell(const std::array<int, 3>& ijk);

    double thickness() const;

private:
    using SizeType = std::vector<double>::size_type;

    const std::array<int, 3>    cdim_;
    const std::vector<double>&  zcorn_;
    const std::vector<SizeType> off_;

    std::vector<ArithmeticAverage> horizon_;
};

ModelThickness::ModelThickness(const std::array<int, 3>&  cdim,
                               const std::vector<double>& zcorn)
    : cdim_   (cdim)
    , zcorn_  (zcorn)
    , off_    (zcorn_offset(cdim[0], cdim[1]))
    , horizon_(2 * cdim[2])
{}

void ModelThickness::activateCell(const std::array<int, 3>& ijk)
{
    const auto  start = zcorn_start(cdim_[0], cdim_[1], ijk);
    auto* const layer = horizon_.data() +
        static_cast<SizeType>(2 * ijk[2]);

    const auto corners_per_layer = off_.size() / 2;

    auto ncorners = SizeType{0};

    for (const auto& corner : off_) {
        layer[ncorners++ / corners_per_layer]
            .add(zcorn_[start + corner]);
    }
}

double ModelThickness::thickness() const
{
    const auto& top = horizon_.front();
    const auto& bot = horizon_.back();

    return bot.value() - top.value();
}

namespace Opm {

static const std::vector<size_t> voigt_idx_tab = {0,4,8,5,2,1,7,6,3}; //!< Voigt-to-C index table

// Assumes that permtensor_t use C ordering.
double getVoigtValue(const SinglePhaseUpscaler::permtensor_t& K, int voigt_idx)
{
#if !defined(NDEBUG)
    OPM_ERROR_IF(not ((K.numRows() == 3) && (K.numCols() == 3)),
                 "Function getVoigtValue() is only supported "
                 "for 3-by-3 tensors");
#endif

    if (voigt_idx < 0 || voigt_idx > 8) {
        std::cerr << "Voigt index out of bounds (only 0-8 allowed)" << std::endl;
        throw std::exception();
    }

    return K.data()[voigt_idx_tab[voigt_idx]];
}

// Assumes that permtensor_t use C ordering.
void setVoigtValue(SinglePhaseUpscaler::permtensor_t& K, int voigt_idx, double val)
{
#if !defined(NDEBUG)
    OPM_ERROR_IF(not ((K.numRows() == 3) && (K.numCols() == 3)),
                 "Function setVoigtValue() is only supported "
                 "for 3-by-3 tensors.");
#endif

    if (voigt_idx < 0 || voigt_idx > 8) {
        std::cerr << "Voigt index out of bounds (only 0-8 allowed)" << std::endl;
        throw std::exception();
    }

    K.data()[voigt_idx_tab[voigt_idx]] = val;
}

RelPermUpscaleHelper::RelPermUpscaleHelper(int mpi_rank,
                                           std::map<std::string,std::string>& options_)
    : isMaster         (mpi_rank == 0)
    , anisotropic_input(false)
    , permTensor       (3, 3, nullptr)
    , tesselatedCells  (0)
    , permTensorInv    (3, 3, nullptr)
    , options          (options_)
{
    {
        const auto& fluids = options["fluids"];

        if ((fluids == "ow") || (fluids == "wo")) {
            saturationstring = "Sw";
        }
        else if ((fluids == "go") || (fluids == "og")) {
            saturationstring = "Sg";
        }
        else {
            std::stringstream str;
            str << "Fluidsystem " << fluids
                << " not valid (-fluids option). Should be ow or go";

            throw std::runtime_error(str.str());
        }
    }

    critRelpThresh = to_double(options["critRelpermThresh"]);
    doEclipseCheck = options["doEclipseCheck"] == "true";
}

void RelPermUpscaleHelper::collectResults()
{
#ifdef HAVE_MPI
   /* Step 8b: Transfer all computed data to master node.
      Master node should post a receive for all values missing,
      other nodes should post a send for all the values they have.
    */
   int mpi_rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
   MPI_Barrier(MPI_COMM_WORLD); // Not strictly necessary.
   if (isMaster) {
       // Loop over all values, receive data and put into local data structure
       for (int idx=0; idx < points; ++idx) {
           if (node_vs_pressurepoint[idx] != 0) {
               // Receive data
               if (upscaleBothPhases) {
                  std::vector<double> recvbuffer(2+2*tensorElementCount);
                   MPI_Recv(recvbuffer.data(), recvbuffer.size(), MPI_DOUBLE,
                            node_vs_pressurepoint[idx], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                   // Put received data into correct place.
                   WaterSaturation[(int)recvbuffer[0]] = recvbuffer[1];
                   for (int voigtIdx=0; voigtIdx < tensorElementCount; ++voigtIdx) {
                       PhasePerm[0][(int)recvbuffer[0]][voigtIdx] = recvbuffer[2+voigtIdx];
                   }
                   for (int voigtIdx=0; voigtIdx < tensorElementCount; ++voigtIdx) {
                       PhasePerm[1][(int)recvbuffer[0]][voigtIdx] = recvbuffer[2+tensorElementCount+voigtIdx];
                   }
               }
               else {
                   std::vector<double> recvbuffer(2+tensorElementCount);
                   MPI_Recv(recvbuffer.data(), recvbuffer.size(), MPI_DOUBLE,
                            node_vs_pressurepoint[idx], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                   // Put received data into correct place.
                   WaterSaturation[(int)recvbuffer[0]] = recvbuffer[1];
                   for (int voigtIdx=0; voigtIdx < tensorElementCount; ++voigtIdx) {
                       PhasePerm[0][(int)recvbuffer[0]][voigtIdx] = recvbuffer[2+voigtIdx];
                   }
               }
           }
       }
   }
   else {
       for (int idx=0; idx < points; ++idx) {
           if (node_vs_pressurepoint[idx] == mpi_rank) {
               // Pack and send data. C-style.
               if (upscaleBothPhases) {
                   std::vector<double> sendbuffer(2+2*tensorElementCount);
                   sendbuffer[0] = (double)idx;
                   sendbuffer[1] = WaterSaturation[idx];
                   for (int voigtIdx=0; voigtIdx < tensorElementCount; ++voigtIdx) {
                       sendbuffer[2+voigtIdx] = PhasePerm[0][idx][voigtIdx];
                   }
                   for (int voigtIdx=0; voigtIdx < tensorElementCount; ++voigtIdx) {
                       sendbuffer[2+tensorElementCount+voigtIdx] = PhasePerm[1][idx][voigtIdx];
                   }
                   MPI_Send(sendbuffer.data(), sendbuffer.size(), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
               }
               else {
                   std::vector<double> sendbuffer(2+tensorElementCount);
                   sendbuffer[0] = (double)idx;
                   sendbuffer[1] = WaterSaturation[idx];
                   for (int voigtIdx=0; voigtIdx < tensorElementCount; ++voigtIdx) {
                       sendbuffer[2+voigtIdx] = PhasePerm[0][idx][voigtIdx];
                   }
                   MPI_Send(sendbuffer.data(), sendbuffer.size(), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
               }
           }
       }
   }
#endif
}

std::vector<std::vector<double>> RelPermUpscaleHelper::getRelPerm(int phase) const
{
    SinglePhaseUpscaler::permtensor_t zeroMatrix(3,3,(double*)0);
    zero(zeroMatrix);
    std::vector<std::vector<double>> RelPermValues;
    if (isMaster) {
        RelPermValues.resize(tensorElementCount);
        // Loop over all pressure points
        for (int idx=0; idx < points; ++idx) {
            SinglePhaseUpscaler::permtensor_t phasePermTensor = zeroMatrix;
            zero(phasePermTensor);
            for (int voigtIdx = 0; voigtIdx < tensorElementCount; ++voigtIdx) {
                setVoigtValue(phasePermTensor, voigtIdx, PhasePerm[phase][idx][voigtIdx]);
            }
            SinglePhaseUpscaler::permtensor_t relPermTensor = zeroMatrix;
            prod(phasePermTensor, permTensorInv, relPermTensor);
            for (int voigtIdx = 0; voigtIdx < tensorElementCount; ++voigtIdx) {
                RelPermValues[voigtIdx].push_back(getVoigtValue(relPermTensor, voigtIdx));
            }
        }
        // If doEclipseCheck, critical saturation points should be specified by 0 relperm
        // Numerical errors and maxpermcontrast violate this even if the input has specified
        // these points
        if (doEclipseCheck) {
            for (int voigtIdx = 0; voigtIdx < tensorElementCount; ++voigtIdx) {
                int minidx;
                if (RelPermValues[voigtIdx][0] < RelPermValues[voigtIdx][points-1]) minidx = 0; else minidx = points-1;
                if (RelPermValues[voigtIdx][minidx] < critRelpThresh) {
                    RelPermValues[voigtIdx][minidx] = 0.0;
                }
                else {
                    std::stringstream str;
                    str << "Minimum upscaled relperm value for phase " << phase+1 << " is "
                        << RelPermValues[voigtIdx][minidx] << ", larger than critRelpermThresh." << std::endl
                        << " (voigtidx = " << voigtIdx << ")";
                    throw std::runtime_error(str.str());
                }
            }
        }
    }

    return RelPermValues;
}

void RelPermUpscaleHelper::upscaleSinglePhasePermeability()
{
    if (isMaster) {
        auto cellperm = SinglePhaseUpscaler::permtensor_t(3, 3, nullptr);

        const auto& ecl_idx = upscaler.grid().globalCell();

        for (decltype(ecl_idx.size())
                 i = 0, n = ecl_idx.size(); i < n; ++i)
        {
            const auto cell_idx = ecl_idx[i];

            zero(cellperm);

            if (! anisotropic_input) {
                const auto kval = std::max(perms[0][cell_idx], minSinglePhasePerm);

                cellperm(0,0) = kval;
                cellperm(1,1) = kval;
                cellperm(2,2) = kval;
            }
            else {
                cellperm(0,0) = std::max(minSinglePhasePerm, perms[0][cell_idx]);
                cellperm(1,1) = std::max(minSinglePhasePerm, perms[1][cell_idx]);
                cellperm(2,2) = std::max(minSinglePhasePerm, perms[2][cell_idx]);
            }

            upscaler.setPermeability(i, cellperm);
        }

        permTensor    = upscaler.upscaleSinglePhase();
        permTensorInv = permTensor;

        invert(permTensorInv);
    }
}

void RelPermUpscaleHelper::sanityCheckInput(const Opm::Deck& deck,
                                            const double      minPerm,
                                            const double      maxPerm,
                                            const double      minPoro)
{
    {
        using kw_poro   = ParserKeywords::PORO;
        using kw_permx  = ParserKeywords::PERMX;
        using kw_permy  = ParserKeywords::PERMY;
        using kw_permz  = ParserKeywords::PERMZ;
        using kw_satnum = ParserKeywords::SATNUM;

        // Check that we have the information we need from the eclipse file:
        if (! (deck.hasKeyword<kw_poro >() &&
               deck.hasKeyword<kw_permx>()))
        {
            throw std::runtime_error("Error: Did not find complete set of "
                                     "PORO and PERMX in ECLIPSE file.");
        }

        poros    = deck.get<kw_poro >().back().getSIDoubleData();
                        perms[0] = deck.get<kw_permx>().back().getSIDoubleData();

        EclipseGrid eg(deck);
        zcorns = eg.getZCORN();

        // Load anisotropic (only diagonal supported) input if present in grid
        if (deck.hasKeyword<kw_permy>() && deck.hasKeyword<kw_permz>()) {
            anisotropic_input = true;

            perms[1] = deck.get<kw_permy>().back().getSIDoubleData();
            perms[2] = deck.get<kw_permz>().back().getSIDoubleData();

            if (isMaster) {
                std::cout << "Info: PERMY and PERMZ present, going into "
                          << "anisotropic input mode, no J-functions\n"
                          << "  Options -relPermCurve and -jFunctionCurve "
                          << "are meaningless.\n";
            }
        }

        if (deck.hasKeyword<kw_satnum>()) {
            satnums = deck.get<kw_satnum>().back().getIntData();
        }
        else {
            if (isMaster) {
                std::cout << "SATNUM not found in input file, "
                          << "assuming only one rocktype" << std::endl;
            }

            // Initialize a default satnums-vector with only "ones" (meaning
            // only one rocktype).
            satnums.assign(poros.size(), 1);
        }
    }

    // Sanity check/fix on input for each cell:
    //
    // - Check that SATNUM are set sensibly, that is => 0 and < 1000, error
    //   if not.
    //
    // - Check that porosity is between 0 and 1, error if not.  Set to
    //   minPoro if zero or less than minPoro (due to pcmin/max computation)
    //
    // - Check that permeability is zero or positive. Error if negative.
    //   Set to minPerm if zero or less than minPerm.
    //
    // - Check maximum number of SATNUM values (can be number of rock types
    //   present)

    {
        auto it  = std::find_if(satnums.begin(), satnums.end(), [](int a) { return a < 0; });
        auto it2 = std::find_if(satnums.begin(), satnums.end(), [](int a) { return a > 1000;});

        if (it != satnums.end() || it2 != satnums.end()) {
            std::stringstream str;
            str << "satnums[";

            if (it == satnums.end()) {
                str << it2 - satnums.begin() << "] = " << *it2;
            }
            else {
                str << it - satnums.begin() << "] = " << *it;
            }

            str << ", not sane, quitting.";

            throw std::runtime_error(str.str());
        }
    }

    {
        auto it3 = std::find_if(poros.begin(), poros.end(),
                                [](double value)
                                {
                                    return (value < 0) || (value > 1);
                                });

        if (it3 != poros.end()) {
            std::stringstream str;

            str << "poros[" << it3 - poros.begin() << "] = "
                << *it3 << ", not sane, quitting.";

            throw std::runtime_error(str.str());
        }
    }

    {
        auto name       = std::string("permx");
        auto check_perm = [&name](const std::vector<double>& perm)
        {
            auto it = std::find_if(perm.begin(), perm.end(),
                                   [](double value) { return value < 0; });

            if (it != perm.end()) {
                std::stringstream str;

                str << name << "[" << it - perm.begin() << "] = "
                    << *it << ", not sane, quitting.";

                throw std::runtime_error(str.str());
            }
        };

        check_perm(perms[0]);

        if (anisotropic_input) {
            name = "permy";   check_perm(perms[1]);
            name = "permz";   check_perm(perms[2]);
        }
    }

    auto cells_truncated_from_below_poro = 0;
    {
        for (auto& poro : poros) {
            if ((! (poro < 0)) && (poro < minPoro)) {
                poro = minPoro;

                ++ cells_truncated_from_below_poro;
            }
        }
    }

    auto cells_truncated_from_below_permx = 0;
    auto cells_truncated_from_above_permx = 0;
    {
        const auto lowPerm =
            unit::convert::from(minPerm, prefix::milli*unit::darcy);

        const auto highPerm =
            unit::convert::from(maxPerm, prefix::milli*unit::darcy);

        for (auto& kx : perms[0]) {
            if ((! (kx < 0)) && (kx < lowPerm)) { // Truncate permeability from below
                ++cells_truncated_from_below_permx;

                kx = lowPerm;
            }
            else if (kx > highPerm) { // Truncate permeability from above
                ++ cells_truncated_from_above_permx;

                kx = highPerm;
            }
        }

        for (decltype(satnums.size())
                 i = 0, n = satnums.size(); i < n; ++i)
        {
            // Explicitly handle "no rock" cells, set them to minimum perm
            // and zero porosity.
            if (satnums[i] == 0) {
                perms[0][i] = lowPerm;

                if (anisotropic_input) {
                    perms[1][i] = lowPerm;
                    perms[2][i] = lowPerm;
                }

                // zero poro is fine for these cells, as they are not used
                // in pcmin/max computation.
                poros[i] = 0;
            }
        }
    }

    if (isMaster) {
        if (cells_truncated_from_below_poro > 0) {
            std::cout << "Cells with truncated porosity: "
                      << cells_truncated_from_below_poro << std::endl;
        }

        if (cells_truncated_from_below_permx > 0) {
            std::cout << "Cells with permx truncated from below: "
                      << cells_truncated_from_below_permx << std::endl;
        }

        if (cells_truncated_from_above_permx > 0) {
            std::cout << "Cells with permx truncated from above: "
                      << cells_truncated_from_above_permx << std::endl;
        }
    }
}

bool RelPermUpscaleHelper::checkCurve(MonotCubicInterpolator& func)
{
    double minrelp = func.getMinimumF().second;
    if (minrelp == 0)
        return true;
    else if (minrelp < critRelpThresh) {
        // set to 0
        std::vector<double> svec = func.get_xVector();
        std::vector<double> kvec = func.get_fVector();
        if (kvec[0] < critRelpThresh) {
            kvec[0] = 0.0;
        }
        else if (kvec[kvec.size()-1] < critRelpThresh) {
            kvec[kvec.size()-1] = 0.0;
        }
        func = MonotCubicInterpolator(svec, kvec);
        return true;
    }
    else
        return false;
}

void RelPermUpscaleHelper::checkCriticalSaturations()
{
    for (size_t i=0 ; i < Krfunctions[0][0].size(); ++i) {
        for (size_t j=0;j<(anisotropic_input?3:1);++j) {
            for (size_t k=0;k<(upscaleBothPhases?2:1);++k) {
                if (!checkCurve(Krfunctions[j][k][i])) {
                    std::stringstream str;
                    // Error message
                    str << "Relperm curve for rock " << i
                        << " does not specify critical saturation." << std::endl
                        << "Minimum relperm value is "
                        << Krfunctions[j][k][i].getMinimumF().second
                        << ", critRelpermThresh is " << critRelpThresh;
                    throw std::runtime_error(str.str());
                }
            }
        }
    }
}

void RelPermUpscaleHelper::setupBoundaryConditions()
{
    static const std::map<char, std::pair<SinglePhaseUpscaler::BoundaryConditionType, size_t>> bcmap =
        {{'f', {SinglePhaseUpscaler::Fixed,    3}},
         {'l', {SinglePhaseUpscaler::Linear,   9}},
         {'p', {SinglePhaseUpscaler::Periodic, 9}}};

    auto it = bcmap.find(options["bc"][0]);
    if (it != bcmap.end()) {
        boundaryCondition  = it->second.first;
        tensorElementCount = it->second.second;
    } else
        throw std::runtime_error("Invalid boundary condition. Only one of the letters f, l or p are allowed.");
}

double RelPermUpscaleHelper::tesselateGrid(const Opm::Deck& deck)
{
    const auto linsolver_tolerance         = to_double(options["linsolver_tolerance"]);
    const auto linsolver_verbosity         = to_int   (options["linsolver_verbosity"]);
    const auto linsolver_type              = to_int   (options["linsolver_type"]);
    const auto twodim_hack                 = false;
    const auto linsolver_maxit             = to_int   (options["linsolver_max_iterations"]);
    const auto smooth_steps                = to_int   (options["linsolver_smooth_steps"]);
    const auto linsolver_prolongate_factor = to_double(options["linsolver_prolongate_factor"]);
    const auto minPerm                     = to_double(options["minPerm"]);
    const auto gravity                     = to_double(options["gravity"]);

    if (isMaster) {
        std::cout << "Tesselating grid... ";
    }

    std::flush(std::cout);

    const auto start = clock();

    upscaler.init(deck, boundaryCondition,
                  unit::convert::from(minPerm, prefix::milli*unit::darcy),
                  linsolver_tolerance, linsolver_verbosity, linsolver_type,
                  twodim_hack, linsolver_maxit, linsolver_prolongate_factor,
                  smooth_steps, gravity);

    const auto finish = clock();
    const auto timeused_tesselation =
        (static_cast<double>(finish) - start) / CLOCKS_PER_SEC;

    if (isMaster) {
        std::cout << " (" << timeused_tesselation << " secs)" << std::endl;
    }

    return timeused_tesselation;
}

void RelPermUpscaleHelper::calculateCellPressureGradients()
{
    const auto& grid = upscaler.grid();

    auto celldepth = std::vector<double>(grid.numCells(), 0.0);
    auto mt        = ModelThickness{grid.logicalCartesianSize(), zcorns};
    {
        auto ijk = std::array<int,3>{{0, 0, 0}};

        for (auto c  = grid.leafbegin<0>(),
                 end = grid.leafend  <0>(); c != end; ++c)
        {
            const auto ix = c->index();

            grid.getIJK(ix, ijk);
            mt.activateCell(ijk);

            const auto& cc = c->geometry().center();
            celldepth[ix]  = cc[ cc.size() - 1 ];
        }
    }

    {
        const auto gravity      = to_double(options["gravity"]);
        const auto waterDensity = to_double(options["waterDensity"]);
        const auto oilDensity   = to_double(options["oilDensity"]);

        // Input water and oil density is given in g/cm3.
        const auto dRho = unit::convert::
            from(waterDensity - oilDensity,
                 prefix::milli*unit::kilogram /* 'g' missing from set */
                 / unit::cubic(prefix::centi * unit::meter));

        const auto thick = mt.thickness();

        dP.clear();  dP.reserve(celldepth.size());

        for (const auto& depth : celldepth) {
            dP.push_back(dRho * gravity * (depth - (thick / 2.0)));
        }
    }
}

void RelPermUpscaleHelper::calculateMinMaxCapillaryPressure()
{
    const auto maxPermContrast     = to_double(options["maxPermContrast"]);
    const auto minPerm             = to_double(options["minPerm"]);
    const auto gravity             = to_double(options["gravity"]);
    const auto linsolver_tolerance = to_double(options["linsolver_tolerance"]);
    const auto includeGravity      =
        (std::fabs(gravity) > std::numeric_limits<double>::min()); // true for non-zero gravity

    if (maxPermContrast == 0) {
        throw std::runtime_error("Illegal contrast value");
    }

    auto cellVolumes = std::vector<double>(satnums.size(), 0.0);
    cellPoreVolumes.resize(satnums.size(), 0.0);

    auto dPmin =  std::numeric_limits<double>::max();
    auto dPmax = -std::numeric_limits<double>::max();
    if (!dP.empty()) {
        const auto m = std::minmax_element(dP.begin(), dP.end());

        dPmax = *m.second;
        dPmin = *m.first;
    }

    // Find minimium and maximum capillary pressure values in each
    // cell, and use the global min/max as the two initial pressure
    // points for computations.
    //
    // Also find max single-phase permeability, used to obey the
    // maxPermContrast option.
    //
    // Also find properly upscaled saturation endpoints, these are
    // printed out to stdout for reference during computations, but will
    // automatically appear as the lowest and highest saturation points
    // in finished output.
    Pcmin =  std::numeric_limits<double>::max();
    Pcmax = -std::numeric_limits<double>::max();

    auto maxSinglePhasePerm = 0.0;
    auto Swirvolume         = 0.0;
    auto Sworvolume         = 0.0;

    const auto& ecl_idx = upscaler.grid().globalCell();

    for (auto c  = upscaler.grid().leafbegin<0>(),
             end = upscaler.grid().leafend<0>(); c != end; ++c)
    {
        const auto cell_idx = ecl_idx[c->index()];

        if (satnums[cell_idx] > 0) { // Satnum zero is "no rock"

            cellVolumes    [cell_idx] = c->geometry().volume();
            cellPoreVolumes[cell_idx] = cellVolumes[cell_idx] * poros[cell_idx];

            double Pcmincandidate, Pcmaxcandidate, minSw, maxSw;

            if (! anisotropic_input) {
                const auto& invJ  = InvJfunctions[satnums[cell_idx] - 1];
                const auto  kx    = perms[0][cell_idx];
                const auto  denom = std::sqrt(kx / poros[cell_idx]);

                Pcmincandidate = invJ.getMinimumX().first / denom;
                Pcmaxcandidate = invJ.getMaximumX().first / denom;

                minSw = invJ.getMinimumF().second;
                maxSw = invJ.getMaximumF().second;
            }
            else { // anisotropic input, we do not to J-function scaling
                const auto& satfunc = SwPcfunctions[satnums[cell_idx] - 1];

                Pcmincandidate = satfunc.getMinimumX().first;
                Pcmaxcandidate = satfunc.getMaximumX().first;

                minSw = satfunc.getMinimumF().second;
                maxSw = satfunc.getMaximumF().second;
            }

            Pcmin = std::min(Pcmincandidate, Pcmin);
            Pcmax = std::max(Pcmaxcandidate, Pcmax);

            maxSinglePhasePerm = std::max(maxSinglePhasePerm, perms[0][cell_idx]);

            // Add irreducible water saturation volume
            Swirvolume += minSw * cellPoreVolumes[cell_idx];
            Sworvolume += maxSw * cellPoreVolumes[cell_idx];
        }

        ++tesselatedCells; // keep count.
    }

    minSinglePhasePerm =
        std::max(maxSinglePhasePerm / maxPermContrast,
                 unit::convert::from(minPerm, prefix::milli*unit::darcy));

    if (includeGravity) {
        Pcmin -= dPmax;
        Pcmax -= dPmin;
    }

    if (isMaster) {
        std::cout << "Pcmin:    "
                  << unit::convert::to(Pcmin, unit::barsa)
                  << " [bar]\n"
                  << "Pcmax:    "
                  << unit::convert::to(Pcmax, unit::barsa)
                  << " [bar]" << std::endl;
    }

    if (Pcmin > Pcmax) {
        throw std::runtime_error("ERROR: No legal capillary pressures "
                                 "found for this system. Exiting...");
    }

    // Total porevolume and total volume -> upscaled porosity:
    poreVolume = std::accumulate(cellPoreVolumes.begin(), cellPoreVolumes.end(), 0.0);
    volume     = std::accumulate(cellVolumes.begin()    , cellVolumes.end()    , 0.0);

    Swir = Swirvolume / poreVolume;
    Swor = Sworvolume / poreVolume;

    if (isMaster) {
        std::cout << "LF Pore volume:    " << poreVolume << '\n'
                  << "LF Volume:         " << volume << '\n'
                  << "Upscaled porosity: " << (poreVolume / volume) << '\n'
                  << "Upscaled "           << saturationstring << "ir:     " << Swir << '\n'
                  << "Upscaled "           << saturationstring << "max:    " << Swor << '\n' //Swor=1-Swmax
                  << "Saturation points to be computed: "
                  << points << std::endl;
    }

    // Sometimes, if Swmax=1 or Swir=0 in the input tables, the upscaled
    // values can be a little bit larger (within machine precision) and the
    // check below fails. Hence, check if these values are within the [0,1]
    // interval within some precision (use linsolver_tolerance)
    if ((Swor > 1.0) && (Swor - linsolver_tolerance < 1.0)) {
        Swor = 1.0;
    }

    if ((Swir < 0.0) && (Swir + linsolver_tolerance > 0.0)) {
        Swir = 0.0;
    }

    if ((Swir < 0) || (Swir > 1) || (Swor < 0) || (Swor > 1)) {
        std::stringstream str;

        str << "ERROR: "        // e.g., Swir/Swor nonsensical...
            << saturationstring << "ir/"
            << saturationstring << "or "
            << "nonsensical. Check your input. Exiting";

        throw std::runtime_error(str.str());
    }
}

void RelPermUpscaleHelper::upscaleCapillaryPressure()
{
    const auto saturationThreshold =
        to_double(options["saturationThreshold"]);

    auto largestSaturationInterval = Swor - Swir;
    decltype(Pcmax) Ptestvalue;

    std::stringstream errstr;

    const auto& ecl_idx = upscaler.grid().globalCell();

    while (largestSaturationInterval > (Swor - Swir) / 500.0) {
        if (Pcmax == Pcmin) {
            // This is a dummy situation, we go through once and then we are
            // finished (this will be triggered by zero permeability)
            Ptestvalue = Pcmin;
            largestSaturationInterval = 0;
        }
        else if (WaterSaturationVsCapPressure.getSize() == 0) {
            // No data values previously computed
            Ptestvalue = Pcmax;
        }
        else if (WaterSaturationVsCapPressure.getSize() == 1) {
            // If only one point has been computed, it was for Pcmax. So now
            // do Pcmin.
            Ptestvalue = Pcmin;
        }
        else {
            // Search for largest saturation interval in which there are no
            // computed saturation points (and estimate the capillary
            // pressure that will fall in the center of this saturation
            // interval)
            const auto& SatDiff = WaterSaturationVsCapPressure.getMissingX();

            Ptestvalue                = SatDiff.first;
            largestSaturationInterval = SatDiff.second;
        }

        // Check for saneness of Ptestvalue:
        if (std::isnan(Ptestvalue) || std::isinf(Ptestvalue)) {
            errstr << "ERROR: Ptestvalue was inf or nan\n";
            break; // Jump out of while-loop, just print out the results
                   // up to now and exit the program
        }

        auto waterVolume = 0.0;
        for (decltype(ecl_idx.size())
                 i = 0, n = ecl_idx.size(); i < n; ++i)
        {
            const auto cell_idx = ecl_idx[i];

            double waterSaturationCell = 0.0;

            if (satnums[cell_idx] > 0) { // handle "no rock" cells with satnum zero
                const auto ix = satnums[cell_idx] - 1;

                auto PtestvalueCell = Ptestvalue;
                if (!dP.empty()) {
                    PtestvalueCell -= dP[cell_idx];
                }

                if (!anisotropic_input) {
                    const auto arg = perms[0][cell_idx] / poros[cell_idx];

                    const auto Jvalue = std::sqrt(arg) * PtestvalueCell;

                    waterSaturationCell = InvJfunctions[ix].evaluate(Jvalue);
                }
                else {
                    // anisotropic_input, then we do not do J-function-scaling
                    waterSaturationCell = SwPcfunctions[ix].evaluate(PtestvalueCell);
                }
            }

            waterVolume += waterSaturationCell * cellPoreVolumes[cell_idx];
        }

        WaterSaturationVsCapPressure.addPair(Ptestvalue, waterVolume / poreVolume);
    }

    // Now, it may happen that we have a large number of cells, and
    // some cells with near zero poro and perm. This may cause that
    // Pcmax has been estimated so high that it does not affect Sw
    // within machine precision, and then we need to truncate the
    // largest Pc values:
    WaterSaturationVsCapPressure.chopFlatEndpoints(saturationThreshold);

    // Now we can also invert the upscaled water saturation
    // (it should be monotonic)
    if (! WaterSaturationVsCapPressure.isStrictlyMonotone()) {
        errstr << "Error: Upscaled water saturation not strictly "
               << "monotone in capillary pressure.\n"
               << "       Unphysical input data, exiting.\n"
               << "       Trying to dump " << saturationstring
               << " vs Pc to file swvspc_debug.txt for inspection";

        if (isMaster) {
            std::ofstream outfile("swvspc_debug.txt");

            outfile << "# Pc      " << saturationstring << std::endl;
            outfile << WaterSaturationVsCapPressure.toString();
        }

        throw std::runtime_error(errstr.str());
    }
}

std::tuple<double, double>
RelPermUpscaleHelper::upscalePermeability(int mpi_rank)
{
    const auto minPerm         = to_double(options["minPerm"]);
    const auto maxPermContrast = to_double(options["maxPermContrast"]);

     // Put correct number of zeros in, just to be able to access RelPerm[index] later
    WaterSaturation.resize(points);
    for (size_t i = 0; i < (upscaleBothPhases ? 2 : 1); ++i) {
        PhasePerm[i].resize(points, std::vector<double>(tensorElementCount));
    }

    // Make vector of capillary pressure points corresponding to uniformly
    // distribued saturation points between Swor and Swir.
    {
        auto CapPressureVsWaterSaturation =
            MonotCubicInterpolator(WaterSaturationVsCapPressure.get_fVector(),
                                   WaterSaturationVsCapPressure.get_xVector());

        for (int pointidx = 1; pointidx <= points; ++pointidx) {
            // pointidx=1 corresponds to Swir, pointidx=points to Swor.
            double saturation = Swir + (Swor-Swir)/(points-1)*(pointidx-1);
            pressurePoints.push_back(CapPressureVsWaterSaturation.evaluate(saturation));
        }
    }

    // Preserve max and min pressures
    pressurePoints.front() = Pcmax;
    pressurePoints.back()  = Pcmin;

    // Fill with zeros initially (in case of non-mpi)
    node_vs_pressurepoint.resize(points);

#if defined(HAVE_MPI) && HAVE_MPI
    // Distribute work load over mpi nodes.
    int mpi_nodecount;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_nodecount);

    for (int idx = 0; idx < points; ++idx) {
        // Ensure master node gets equal or less work than the other nodes,
        // since master node also computes single phase perm.
        node_vs_pressurepoint[idx] = (mpi_nodecount-1) - idx % mpi_nodecount;
    }
#endif

    const auto& ecl_idx = upscaler.grid().globalCell();

    const auto start_upscale_wallclock = clock();

    double waterVolumeLF = 0.0;

    // Now loop through the vector of capillary pressure points that this
    // node should compute.
    for (int pointidx = 0; pointidx < points; ++pointidx) {

        // Should "I" (mpi-wise) compute this pressure point?
        if (node_vs_pressurepoint[pointidx] == mpi_rank) {

            const auto Ptestvalue = pressurePoints[pointidx];

            std::array<double,2> maxPhasePerm{{0.0, 0.0}};
            std::array<std::vector<double>,2> phasePermValues;
            std::array<std::vector<std::vector<double>>,2> phasePermValuesDiag;
            std::array<double,2> minPhasePerm;
            std::array<SinglePhaseUpscaler::permtensor_t,2> phasePermTensor;

            for (size_t p = 0; p < (upscaleBothPhases ? 2 : 1); ++p) {
                waterVolumeLF = 0.0;
                phasePermValues    [p].resize(satnums.size());
                phasePermValuesDiag[p].resize(satnums.size());

                for (decltype(ecl_idx.size())
                         i = 0, n = ecl_idx.size(); i < n; ++i)
                {
                    const auto cell_idx = ecl_idx[i];
                    double cellPhasePerm =
                        unit::convert::from(minPerm, prefix::milli*unit::darcy);

                    auto cellPhasePermDiag =
                        std::vector<double>(3, cellPhasePerm);

                    const auto ix = satnums[cell_idx] - 1;
                    const auto kx = perms[0][cell_idx];

                    if (satnums[cell_idx] > 0) { // handle "no rock" cells with satnum zero
                        auto PtestvalueCell = Ptestvalue;
                        if (!dP.empty()) {
                            PtestvalueCell -= dP[cell_idx];
                        }

                        if (!anisotropic_input) {
                            const auto Jvalue =
                                std::sqrt(kx / poros[cell_idx]) * PtestvalueCell;

                            const auto WaterSaturationCell =
                                InvJfunctions[ix].evaluate(Jvalue);

                            waterVolumeLF += WaterSaturationCell * cellPoreVolumes[cell_idx];

                            // Compute cell relative permeability. We use a lower cutoff-value as we
                            // easily divide by zero here.  When water saturation is
                            // zero, we get 'inf', which is circumvented by the cutoff value.
                            cellPhasePerm =
                                Krfunctions[0][p][ix].evaluate(WaterSaturationCell) * kx;
                        }
                        else {
                            const auto WaterSaturationCell =
                                SwPcfunctions[ix].evaluate(PtestvalueCell);

                            waterVolumeLF += WaterSaturationCell * cellPoreVolumes[cell_idx];

                            cellPhasePermDiag[0] =
                                Krfunctions[0][p][ix].evaluate(WaterSaturationCell) * kx;

                            cellPhasePermDiag[1] =
                                Krfunctions[1][p][ix].evaluate(WaterSaturationCell) * perms[1][cell_idx];

                            cellPhasePermDiag[2] =
                                Krfunctions[2][p][ix].evaluate(WaterSaturationCell) * perms[2][cell_idx];
                        }

                        phasePermValues    [p][cell_idx] = cellPhasePerm;
                        phasePermValuesDiag[p][cell_idx] = cellPhasePermDiag;

                        maxPhasePerm[p] = std::max(maxPhasePerm[p], cellPhasePerm);
                        maxPhasePerm[p] = std::max(maxPhasePerm[p],
                                                   *std::max_element(cellPhasePermDiag.begin(),
                                                                     cellPhasePermDiag.end()));
                    }
                }

                // Now we can determine the smallest permitted permeability
                // we can calculate for We have both a fixed bottom limit,
                // as well as a possible higher limit determined by a
                // maximum allowable permeability.
                minPhasePerm[p] = std::max(maxPhasePerm[p] / maxPermContrast,
                                           unit::convert::from(minPerm, prefix::milli*unit::darcy));

                // Now remodel the phase permeabilities obeying minPhasePerm
                SinglePhaseUpscaler::permtensor_t cellperm(3, 3, nullptr);
                for (decltype(ecl_idx.size())
                         i = 0, n = ecl_idx.size(); i < n; ++i)
                {
                    const auto cell_idx = ecl_idx[i];
                    zero(cellperm);

                    if (!anisotropic_input) {
                        const auto cellPhasePerm =
                            std::max(minPhasePerm[p], phasePermValues[p][cell_idx]);

                        const auto kval = std::max(minPhasePerm[p], cellPhasePerm);

                        cellperm(0,0) = kval;
                        cellperm(1,1) = kval;
                        cellperm(2,2) = kval;
                    }
                    else { // anisotropic_input
                        // Truncate values lower than minPhasePerm upwards.
                        auto& k = phasePermValuesDiag[p][cell_idx];

                        cellperm(0,0) = k[0] = std::max(minPhasePerm[p], k[0]);
                        cellperm(1,1) = k[1] = std::max(minPhasePerm[p], k[1]);
                        cellperm(2,2) = k[2] = std::max(minPhasePerm[p], k[2]);
                    }

                    upscaler.setPermeability(i, cellperm);
                }

                //  Call single-phase upscaling code
                phasePermTensor[p] = upscaler.upscaleSinglePhase();
            }

            // Here we recalculate the upscaled water saturation,
            // although it is already known when we asked for the
            // pressure point to compute for. Nonetheless, we
            // recalculate here to avoid any minor roundoff-error and
            // interpolation error (this means that the saturation
            // points are not perfectly uniformly distributed)
            WaterSaturation[pointidx] =  waterVolumeLF/poreVolume;

#if defined(HAVE_MPI) && HAVE_MPI
            std::cout << "Rank " << mpi_rank << ": ";
#endif  // HAVE_MPI

            std::cout << Ptestvalue << "\t" << WaterSaturation[pointidx];

            // Store and print phase-perm-result
            for (int voigtIdx=0; voigtIdx < tensorElementCount; ++voigtIdx) {
                for (size_t p = 0; p < (upscaleBothPhases ? 2 : 1); ++p) {
                    PhasePerm[p][pointidx][voigtIdx] =
                        getVoigtValue(phasePermTensor[p], voigtIdx);

                    std::cout << "\t" << PhasePerm[p][pointidx][voigtIdx];
                }
            }

            std::cout << '\n';
        }
    }

    clock_t finish_upscale_wallclock = clock();
    double timeused_upscale_wallclock =
        (double(finish_upscale_wallclock) -
         double(start_upscale_wallclock)) / CLOCKS_PER_SEC;

    collectResults();

    // Average time pr. upscaling point:
#if defined(HAVE_MPI) && HAVE_MPI

    // Sum the upscaling time used by all processes
    double timeused_total;
    MPI_Reduce(&timeused_upscale_wallclock, &timeused_total, 1,
               MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    const auto avg_upscaling_time_pr_point = timeused_total / points;

#else

    const auto avg_upscaling_time_pr_point =
        timeused_upscale_wallclock / points;

#endif  // HAVE_MPI

    return std::make_tuple(timeused_upscale_wallclock,
                           avg_upscaling_time_pr_point);
}

}
