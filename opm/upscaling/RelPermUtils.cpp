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
#include <algorithm>
#include <iostream>
#include <exception>
#include <sstream>

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
                                           std::map<std::string,std::string>& options) :
    isMaster(mpi_rank == 0),
    anisotropic_input(false),
    permTensor(3,3,nullptr),
    permTensorInv(3,3,nullptr),
    tesselatedCells(0)
{
    if (options["fluids"] == "ow" || options["fluids"] == "wo")
        saturationstring = "Sw";
    else if (options["fluids"] == "go" || options["fluids"] == "og")
        saturationstring = "Sg";
    else {
        std::stringstream str;
        str << "Fluidsystem " << options["fluids"] << " not valid (-fluids option). Should be ow or go";
        throw std::runtime_error(str.str());
    }
}

void RelPermUpscaleHelper::collectResults()
{
#ifdef HAVE_MPI
   /* Step 8b: Transfer all computed data to master node.
      Master node should post a receive for all values missing,
      other nodes should post a send for all the values they have.
    */
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
                       PhasePerm[(int)recvbuffer[0]][voigtIdx] = recvbuffer[2+voigtIdx];
                   }
                   for (int voigtIdx=0; voigtIdx < tensorElementCount; ++voigtIdx) {
                       Phase2Perm[(int)recvbuffer[0]][voigtIdx] = recvbuffer[2+tensorElementCount+voigtIdx];
                   }
               }
               else {
                   std::vector<double> recvbuffer(2+tensorElementCount);
                   MPI_Recv(recvbuffer.data(), recvbuffer.size(), MPI_DOUBLE,
                            node_vs_pressurepoint[idx], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                   // Put received data into correct place.
                   WaterSaturation[(int)recvbuffer[0]] = recvbuffer[1];
                   for (int voigtIdx=0; voigtIdx < tensorElementCount; ++voigtIdx) {
                       PhasePerm[(int)recvbuffer[0]][voigtIdx] = recvbuffer[2+voigtIdx];
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
                       sendbuffer[2+voigtIdx] = PhasePerm[idx][voigtIdx];
                   }
                   for (int voigtIdx=0; voigtIdx < tensorElementCount; ++voigtIdx) {
                       sendbuffer[2+tensorElementCount+voigtIdx] = Phase2Perm[idx][voigtIdx];
                   }
                   MPI_Send(sendbuffer.data(), sendbuffer.size(), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
               }
               else {
                   std::vector<double> sendbuffer(2+tensorElementCount);
                   sendbuffer[0] = (double)idx;
                   sendbuffer[1] = WaterSaturation[idx];
                   for (int voigtIdx=0; voigtIdx < tensorElementCount; ++voigtIdx) {
                       sendbuffer[2+voigtIdx] = PhasePerm[idx][voigtIdx];
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
        SinglePhaseUpscaler::permtensor_t cellperm(3,3,nullptr);
        zero(cellperm);
        const std::vector<int>& ecl_idx = upscaler.grid().globalCell();
        for (unsigned int i = 0; i < ecl_idx.size(); ++i) {
            unsigned int cell_idx = ecl_idx[i];
            zero(cellperm);
            if (! anisotropic_input) {
                double kval = std::max(perms[0][cell_idx], minSinglePhasePerm);
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
        permTensor = upscaler.upscaleSinglePhase();
        permTensorInv = permTensor;
        invert(permTensorInv);
    }
}

void RelPermUpscaleHelper::sanityCheckInput(Opm::DeckConstPtr deck,
                                            double minPerm,
                                            double maxPerm,
                                            double minPoro)

{
    // Check that we have the information we need from the eclipse file:
    if (! (deck->hasKeyword("SPECGRID") && deck->hasKeyword("COORD") && deck->hasKeyword("ZCORN")
           && deck->hasKeyword("PORO") && deck->hasKeyword("PERMX"))) {
      throw std::runtime_error("Error: Did not find SPECGRID, COORD, ZCORN, PORO and PERMX in Eclipse file.");
    }

    poros    = deck->getKeyword("PORO")->getRawDoubleData();
    perms[0] = deck->getKeyword("PERMX")->getRawDoubleData();
    zcorns   = deck->getKeyword("ZCORN")->getRawDoubleData();

    // Load anisotropic (only diagonal supported) input if present in grid
    if (deck->hasKeyword("PERMY") && deck->hasKeyword("PERMZ")) {
        anisotropic_input = true;
        perms[1] = deck->getKeyword("PERMY")->getRawDoubleData();
        perms[2] = deck->getKeyword("PERMZ")->getRawDoubleData();
        if (isMaster)
          std::cout << "Info: PERMY and PERMZ present, going into anisotropic input mode, no J-functions\n"
                    << "      Options -relPermCurve and -jFunctionCurve is meaningless.\n";
    }


    /* Initialize a default satnums-vector with only "ones" (meaning only one rocktype) */
    satnums.resize(poros.size(), 1);

    if (deck->hasKeyword("SATNUM")) {
        satnums = deck->getKeyword("SATNUM")->getIntData();
    }
    else {
        if (isMaster)
          std::cout << "SATNUM not found in input file, assuming only one rocktype" << std::endl;
    }

    /* Sanity check/fix on input for each cell:
       - Check that SATNUM are set sensibly, that is => 0 and < 1000, error if not.
       - Check that porosity is between 0 and 1, error if not.
       Set to minPoro if zero or less than minPoro (due to pcmin/max computation)
       - Check that permeability is zero or positive. Error if negative.
       Set to minPerm if zero or less than minPerm.
       - Check maximum number of SATNUM values (can be number of rock types present)
    */
    int cells_truncated_from_below_poro = 0;
    int cells_truncated_from_below_permx = 0;
    int cells_truncated_from_above_permx = 0;
    auto it = std::find_if(satnums.begin(), satnums.end(), [](int a) { return a < 0; });
    auto it2 = std::find_if(satnums.begin(), satnums.end(), [](int a) { return a > 1000;});
    if (it != satnums.end() || it2 != satnums.end()) {
        std::stringstream str;
        str << "satnums[";
        if (it == satnums.end())
            str << it2-satnums.begin() << "] = " << *it2;
        else
            str << it-satnums.begin() << "] = " << *it;
        str << ", not sane, quitting.";
        throw std::runtime_error(str.str());
    }

    auto&& find_error_both  = [](double value) { return value < 0 || value > 1; };
    auto it3 = std::find_if(poros.begin(), poros.end(), find_error_both);
    if (it3 != poros.end()) {
        std::stringstream str;
        str << "poros[" << it3-poros.begin() <<"] = " << *it << ", not sane, quitting.";
        throw std::runtime_error(str.str());
    }
    auto&& find_error_below = [](double value) { return value < 0; };
    std::string name = "permx";
    auto&& check_perm = [name,&find_error_below](const std::vector<double>& perm)
                        {
                            auto it = std::find_if(perm.begin(), perm.end(), find_error_below);
                            if (it != perm.end()) {
                                std::stringstream str;
                                str << name <<"[" << it-perm.begin() <<"] = " << *it << ", not sane, quitting.";
                                throw std::runtime_error(str.str());
                            }
                        };
    check_perm(perms[0]);
    if (anisotropic_input) {
        name ="permy";
        check_perm(perms[1]);
        name ="permz";
        check_perm(perms[2]);
    }
    std::transform(poros.begin(), poros.end(), poros.begin(),
                   [minPoro,&cells_truncated_from_below_poro](double value)
                   {
                       if (value >= 0 && value < minPoro) { // Truncate porosity from below
                           ++cells_truncated_from_below_poro;
                           return minPoro;
                       }
                       return value;
                   });

    std::transform(perms[0].begin(), perms[0].end(), perms[0].begin(),
                   [minPerm, maxPerm, &cells_truncated_from_below_permx,
                    &cells_truncated_from_above_permx](double value)
                   {
                       if ((value >= 0) && (value < minPerm)) { // Truncate permeability from below
                           ++cells_truncated_from_below_permx;
                           return minPerm;
                       }
                       if (value > maxPerm) { // Truncate permeability from above
                           ++cells_truncated_from_above_permx;
                           return maxPerm;
                       }

                       return value;
                   });

    for (unsigned int i = 0; i < satnums.size(); ++i) {
        // Explicitly handle "no rock" cells, set them to minimum perm and zero porosity.
        if (satnums[i] == 0) {
            perms[0][i] = minPerm;
            if (anisotropic_input) {
                perms[1][i] = minPerm;
                perms[2][i] = minPerm;
            }
            poros[i] = 0; // zero poro is fine for these cells, as they are not
            // used in pcmin/max computation.
        }
    }
    if (isMaster && cells_truncated_from_below_poro > 0)
        std::cout << "Cells with truncated porosity: " << cells_truncated_from_below_poro << std::endl;
    if (isMaster && cells_truncated_from_below_permx > 0)
        std::cout << "Cells with permx truncated from below: " << cells_truncated_from_below_permx << std::endl;
    if (isMaster && cells_truncated_from_above_permx > 0)
        std::cout << "Cells with permx truncated from above: " << cells_truncated_from_above_permx << std::endl;
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

void RelPermUpscaleHelper::setupBoundaryConditions(std::map<std::string,std::string>& options)
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

double RelPermUpscaleHelper::tesselateGrid(Opm::DeckConstPtr deck,
                                         std::map<std::string,std::string>& options)
{
  double linsolver_tolerance = atof(options["linsolver_tolerance"].c_str());
  int linsolver_verbosity = atoi(options["linsolver_verbosity"].c_str());
  int linsolver_type = atoi(options["linsolver_type"].c_str());
  int linsolver_maxit = atoi(options["linsolver_max_iterations"].c_str());
  int smooth_steps = atoi(options["linsolver_smooth_steps"].c_str());
  double linsolver_prolongate_factor = atof(options["linsolver_prolongate_factor"].c_str());
  const double minPerm = atof(options["minPerm"].c_str());

   if (isMaster)
     std::cout << "Tesselating grid... ";

   flush(std::cout);
   clock_t start = clock();

  upscaler.init(deck, boundaryCondition,
                Opm::unit::convert::from(minPerm, Opm::prefix::milli*Opm::unit::darcy),
                linsolver_tolerance, linsolver_verbosity, linsolver_type, false,
                linsolver_maxit, linsolver_prolongate_factor, smooth_steps);

   clock_t finish = clock();
   double timeused_tesselation = (double(finish)-double(start))/CLOCKS_PER_SEC;
   if (isMaster)
     std::cout << " (" << timeused_tesselation <<" secs)" << std::endl;

   return timeused_tesselation;
}

std::vector<double>
        RelPermUpscaleHelper::calculateCellPressureGradients(const std::array<int,3>& res,
                                                             std::map<std::string,std::string>& options)
{
    const double gravity            = atof(options["gravity"].c_str());
    const double waterDensity       = atof(options["waterDensity"].c_str());
    const double oilDensity         = atof(options["oilDensity"].c_str());

    // height of model is calculated as the average of the z-values at the top layer
    // This calculation makes assumption on the indexing of cells in the grid, going from bottom to top.
    double modelHeight = 0;
    for (size_t zIdx = (4 * res[0] * res[1] * (2*res[2]-1)); zIdx < zcorns.size(); ++zIdx)
        modelHeight += zcorns[zIdx] / (4*res[0]*res[1]);

    // We assume that the spatial units in the grid file is in centimetres,
    // so we divide by 100 to get to metres.
    modelHeight /= 100.0;

    // Input water and oil density is given in g/cm3, we convert it to kg/m3 (SI)
    // by multiplying with 1000.
    double dRho = (waterDensity-oilDensity) * 1000; // SI unit (kg/m3)

    // Calculating difference in capillary pressure for all cells
    std::vector<double> dP(satnums.size(), 0);
    for (size_t cellIdx = 0; cellIdx < satnums.size(); ++cellIdx) {
        int i,j,k; // Position of cell in cell hierarchy
        std::vector<int> zIndices(8,0); // 8 corners with 8 heights
        int horIdx = (cellIdx+1) - int(std::floor(((double)(cellIdx+1))/((double)(res[0]*res[1]))))*res[0]*res[1]; // index in the corresponding horizon
        if (horIdx == 0)
            horIdx = res[0]*res[1];

        i = horIdx - int(std::floor(((double)horIdx)/((double)res[0])))*res[0];

        if (i == 0)
            i = res[0];

        j = (horIdx-i)/res[0]+1;
        k = ((cellIdx+1)-res[0]*(j-1)-1)/(res[0]*res[1])+1;
        int zBegin = 8*res[0]*res[1]*(k-1); // indices of Z-values of bottom
        int level2 = 4*res[0]*res[1]; // number of z-values in one horizon
        zIndices[0] = zBegin + 4*res[0]*(j-1)+2*i-1;
        zIndices[1] = zBegin + 4*res[0]*(j-1)+2*i;
        zIndices[2] = zBegin + 2*res[0]*(2*j-1)+2*i;
        zIndices[3] = zBegin + 2*res[0]*(2*j-1)+2*i-1;
        zIndices[4] = zBegin + level2 + 4*res[0]*(j-1)+2*i-1;
        zIndices[5] = zBegin + level2 + 4*res[0]*(j-1)+2*i;
        zIndices[6] = zBegin + level2 + 2*res[0]*(2*j-1)+2*i;
        zIndices[7] = zBegin + level2 + 2*res[0]*(2*j-1)+2*i-1;

        double cellDepth = 0;
        for (size_t corner = 0; corner < 8; ++corner)
            cellDepth += zcorns[zIndices[corner]-1] / 8.0;

        // cellDepth is in cm, convert to m by dividing by 100
        cellDepth /= 100.0;
        dP[cellIdx] = dRho * gravity * (cellDepth-modelHeight/2.0);
    }

    return dP;
}

void RelPermUpscaleHelper::calculateMinMaxCapillaryPressure(double dPmin, double dPmax,
                                                            std::map<std::string,std::string>& options)
{
    const double maxPermContrast = atof(options["maxPermContrast"].c_str());
    const double minPerm         = atof(options["minPerm"].c_str());
    const double gravity         = atof(options["gravity"].c_str());
    double linsolver_tolerance   = atof(options["linsolver_tolerance"].c_str());
    const bool includeGravity    = (fabs(gravity) > DBL_MIN); // true for non-zero gravity
    const double milliDarcyToSqMetre =
                        Opm::unit::convert::to(1.0*Opm::prefix::milli*Opm::unit::darcy,
                                               Opm::unit::square(Opm::unit::meter));
    if (maxPermContrast == 0)
        throw std::runtime_error("Illegal contrast value");

    std::vector<double> cellVolumes(satnums.size(), 0.0);
    cellPoreVolumes.resize(satnums.size(), 0.0);

    /* Find minimium and maximum capillary pressure values in each
       cell, and use the global min/max as the two initial pressure
       points for computations.

       Also find max single-phase permeability, used to obey the
       maxPermContrast option.

       Also find properly upscaled saturation endpoints, these are
       printed out to stdout for reference during computations, but will
       automatically appear as the lowest and highest saturation points
       in finished output.
    */
    Pcmax = -DBL_MAX, Pcmin = DBL_MAX;
    double maxSinglePhasePerm = 0;
    double Swirvolume = 0;
    double Sworvolume = 0;
    // cell_idx is the eclipse index.
    const std::vector<int>& ecl_idx = upscaler.grid().globalCell();
    Dune::CpGrid::Codim<0>::LeafIterator c = upscaler.grid().leafbegin<0>();
    for (; c != upscaler.grid().leafend<0>(); ++c) {
        unsigned int cell_idx = ecl_idx[c->index()];
        if (satnums[cell_idx] > 0) { // Satnum zero is "no rock"

            cellVolumes[cell_idx] = c->geometry().volume();
            cellPoreVolumes[cell_idx] = cellVolumes[cell_idx] * poros[cell_idx];

            double Pcmincandidate, Pcmaxcandidate, minSw, maxSw;

            if (! anisotropic_input) {
                Pcmincandidate = InvJfunctions[int(satnums[cell_idx])-1].getMinimumX().first
                    / sqrt(perms[0][cell_idx] * milliDarcyToSqMetre / poros[cell_idx]);
                Pcmaxcandidate = InvJfunctions[int(satnums[cell_idx])-1].getMaximumX().first
                    / sqrt(perms[0][cell_idx] * milliDarcyToSqMetre/poros[cell_idx]);
                minSw = InvJfunctions[int(satnums[cell_idx])-1].getMinimumF().second;
                maxSw = InvJfunctions[int(satnums[cell_idx])-1].getMaximumF().second;
            }
            else { // anisotropic input, we do not to J-function scaling
                Pcmincandidate = SwPcfunctions[int(satnums[cell_idx])-1].getMinimumX().first;
                Pcmaxcandidate = SwPcfunctions[int(satnums[cell_idx])-1].getMaximumX().first;

                minSw = SwPcfunctions[int(satnums[cell_idx])-1].getMinimumF().second;
                maxSw = SwPcfunctions[int(satnums[cell_idx])-1].getMaximumF().second;
            }
            Pcmin = std::min(Pcmincandidate, Pcmin);
            Pcmax = std::max(Pcmaxcandidate, Pcmax);

            maxSinglePhasePerm = std::max( maxSinglePhasePerm, perms[0][cell_idx]);

            // Add irreducible water saturation volume
            Swirvolume += minSw * cellPoreVolumes[cell_idx];
            Sworvolume += maxSw * cellPoreVolumes[cell_idx];
        }
        ++tesselatedCells; // keep count.
    }
    minSinglePhasePerm = std::max(maxSinglePhasePerm/maxPermContrast, minPerm);

    if (includeGravity) {
        Pcmin -= dPmax;
        Pcmax -= dPmin;
    }

    if (isMaster) {
        std::cout << "Pcmin:    " << Pcmin << std::endl;
        std::cout << "Pcmax:    " << Pcmax << std::endl;
    }

    if (Pcmin > Pcmax)
        throw std::runtime_error("ERROR: No legal capillary pressures found for this system. Exiting...");

    // Total porevolume and total volume -> upscaled porosity:
    poreVolume = std::accumulate(cellPoreVolumes.begin(), cellPoreVolumes.end(), 0.0);
    volume = std::accumulate(cellVolumes.begin(), cellVolumes.end(), 0.0);

    Swir = Swirvolume/poreVolume;
    Swor = Sworvolume/poreVolume;

    if (isMaster) {
        std::cout << "LF Pore volume:    " << poreVolume << std::endl;
        std::cout << "LF Volume:         " << volume << std::endl;
        std::cout << "Upscaled porosity: " << poreVolume/volume << std::endl;
        std::cout << "Upscaled " << saturationstring << "ir:     " << Swir << std::endl;
        std::cout << "Upscaled " << saturationstring << "max:    " << Swor << std::endl; //Swor=1-Swmax
        std::cout << "Saturation points to be computed: " << points << std::endl;
    }

    // Sometimes, if Swmax=1 or Swir=0 in the input tables, the upscaled
    // values can be a little bit larger (within machine precision) and
    // the check below fails. Hence, check if these values are within the
    // the [0 1] interval within some precision (use linsolver_precision)
    if (Swor > 1.0 && Swor - linsolver_tolerance < 1.0) {
        Swor = 1.0;
    }
    if (Swir < 0.0 && Swir + linsolver_tolerance > 0.0) {
        Swir = 0.0;
    }
    if (Swir < 0 || Swir > 1 || Swor < 0 || Swor > 1) {
        std::stringstream str;
        str << "ERROR: " << saturationstring << "ir/" << saturationstring << "or unsensible. Check your input. Exiting";
        throw std::runtime_error(str.str());
    }
}

void RelPermUpscaleHelper::upscaleCapillaryPressure(std::map<std::string,std::string>& options,
                                                    const std::vector<double>& dP)
{
    const double saturationThreshold = atof(options["saturationThreshold"].c_str());
    double largestSaturationInterval = Swor-Swir;
    double Ptestvalue = Pcmax;
    std::stringstream errstr;
    const std::vector<int>& ecl_idx = upscaler.grid().globalCell();
    const double milliDarcyToSqMetre =
                        Opm::unit::convert::to(1.0*Opm::prefix::milli*Opm::unit::darcy,
                                               Opm::unit::square(Opm::unit::meter));

    while (largestSaturationInterval > (Swor-Swir)/500.0) {
        if (Pcmax == Pcmin) {
            // This is a dummy situation, we go through once and then
            // we are finished (this will be triggered by zero permeability)
            Ptestvalue = Pcmin;
            largestSaturationInterval = 0;
        }
        else if (WaterSaturationVsCapPressure.getSize() == 0) {
            /* No data values previously computed */
            Ptestvalue = Pcmax;
        }
        else if (WaterSaturationVsCapPressure.getSize() == 1) {
            /* If only one point has been computed, it was for Pcmax. So now
               do Pcmin */
            Ptestvalue = Pcmin;
        }
        else {
            /* Search for largest saturation interval in which there are no
               computed saturation points (and estimate the capillary pressure
               that will fall in the center of this saturation interval)
               */
            std::pair<double,double> SatDiff = WaterSaturationVsCapPressure.getMissingX();
            Ptestvalue = SatDiff.first;
            largestSaturationInterval = SatDiff.second;
        }

        // Check for saneness of Ptestvalue:
        if (std::isnan(Ptestvalue) || std::isinf(Ptestvalue)) {
            errstr << "ERROR: Ptestvalue was inf or nan" << std::endl;
            break; // Jump out of while-loop, just print out the results
                   // up to now and exit the program
        }

        double waterVolume = 0.0;
        for (size_t i = 0; i < ecl_idx.size(); ++i) {
            unsigned int cell_idx = ecl_idx[i];
            double waterSaturationCell = 0.0;
            if (satnums[cell_idx] > 0) { // handle "no rock" cells with satnum zero
                double PtestvalueCell = Ptestvalue;
                if (!dP.empty())
                    PtestvalueCell -= dP[cell_idx];

                if (!anisotropic_input) {
                    double Jvalue = sqrt(perms[0][cell_idx] * milliDarcyToSqMetre / poros[cell_idx]) * PtestvalueCell;
                    waterSaturationCell = InvJfunctions[int(satnums[cell_idx])-1].evaluate(Jvalue);
                }
                else // anisotropic_input, then we do not do J-function-scaling
                    waterSaturationCell = SwPcfunctions[int(satnums[cell_idx])-1].evaluate(PtestvalueCell);
            }
            waterVolume += waterSaturationCell  * cellPoreVolumes[cell_idx];
        }
        WaterSaturationVsCapPressure.addPair(Ptestvalue, waterVolume/poreVolume);
    }

    // Now, it may happen that we have a large number of cells, and
    // some cells with near zero poro and perm. This may cause that
    // Pcmax has been estimated so high that it does not affect Sw
    // within machine precision, and then we need to truncate the
    // largest Pc values:
    WaterSaturationVsCapPressure.chopFlatEndpoints(saturationThreshold);

    // Now we can also invert the upscaled water saturation
    // (it should be monotonic)
    if (!WaterSaturationVsCapPressure.isStrictlyMonotone()) {
        errstr <<  "Error: Upscaled water saturation not strictly monotone in capillary pressure." << std::endl
               << "       Unphysical input data, exiting." << std::endl
               << "       Trying to dump " << saturationstring << " vs Pc to file swvspc_debug.txt for inspection";
        if (isMaster) {
            std::ofstream outfile;
            outfile.open("swvspc_debug.txt", std::ios::out | std::ios::trunc);
            outfile << "# Pc      " << saturationstring << std::endl;
            outfile << WaterSaturationVsCapPressure.toString();
            outfile.close();
        }
        throw std::runtime_error(errstr.str());
    }
}

std::tuple<double, double>
    RelPermUpscaleHelper::upscalePermeability(std::map<std::string,std::string>& options,
                                              const std::vector<double>& dP,
                                              int mpi_rank)
{
    const double minPerm = atof(options["minPerm"].c_str());
    const double maxPermContrast = atof(options["maxPermContrast"].c_str());
    const double milliDarcyToSqMetre =
                        Opm::unit::convert::to(1.0*Opm::prefix::milli*Opm::unit::darcy,
                                               Opm::unit::square(Opm::unit::meter));

     // Put correct number of zeros in, just to be able to access RelPerm[index] later
    WaterSaturation.resize(points);
    for (size_t i = 0; i < (upscaleBothPhases?2:1); ++i)
        PhasePerm[i].resize(points, std::vector<double>(tensorElementCount));

    // Make vector of capillary pressure points corresponding to uniformly distribued
    // saturation points between Swor and Swir.
    MonotCubicInterpolator CapPressureVsWaterSaturation(WaterSaturationVsCapPressure.get_fVector(),
                                                        WaterSaturationVsCapPressure.get_xVector());

    for (int pointidx = 1; pointidx <= points; ++pointidx) {
        // pointidx=1 corresponds to Swir, pointidx=points to Swor.
        double saturation = Swir + (Swor-Swir)/(points-1)*(pointidx-1);
        pressurePoints.push_back(CapPressureVsWaterSaturation.evaluate(saturation));
    }
    // Preserve max and min pressures
    pressurePoints.front() = Pcmax;
    pressurePoints.back()  = Pcmin;

    // Fill with zeros initially (in case of non-mpi)
    node_vs_pressurepoint.resize(points);

#if HAVE_MPI
    // Distribute work load over mpi nodes.
    for (int idx=0; idx < points; ++idx) {
        // Ensure master node gets equal or less work than the other nodes, since
        // master node also computes single phase perm.
        node_vs_pressurepoint[idx] = (mpi_nodecount-1) - idx % mpi_nodecount;
    }
#endif

    const std::vector<int>& ecl_idx = upscaler.grid().globalCell();
    clock_t start_upscale_wallclock = clock();

    double waterVolumeLF;
    // Now loop through the vector of capillary pressure points that
    // this node should compute.
    for (int pointidx = 0; pointidx < points; ++pointidx) {

        // Should "I" (mpi-wise) compute this pressure point?
        if (node_vs_pressurepoint[pointidx] == mpi_rank) {

            double Ptestvalue = pressurePoints[pointidx];

            double accPhasePerm = 0.0;
            double accPhase2Perm = 0.0;

            double maxPhasePerm = 0.0;
            double maxPhase2Perm = 0.0;

            std::vector<double> phasePermValues, phase2PermValues;
            std::vector<std::vector<double> > phasePermValuesDiag, phase2PermValuesDiag;
            phasePermValues.resize(satnums.size());
            phasePermValuesDiag.resize(satnums.size());
            if (upscaleBothPhases) {
                phase2PermValues.resize(satnums.size());
                phase2PermValuesDiag.resize(satnums.size());
            }
            waterVolumeLF = 0.0;
            for (size_t i = 0; i < ecl_idx.size(); ++i) {
                unsigned int cell_idx = ecl_idx[i];
                double cellPhasePerm = minPerm;
                double cellPhase2Perm = minPerm;
                std::vector<double>  cellPhasePermDiag, cellPhase2PermDiag;
                cellPhasePermDiag.resize(3, minPerm);
                if (upscaleBothPhases)
                    cellPhase2PermDiag.resize(3, minPerm);

                if (satnums[cell_idx] > 0) { // handle "no rock" cells with satnum zero
                    double PtestvalueCell = Ptestvalue;
                    if (!dP.empty())
                        PtestvalueCell -= dP[cell_idx];

                    if (!anisotropic_input) {
                        double Jvalue = sqrt(perms[0][cell_idx] * milliDarcyToSqMetre/poros[cell_idx]) * PtestvalueCell;
                        double WaterSaturationCell
                            = InvJfunctions[int(satnums[cell_idx])-1].evaluate(Jvalue);
                        waterVolumeLF += WaterSaturationCell * cellPoreVolumes[cell_idx];

                        // Compute cell relative permeability. We use a lower cutoff-value as we
                        // easily divide by zero here.  When water saturation is
                        // zero, we get 'inf', which is circumvented by the cutoff value.
                        cellPhasePerm =
                            Krfunctions[0][0][int(satnums[cell_idx])-1].evaluate(WaterSaturationCell) *
                            perms[0][cell_idx];
                        if (upscaleBothPhases) {
                            cellPhase2Perm =
                                Krfunctions[0][1][int(satnums[cell_idx])-1].evaluate(WaterSaturationCell) *
                                perms[0][cell_idx];
                        }
                    }
                    else {
                        double WaterSaturationCell = SwPcfunctions[int(satnums[cell_idx])-1].evaluate(PtestvalueCell);
                        waterVolumeLF += WaterSaturationCell * cellPoreVolumes[cell_idx];

                        cellPhasePermDiag[0] = Krfunctions[0][0][int(satnums[cell_idx])-1].evaluate(WaterSaturationCell) *
                            perms[0][cell_idx];
                        cellPhasePermDiag[1] = Krfunctions[1][0][int(satnums[cell_idx])-1].evaluate(WaterSaturationCell) *
                            perms[1][cell_idx];
                        cellPhasePermDiag[2] = Krfunctions[2][0][int(satnums[cell_idx])-1].evaluate(WaterSaturationCell) *
                            perms[2][cell_idx];
                        if (upscaleBothPhases) {
                            cellPhase2PermDiag[0] = Krfunctions[0][1][int(satnums[cell_idx])-1].evaluate(WaterSaturationCell) *
                                perms[0][cell_idx];
                            cellPhase2PermDiag[1] = Krfunctions[1][1][int(satnums[cell_idx])-1].evaluate(WaterSaturationCell) *
                                perms[1][cell_idx];
                            cellPhase2PermDiag[2] = Krfunctions[2][1][int(satnums[cell_idx])-1].evaluate(WaterSaturationCell) *
                                perms[2][cell_idx];
                        }
                    }

                    phasePermValues[cell_idx] = cellPhasePerm;
                    phasePermValuesDiag[cell_idx] = cellPhasePermDiag;
                    maxPhasePerm = std::max(maxPhasePerm, cellPhasePerm);
                    maxPhasePerm = std::max(maxPhasePerm, *std::max_element(cellPhasePermDiag.begin(),
                                                                            cellPhasePermDiag.end()));
                    if (upscaleBothPhases) {
                        phase2PermValues[cell_idx] = cellPhase2Perm;
                        phase2PermValuesDiag[cell_idx] = cellPhase2PermDiag;
                        maxPhase2Perm = std::max(maxPhase2Perm, cellPhase2Perm);
                        maxPhase2Perm = std::max(maxPhase2Perm, *std::max_element(cellPhase2PermDiag.begin(),
                                                                                  cellPhase2PermDiag.end()));
                    }
                }
            }
            // Now we can determine the smallest permitted permeability we can calculate for

            // We have both a fixed bottom limit, as well as a possible higher limit determined
            // by a maximum allowable permeability.
            double minPhasePerm = std::max(maxPhasePerm/maxPermContrast, minPerm);
            double minPhase2Perm;
            if (upscaleBothPhases)
                minPhase2Perm = std::max(maxPhase2Perm/maxPermContrast, minPerm);

            // Now remodel the phase permeabilities obeying minPhasePerm
            SinglePhaseUpscaler::permtensor_t cellperm(3,3,nullptr);
            zero(cellperm);
            for (size_t i = 0; i < ecl_idx.size(); ++i) {
                unsigned int cell_idx = ecl_idx[i];
                zero(cellperm);
                if (!anisotropic_input) {
                    double cellPhasePerm = std::max(minPhasePerm, phasePermValues[cell_idx]);
                    accPhasePerm += cellPhasePerm;
                    double kval = std::max(minPhasePerm, cellPhasePerm);
                    cellperm(0,0) = kval;
                    cellperm(1,1) = kval;
                    cellperm(2,2) = kval;
                }
                else { // anisotropic_input
                    // Truncate values lower than minPhasePerm upwards.
                    phasePermValuesDiag[cell_idx][0] = std::max(minPhasePerm, phasePermValuesDiag[cell_idx][0]);
                    phasePermValuesDiag[cell_idx][1] = std::max(minPhasePerm, phasePermValuesDiag[cell_idx][1]);
                    phasePermValuesDiag[cell_idx][2] = std::max(minPhasePerm, phasePermValuesDiag[cell_idx][2]);
                    accPhasePerm += phasePermValuesDiag[cell_idx][0]; // not correct anyway
                    cellperm(0,0) = phasePermValuesDiag[cell_idx][0];
                    cellperm(1,1) = phasePermValuesDiag[cell_idx][1];
                    cellperm(2,2) = phasePermValuesDiag[cell_idx][2];
                }
                upscaler.setPermeability(i, cellperm);
            }

            // Output average phase perm, this is just a reality check so that we are not way off.
            //cout << ", Arith. mean phase perm = " << accPhasePerm/float(tesselatedCells) << " mD, ";

            //  Call single-phase upscaling code
            SinglePhaseUpscaler::permtensor_t phasePermTensor = upscaler.upscaleSinglePhase();

            // Now upscale phase permeability for phase 2
            SinglePhaseUpscaler::permtensor_t phase2PermTensor;
            if (upscaleBothPhases) {
                zero(cellperm);
                for (size_t i = 0; i < ecl_idx.size(); ++i) {
                    unsigned int cell_idx = ecl_idx[i];
                    zero(cellperm);
                    if (!anisotropic_input) {
                        double cellPhase2Perm = std::max(minPhase2Perm, phase2PermValues[cell_idx]);
                        accPhase2Perm += cellPhase2Perm;
                        double kval = std::max(minPhase2Perm, cellPhase2Perm);
                        cellperm(0,0) = kval;
                        cellperm(1,1) = kval;
                        cellperm(2,2) = kval;
                    }
                    else { // anisotropic_input
                        // Truncate values lower than minPhasePerm upwards.
                        phase2PermValuesDiag[cell_idx][0] = std::max(minPhase2Perm, phase2PermValuesDiag[cell_idx][0]);
                        phase2PermValuesDiag[cell_idx][1] = std::max(minPhase2Perm, phase2PermValuesDiag[cell_idx][1]);
                        phase2PermValuesDiag[cell_idx][2] = std::max(minPhase2Perm, phase2PermValuesDiag[cell_idx][2]);
                        accPhase2Perm += phase2PermValuesDiag[cell_idx][0]; // not correct anyway
                        cellperm(0,0) = phase2PermValuesDiag[cell_idx][0];
                        cellperm(1,1) = phase2PermValuesDiag[cell_idx][1];
                        cellperm(2,2) = phase2PermValuesDiag[cell_idx][2];
                    }
                    upscaler.setPermeability(i, cellperm);
                }
                phase2PermTensor = upscaler.upscaleSinglePhase();
            }

            // Here we recalculate the upscaled water saturation,
            // although it is already known when we asked for the
            // pressure point to compute for. Nonetheless, we
            // recalculate here to avoid any minor roundoff-error and
            // interpolation error (this means that the saturation
            // points are not perfectly uniformly distributed)
            WaterSaturation[pointidx] =  waterVolumeLF/poreVolume;


#ifdef HAVE_MPI
            std::cout << "Rank " << mpi_rank << ": ";
#endif
            std::cout << Ptestvalue << "\t" << WaterSaturation[pointidx];
            // Store and print phase-perm-result
            for (int voigtIdx=0; voigtIdx < tensorElementCount; ++voigtIdx) {
                PhasePerm[0][pointidx][voigtIdx] = getVoigtValue(phasePermTensor, voigtIdx);
                std::cout << "\t" << getVoigtValue(phasePermTensor, voigtIdx);
                if (upscaleBothPhases){
                    PhasePerm[1][pointidx][voigtIdx] = getVoigtValue(phase2PermTensor, voigtIdx);
                    std::cout << "\t" << getVoigtValue(phase2PermTensor, voigtIdx);
                }
            }
            std::cout << std::endl;
        }
    }

    clock_t finish_upscale_wallclock = clock();
    double timeused_upscale_wallclock = (double(finish_upscale_wallclock)-double(start_upscale_wallclock))/CLOCKS_PER_SEC;

    collectResults();

    // Average time pr. upscaling point:
#ifdef HAVE_MPI
    // Sum the upscaling time used by all processes
    double timeused_total;
    MPI_Reduce(&timeused_upscale_wallclock, &timeused_total, 1, MPI_DOUBLE,
               MPI_SUM, 0, MPI_COMM_WORLD);
    double avg_upscaling_time_pr_point = timeused_total/(double)points;
#else
    double avg_upscaling_time_pr_point = timeused_upscale_wallclock / (double)points;
#endif

    return std::make_tuple(timeused_upscale_wallclock, avg_upscaling_time_pr_point);
}

}
