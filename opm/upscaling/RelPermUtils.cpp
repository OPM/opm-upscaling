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

}
