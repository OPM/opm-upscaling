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
#include <iostream>
#include <exception>

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

}
