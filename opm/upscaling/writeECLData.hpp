/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.
  Copyright 2012 Statoil ASA.

  This file is part of the Open Porous Media project (OPM).

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

#ifndef OPM_WRITEECLDATA_HEADER_INCLUDED
#define OPM_WRITEECLDATA_HEADER_INCLUDED

#include <string>

#include <time.h>

#include <opm/output/data/Solution.hpp>

namespace Opm
{

  // ECLIPSE output for general grids.
  void writeECLData(int nx, int ny, int nz,
                    int number_of_cells,
                    data::Solution,
                    const int current_step,
                    const double current_time,
                    time_t current_posix_time,
                    const std::string& output_dir,
                    const std::string& base_name);

}

#endif
