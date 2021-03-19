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

#if HAVE_CONFIG_H
#include "config.h"
#endif

#include <opm/upscaling/writeECLData.hpp>
#include <opm/parser/eclipse/Units/Units.hpp>

#include <opm/io/eclipse/OutputStream.hpp>
#include <opm/output/eclipse/DoubHEAD.hpp>
#include <opm/output/eclipse/InteHEAD.hpp>
#include <opm/output/eclipse/LogiHEAD.hpp>

#include <time.h>
#include <chrono>
#include <vector>

namespace {
  Opm::RestartIO::InteHEAD::Phases phases()
  {
    auto p = Opm::RestartIO::InteHEAD::Phases{};

    p.oil = p.water = true;
    p.gas = false;

    return p;
  }

  Opm::RestartIO::InteHEAD::TimePoint timeStamp(const time_t time_stamp)
  {
    return Opm::RestartIO::getSimulationTimePoint(time_stamp, 0.0);
  }

  std::vector<int>
  createIntehead(const int nx, const int ny, const int nz, const int nactive,
                 const time_t time_stamp)
  {
    const auto ih = ::Opm::RestartIO::InteHEAD{}
      .dimensions         (nx, ny, nz)
      .numActive          (nactive)
      .unitConventions    (Opm::UnitSystem::newMETRIC())
      .params_NWELZ       (155, 122, 130, 3)
      .wellTableDimensions({ 0, 0, 0, 0, 0, 0, 0 })
      .calendarDate       (timeStamp(time_stamp))
      .activePhases       (phases())
      .variousParam       (201702, 100);

    return ih.data();
  }

  std::vector<double>
  createDoubHead(const time_t time_stamp, const double elapsed)
  {
    const auto dur = std::chrono::duration<
      double, std::chrono::seconds::period>{ elapsed };

    const auto start = std::chrono::system_clock::from_time_t(time_stamp)
      - dur;

    auto ts = Opm::RestartIO::DoubHEAD::TimeStamp {
      std::chrono::time_point_cast<
        std::chrono::system_clock::time_point::duration
      >(start), dur
    };

    const auto dh = ::Opm::RestartIO::DoubHEAD{}
    .timeStamp(ts);

    return dh.data();
  }
}

namespace Opm {

  /*
    This function will write the data solution data in the DataMap
    @data as an ECLIPSE restart file, in addition to the solution
    fields the ECLIPSE restart file will have a minimum (hopefully
    sufficient) amount of header information.

    The ECLIPSE restart files come in two varietes; unified restart
    files which have all the report steps lumped together in one large
    chunk and non-unified restart files which are one file for each
    report step. In addition the files can be either formatted
    (i.e. ASCII) or unformatted (i.e. binary).

    The writeECLData() function has two hardcoded settings:
    'file_type' and 'fmt_file' which regulate which type of files the
    should be created. The extension of the files follow a convention:

      Unified, formatted    : .FUNRST
      Unified, unformatted  : .UNRST
      Multiple, formatted   : .Fnnnn
      Multiple, unformatted : .Xnnnn

    For the multiple files the 'nnnn' part is the report number,
    formatted with '%04d' format specifier. The writeECLData()
    function will use the ecl_util_alloc_filename() function to create
    an ECLIPSE filename according to this conventions.
  */

  void writeECLData(int nx, int ny, int nz, int nactive,
                    data::Solution data,
                    const int current_step,
                    const double current_time,
                    time_t current_posix_time,
                    const std::string& output_dir,
                    const std::string& base_name) {

    EclIO::OutputStream::Restart rstFile {
      EclIO::OutputStream::ResultSet { output_dir, base_name },
      current_step,
      EclIO::OutputStream::Formatted { false },
      EclIO::OutputStream::Unified   { true }
    };

    rstFile.write("INTEHEAD", createIntehead(nx, ny, nz, nactive, current_posix_time));
    rstFile.write("DOUBHEAD", createDoubHead(current_posix_time, current_time));

    rstFile.message("STARTSOL");

    for (const auto&  elm : data) {
      if (elm.second.target != data::TargetType::RESTART_SOLUTION)
        continue;

      const auto& val = elm.second.data;
      rstFile.write(elm.first, std::vector<float>(val.begin(), val.end()));
    }

    rstFile.message("ENDSOL");
  }
}
