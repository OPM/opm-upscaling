/*
  Copyright 2011 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_WELLS_HEADER_INCLUDED
#define OPM_WELLS_HEADER_INCLUDED

#include <dune/common/ErrorMacros.hpp>
#include <dune/common/SparseTable.hpp>
#include <vector>

// Forward declaration.
namespace Dune
{
    class EclipseGridParser;
}


namespace Opm
{



    /// A class designed to encapsulate a set of rate- or
    /// pressure-controlled wells.
    class Wells
    {
    public:
        void init(const Dune::EclipseGridParser& parser);

        // Well-centric interface.
        int numWells() const;
        enum WellType { Injector, Producer };
        WellType type(int wellnum) const;
        enum WellControl { Rate, Pressure };
        WellControl control(int wellnum) const;
        double target(int wellnum) const;
        int numPerforations(int wellnum) const;
        int wellCell(int wellnum, int perfnum) const;
        double wellIndex(int wellnum, int perfnum) const;
        double pressureDelta(int wellnum, int perfnum) const;

        // Cell-centric interface.
        double wellOutflowRate(int cell) const;
    private:
        struct WellData { WellType type; WellControl control; double target; };
        std::vector<WellData> well_data_;
        struct PerfData { int cell; double well_index; double pdelta; };
        Dune::SparseTable<PerfData> perf_data_;
    };


    // ------------ Method implementations --------------

    inline void Wells::init(const Dune::EclipseGridParser& parser)
    {
    }

    inline int Wells::numWells() const
    {
        return well_data_.size();
    }

    inline Wells::WellType Wells::type(int wellnum) const
    {
        return well_data_[wellnum].type;
    }

    inline Wells::WellControl Wells::control(int wellnum) const
    {
        return well_data_[wellnum].control;
    }

    inline double Wells::target(int wellnum) const
    {
        return well_data_[wellnum].target;
    }

    inline int Wells::numPerforations(int wellnum) const
    {
        return perf_data_[wellnum].size();
    }

    inline int Wells::wellCell(int wellnum, int perfnum) const
    {
        return perf_data_[wellnum][perfnum].cell;
    }

    inline double Wells::wellIndex(int wellnum, int perfnum) const
    {
        return perf_data_[wellnum][perfnum].well_index;
    }

    inline double Wells::pressureDelta(int wellnum, int perfnum) const
    {
        return perf_data_[wellnum][perfnum].pdelta;
    }

    inline double Wells::wellOutflowRate(int cell) const
    {
        THROW("Not implemented");
        return 0.0;
    }


} // namespace Opm


#endif // OPM_WELLS_HEADER_INCLUDED
