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

#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/SparseTable.hpp>
#include <vector>

// Forward declaration.
namespace Opm
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
        void init(const Opm::EclipseGridParser& parser);

        // Well-centric interface.
        int numWells() const;
        enum WellType { Injector, Producer };
        WellType type(int wellnum) const;
        enum WellControl { Rate, Pressure };
        WellControl control(int wellnum) const;
        double target(int wellnum) const;
        double referenceDepth(int wellnum) const;
        int numPerforations(int wellnum) const;
        int wellCell(int wellnum, int perfnum) const;
        double wellIndex(int wellnum, int perfnum) const;
        double pressureDelta(int wellnum, int perfnum) const;

        // Updating rates and pressures after pressure solve.
        void update(int num_cells,
                    const std::vector<double>& well_pressures,
                    const std::vector<double>& well_fluxes);

        // Cell-centric interface. Mostly used by transport solver.
        double perforationPressure(int cell) const;
        double wellToReservoirFlux(int cell) const;
        Dune::FieldVector<double, 3> injectionMixture(int cell) const;

    private:
	struct WellData { WellType type; WellControl control; double target; double reference_bhp_depth; };
        std::vector<WellData> well_data_;
        struct PerfData { int cell; double well_index; double pdelta; };
        Opm::SparseTable<PerfData> perf_data_;
    };


    // ------------ Method implementations --------------

    inline void Wells::init(const Opm::EclipseGridParser& parser)
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

    inline double Wells::referenceDepth(int wellnum) const
    {
        return well_data_[wellnum].reference_bhp_depth;
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

    inline void Wells::update(int /*num_cells*/,
                              const std::vector<double>& /*well_pressures*/,
                              const std::vector<double>& /*well_fluxes*/)
    {
    }

    inline double Wells::perforationPressure(int cell) const
    {
        THROW("Not implemented");
        return 0.0;
    }

    inline double Wells::wellToReservoirFlux(int cell) const
    {
        return 0.0;
    }

    inline Dune::FieldVector<double, 3> Wells::injectionMixture(int cell) const
    {
        THROW("Not implemented");
        return Dune::FieldVector<double, 3>(0.0);
    }


} // namespace Opm


#endif // OPM_WELLS_HEADER_INCLUDED
