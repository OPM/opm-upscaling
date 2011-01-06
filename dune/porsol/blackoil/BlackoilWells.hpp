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

#ifndef OPM_BLACKOILWELLS_HEADER_INCLUDED
#define OPM_BLACKOILWELLS_HEADER_INCLUDED

#include <dune/porsol/blackoil/fluid/BlackoilDefs.hpp>
#include <dune/common/ErrorMacros.hpp>
#include <dune/common/SparseTable.hpp>
#include <dune/common/fvector.hh>
#include <vector>

// Forward declaration.
namespace Dune
{
    class EclipseGridParser;
}


namespace Opm
{



    /// A class designed to encapsulate a set of rate- or
    /// pressure-controlled wells in the black-oil setting.
    class BlackoilWells : public BlackoilDefs
    {
    public:
        void init(const Dune::EclipseGridParser& parser);

        // Well-centric interface. Mostly used by pressure solver.
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

        // Updating rates and pressures after pressure solve.
        void update(int num_cells,
                    const std::vector<double>& well_pressures,
                    const std::vector<double>& well_fluxes);

        // Cell-centric interface. Mostly used by transport solver.
        double perforationPressure(int cell) const;
        double wellToReservoirFlux(int cell) const;
        Dune::FieldVector<double, 3> injectionMixture(int cell) const;

    private:
        struct WellData { WellType type; WellControl control; double target; };
        std::vector<WellData> well_data_;
        struct PerfData { int cell; double well_index; double pdelta; };
        Dune::SparseTable<PerfData> perf_data_;
        std::vector<double> well_flux_;
        std::vector<double> well_pressure_;
        Dune::FieldVector<double, 3> injection_mixture_;
    };


    // ------------ Method implementations --------------

    inline void BlackoilWells::init(const Dune::EclipseGridParser& parser)
    {
        // Temporary test hack.
        WellData w1 = { Injector, Pressure, 2e7 };
        WellData w2 = { Producer, Pressure, 1e7 };
        well_data_.push_back(w1);
        well_data_.push_back(w2);
        PerfData p1 = { 0, 1e-11, 0.0 };
        PerfData p2 = { 99, 1e-11, 0.0 };
        perf_data_.appendRow(&p1, &p1 + 1);
        perf_data_.appendRow(&p2, &p2 + 1);
        injection_mixture_ = 0.0;
        injection_mixture_[Gas] = 1.0;
    }

    inline int BlackoilWells::numWells() const
    {
        return well_data_.size();
    }

    inline BlackoilWells::WellType BlackoilWells::type(int wellnum) const
    {
        return well_data_[wellnum].type;
    }

    inline BlackoilWells::WellControl BlackoilWells::control(int wellnum) const
    {
        return well_data_[wellnum].control;
    }

    inline double BlackoilWells::target(int wellnum) const
    {
        return well_data_[wellnum].target;
    }

    inline int BlackoilWells::numPerforations(int wellnum) const
    {
        return perf_data_[wellnum].size();
    }

    inline int BlackoilWells::wellCell(int wellnum, int perfnum) const
    {
        return perf_data_[wellnum][perfnum].cell;
    }

    inline double BlackoilWells::wellIndex(int wellnum, int perfnum) const
    {
        return perf_data_[wellnum][perfnum].well_index;
    }

    inline double BlackoilWells::pressureDelta(int wellnum, int perfnum) const
    {
        return perf_data_[wellnum][perfnum].pdelta;
    }

    inline void BlackoilWells::update(int num_cells,
                                      const std::vector<double>& well_pressures,
                                      const std::vector<double>& well_fluxes)
    {
        // Input is per perforation, data members store for all cells.
        ASSERT(perf_data_.dataSize() == int(well_pressures.size()));
        well_pressure_.resize(num_cells, -1e100);
        well_flux_.resize(num_cells, 0.0);
        int pcount = 0;
        for (int w = 0; w < numWells(); ++w) {
            for (int perf = 0; perf < numPerforations(w); ++perf) {
                int cell = wellCell(w, perf);
                well_pressure_[cell] = well_pressures[pcount];
                well_flux_[cell] = well_fluxes[pcount];
                ++pcount;
            }
        }
        ASSERT(pcount == perf_data_.dataSize());
    }

    inline double BlackoilWells::wellToReservoirFlux(int cell) const
    {
        return well_flux_[cell];
    }

    inline double BlackoilWells::perforationPressure(int cell) const
    {
        return well_pressure_[cell];
    }

    inline Dune::FieldVector<double, 3> BlackoilWells::injectionMixture(int cell) const
    {
        return injection_mixture_;
    }

} // namespace Opm


#endif // OPM_BLACKOILWELLS_HEADER_INCLUDED
