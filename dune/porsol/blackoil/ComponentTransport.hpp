/*
  Copyright 2010 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_COMPONENTTRANSPORT_HEADER_INCLUDED
#define OPM_COMPONENTTRANSPORT_HEADER_INCLUDED

#include <tr1/array>
#include <vector>
#include <dune/porsol/blackoil/fluid/BlackoilDefs.hpp>
#include <dune/porsol/blackoil/BlackoilFluid.hpp>
#include <dune/common/param/ParameterGroup.hpp>

namespace Opm
{


class EquationOfStateBlackOil : public BlackoilDefs
{
public:
    typedef std::tr1::array<CompVec, numPhases> ComponentInPhaseMatrix;

    EquationOfStateBlackOil(const BlackoilFluid& fluid)
        : fluid_(fluid)
    {
    }

    void compute(const PhaseVec& phase_pressure, const CompVec& z)
    {
        BlackoilFluid::FluidState state = fluid_.computeState(phase_pressure, z);
        double total_mobility = 0.0;
        for (int phase = 0; phase < numPhases; ++phase) {
            total_mobility += state.mobility_[phase];
        }
        fractional_flow_ = state.mobility_;
        fractional_flow_ /= total_mobility;
        Dune::SharedFortranMatrix A(numComponents, numPhases, state.phase_to_comp_);
        Dune::SharedCMatrix cbp(numComponents, numPhases, &comp_by_phase_[0][0]);
        cbp = A; // Converting between orderings.
    }

    const PhaseVec& fractionalFlow() const
    {
        return fractional_flow_;
    }

    const ComponentInPhaseMatrix& compositionByPhase()
    {
        return comp_by_phase_;
    }
private:
    const BlackoilFluid& fluid_;
    PhaseVec fractional_flow_;
    ComponentInPhaseMatrix comp_by_phase_;
};


template <class Grid, class Rock, class EquationOfState>
class ExplicitCompositionalTransport : public BlackoilDefs
{
public:
    void init(const Dune::parameter::ParameterGroup& param)
    {
    }

    void transport(const Grid& grid,
                   const Rock& rock,
                   const PhaseVec& external_pressure,
                   const CompVec& external_composition,
                   const std::vector<double>& face_flux,
                   EquationOfState& eos,
                   const std::vector<PhaseVec>& phase_pressure,
                   const double dt,
                   std::vector<CompVec>& comp_amount)
    {
        int num_cells = grid.numCells();
        std::vector<CompVec> comp_change;
        std::vector<double> cell_outflux;
        double cur_time = 0.0;
        while (cur_time < dt) {
            updateFluidProperties(phase_pressure, comp_amount, external_pressure, external_composition, eos);
            computeChange(grid, face_flux, comp_change, cell_outflux);
            double min_time = 1e100;
            for (int cell = 0; cell < num_cells; ++cell) {
                double time = (rock.porosity(cell)*grid.cellVolume(cell))/cell_outflux[cell];
                min_time = std::min(time, min_time);
            }
            double step_time = dt - cur_time;
            if (min_time < step_time) {
                step_time = min_time;
                cur_time += min_time;
            } else {
                cur_time = dt;
            }
            std::cout << "Taking step in explicit transport solver: " << step_time << std::endl;
            for (int cell = 0; cell < num_cells; ++cell) {
                comp_change[cell] *= (step_time/rock.porosity(cell));
                comp_amount[cell] += comp_change[cell];
            }
        }
    }


private: // Data

    PhaseVec bdy_fractional_flow_;
    typename EquationOfState::ComponentInPhaseMatrix bdy_comp_in_phase_;
    std::vector<PhaseVec> fractional_flow_;
    std::vector<typename EquationOfState::ComponentInPhaseMatrix> comp_in_phase_;


private: // Methods

    void updateFluidProperties(const std::vector<PhaseVec>& phase_pressure,
                               const std::vector<CompVec>& comp_amount,
                               const PhaseVec& external_pressure,
                               const CompVec& external_composition,
                               EquationOfState& eos)
    {
        int num_cells = phase_pressure.size();
        fractional_flow_.resize(num_cells);
        comp_in_phase_.resize(num_cells);
        for (int cell = 0; cell < num_cells; ++cell) {
            eos.compute(phase_pressure[cell], comp_amount[cell]);
            fractional_flow_[cell] = eos.fractionalFlow();
            comp_in_phase_[cell] = eos.compositionByPhase();
        }
        eos.compute(external_pressure, external_composition);
        bdy_fractional_flow_ = eos.fractionalFlow();
        bdy_comp_in_phase_ = eos.compositionByPhase();
    }



    void computeChange(const Grid& grid,
                       const std::vector<double>& face_flux,
                       std::vector<CompVec>& comp_change,
                       std::vector<double>& cell_outflux)
    {
        comp_change.clear();
        CompVec zero(0.0);
        comp_change.resize(grid.numCells(), zero);
        cell_outflux.clear();
        cell_outflux.resize(grid.numCells(), 0.0);
        for (int face = 0; face < grid.numFaces(); ++face) {
            int c0 = grid.faceCell(face, 0);
            int c1 = grid.faceCell(face, 1);
            CompVec change(0.0);
            for (int phase = 0; phase < numPhases; ++phase) {
                int upwind_cell = (face_flux[face] > 0.0) ? c0 : c1;
                double ff = upwind_cell < 0 ? bdy_fractional_flow_[phase] : fractional_flow_[upwind_cell][phase];
                double phase_flux = face_flux[face]*ff;
                CompVec z_in_phase = upwind_cell < 0 ? bdy_comp_in_phase_[phase] : comp_in_phase_[upwind_cell][phase];
                z_in_phase *= phase_flux;
                change += z_in_phase;
                if (upwind_cell >= 0) {
                    cell_outflux[upwind_cell] += std::fabs(phase_flux);
                }
            }
            if (c0 >= 0) {
                comp_change[c0] -= change;
            }
            if (c1 >= 0) {
                comp_change[c1] += change;
            }
        }
    }
};


} // namespace Opm


#endif // OPM_COMPONENTTRANSPORT_HEADER_INCLUDED
