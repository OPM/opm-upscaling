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
                   const double external_pressure,
                   const CompVec& external_composition,
                   const std::vector<double>& face_flux,
                   EquationOfState& eos,
                   const std::vector<PhaseVec>& phase_pressure,
                   const double dt,
                   std::vector<CompVec>& comp_amount)
    {
        int num_cells = grid.numCells();
        std::vector<PhaseVec> fractional_flow(num_cells);
        std::vector<typename EquationOfState::ComponentInPhaseMatrix> comp_in_phase(num_cells);
        for (int cell = 0; cell < num_cells; ++cell) {
            eos.compute(phase_pressure[cell], comp_amount[cell]);
            fractional_flow[cell] = eos.fractionalFlow();
            comp_in_phase[cell] = eos.compositionByPhase();
        }
        eos.compute(PhaseVec(external_pressure), external_composition);
        transportImpl(grid, rock, face_flux, fractional_flow, comp_in_phase,
                      eos.fractionalFlow(), eos.compositionByPhase(), dt, comp_amount);
    }



    void transportImpl(const Grid& grid,
                       const Rock& rock,
                       const std::vector<double>& face_flux,
                       const std::vector<PhaseVec>& fractional_flow,
                       const std::vector<typename EquationOfState::ComponentInPhaseMatrix>& comp_in_phase,
                       const PhaseVec& bdy_fractional_flow,
                       const typename EquationOfState::ComponentInPhaseMatrix& bdy_comp_in_phase,
                       const double dt,
                       std::vector<CompVec>& comp_amount)
    {
        std::vector<CompVec> comp_change;
        computeChange(grid, face_flux, fractional_flow, comp_in_phase,
                      bdy_fractional_flow, bdy_comp_in_phase, comp_change);
        for (int cell = 0; cell < grid.numCells(); ++cell) {
            comp_change[cell] *= (dt/rock.porosity(cell));
            comp_amount[cell] += comp_change[cell];
        }
    }


    void computeChange(const Grid& grid,
                       const std::vector<double>& face_flux,
                       const std::vector<PhaseVec>& fractional_flow,
                       const std::vector<typename EquationOfState::ComponentInPhaseMatrix>& comp_in_phase,
                       const PhaseVec& bdy_fractional_flow,
                       const typename EquationOfState::ComponentInPhaseMatrix& bdy_comp_in_phase,
                       std::vector<CompVec>& comp_change)
    {
        comp_change.clear();
        CompVec zero(0.0);
        comp_change.resize(grid.numCells(), zero);
        for (int face = 0; face < grid.numFaces(); ++face) {
            int c0 = grid.faceCell(face, 0);
            int c1 = grid.faceCell(face, 1);
            CompVec change(0.0);
            for (int phase = 0; phase < numPhases; ++phase) {
                int upwind_cell = (face_flux[face] > 0.0) ? c0 : c1;
                double ff = upwind_cell < 0 ? bdy_fractional_flow[phase] : fractional_flow[upwind_cell][phase];
                double phase_flux = face_flux[face]*ff;
                CompVec z_in_phase = upwind_cell < 0 ? bdy_comp_in_phase[phase] : comp_in_phase[upwind_cell][phase];
                z_in_phase *= phase_flux;
                    change += z_in_phase;
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
