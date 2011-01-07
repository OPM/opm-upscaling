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


template <class Grid, class Rock, class Fluid, class Wells>
class ExplicitCompositionalTransport : public BlackoilDefs
{
public:
    /// @brief
    ///    Default constructor. Does nothing.
    ExplicitCompositionalTransport()
        : pgrid_(0), prock_(0), pfluid_(0)
    {
    }

    void init(const Dune::parameter::ParameterGroup& param)
    {
    }

    void setup(const Grid& grid,
               const Rock& rock,
               const Fluid& fluid,
               const Wells& wells)
    {
        pgrid_ = &grid;
        prock_ = &rock;
        pfluid_ = &fluid;
        pwells_ = &wells;
    }


    /// Return value is the time actually used, it may be smaller than dt if
    /// we stop due to unacceptable volume discrepancy.
    double transport(const PhaseVec& external_pressure,
                     const CompVec& external_composition,
                     const std::vector<double>& face_flux,
                     const std::vector<PhaseVec>& cell_pressure,
                     const std::vector<PhaseVec>& face_pressure,
                     const double dt,
                     const double voldisclimit,
                     std::vector<CompVec>& cell_z)
    {
        int num_cells = pgrid_->numCells();
        std::vector<CompVec> comp_change;
        std::vector<double> cell_outflux;
        std::vector<double> cell_max_ff_deriv;
        double cur_time = 0.0;
        updateFluidProperties(cell_pressure, face_pressure, cell_z,
                              external_pressure, external_composition);
        if (!volumeDiscrepancyAcceptable(voldisclimit)) {
            return 0.0;
        }
        std::vector<CompVec> cell_z_start;
        while (cur_time < dt) {
            cell_z_start = cell_z;
            computeChange(face_flux, comp_change, cell_outflux, cell_max_ff_deriv);
            double min_time = 1e100;
            for (int cell = 0; cell < num_cells; ++cell) {
                double time = (prock_->porosity(cell)*pgrid_->cellVolume(cell))/(cell_outflux[cell]*cell_max_ff_deriv[cell]);
                min_time = std::min(time, min_time);
            }
            min_time *= 0.49; // Semi-random CFL factor... \TODO rigorize
            double step_time = dt - cur_time;
            if (min_time < step_time) {
                step_time = min_time;
                cur_time += min_time;
            } else {
                cur_time = dt;
            }
            std::cout << "Taking step in explicit transport solver: " << step_time << std::endl;
            for (int cell = 0; cell < num_cells; ++cell) {
                comp_change[cell] *= (step_time/prock_->porosity(cell));
                cell_z[cell] += comp_change[cell];
            }
            // After changing z, we recompute fluid properties.
            updateFluidProperties(cell_pressure, face_pressure, cell_z,
                                  external_pressure, external_composition);
            bool ok = volumeDiscrepancyAcceptable(voldisclimit);
            if (!ok) {
                // Roll back to last ok step.
                cell_z = cell_z_start;
                cur_time -= step_time;
                return cur_time;
            }
        }
        return dt;
    }


private: // Data
    const Grid* pgrid_;
    const Rock* prock_;
    const Fluid* pfluid_;
    const Wells* pwells_;
    typename Fluid::FluidData fluid_data_;
    struct TransportFluidData
    {
        PhaseVec saturation;
        PhaseVec fractional_flow;
        std::tr1::array<CompVec, numPhases> phase_to_comp;
        PhaseVec relperm;
        PhaseVec viscosity;
    };
    TransportFluidData bdy_;
    std::vector<int> perf_cells_;
    std::vector<double> perf_flow_;
    std::vector<TransportFluidData> perf_props_;

private: // Methods

    TransportFluidData computeProps(const PhaseVec& pressure,
                                    const CompVec& composition)
    {
        BlackoilFluid::FluidState state = pfluid_->computeState(pressure, composition);
        TransportFluidData data;
        data.saturation = state.saturation_;
        double total_mobility = 0.0;
        for (int phase = 0; phase < numPhases; ++phase) {
            total_mobility += state.mobility_[phase];
        }
        data.fractional_flow = state.mobility_;
        data.fractional_flow /= total_mobility;
        std::copy(state.phase_to_comp_, state.phase_to_comp_ + numComponents*numPhases,
                  &data.phase_to_comp[0][0]);
        data.relperm = state.relperm_;
        data.viscosity = state.viscosity_;
        return data;
    }

    void updateFluidProperties(const std::vector<PhaseVec>& cell_pressure,
                               const std::vector<PhaseVec>& face_pressure,
                               const std::vector<CompVec>& cell_z,
                               const PhaseVec& external_pressure,
                               const CompVec& external_composition)
    {
        // Properties in reservoir.
        const double dummy_dt = 1.0;
        fluid_data_.compute(*pgrid_, *prock_, *pfluid_, cell_pressure, face_pressure, cell_z, external_composition, dummy_dt);

        // Properties on boundary. \TODO no need to ever recompute this.
        bdy_ = computeProps(external_pressure, external_composition);

        // Properties at well perforations.
        // \TODO only need to recompute this once per pressure update.
        // No, that is false, at production perforations the cell z is
        // used, which may change every step.
        perf_cells_.clear();
        perf_flow_.clear();
        perf_props_.clear();
        int num_cells = pgrid_->numCells();
        for (int cell = 0; cell < num_cells; ++cell) {
            double flow = pwells_->wellToReservoirFlux(cell);
            if (flow != 0.0) {
                perf_cells_.push_back(cell);
                perf_flow_.push_back(flow);
                // \TODO handle capillary in perforation pressure below?
                PhaseVec well_pressure = flow > 0.0 ? PhaseVec(pwells_->perforationPressure(cell)) : cell_pressure[cell];
                CompVec well_mixture = flow > 0.0 ? pwells_->injectionMixture(cell) : cell_z[cell];
                perf_props_.push_back(computeProps(well_pressure, well_mixture));
            }
        }
    }


    bool volumeDiscrepancyAcceptable(const double voldisclimit) const
    {
        double rel_voldiscr = *std::max_element(fluid_data_.relvoldiscr.begin(), fluid_data_.relvoldiscr.end());
        if (rel_voldiscr > voldisclimit) {
            std::cout << "    Relative volume discrepancy too large: " << rel_voldiscr << std::endl;
            return false;
        } else {
            return true;
        }
    }


    void computeChange(const std::vector<double>& face_flux,
                       std::vector<CompVec>& comp_change,
                       std::vector<double>& cell_outflux,
                       std::vector<double>& cell_max_ff_deriv)
    {
        int num_cells = pgrid_->numCells();
        comp_change.clear();
        CompVec zero(0.0);
        comp_change.resize(num_cells, zero);
        cell_outflux.clear();
        cell_outflux.resize(num_cells, 0.0);
        cell_max_ff_deriv.clear();
        cell_max_ff_deriv.resize(num_cells, 0.0);
        for (int face = 0; face < pgrid_->numFaces(); ++face) {
            // Set up needed quantities.
            int c0 = pgrid_->faceCell(face, 0);
            int c1 = pgrid_->faceCell(face, 1);
            int upwind_cell = (face_flux[face] > 0.0) ? c0 : c1;
            int downwind_cell = (face_flux[face] > 0.0) ? c1 : c0;
            PhaseVec upwind_sat = upwind_cell < 0 ? bdy_.saturation : fluid_data_.saturation[upwind_cell];
            PhaseVec upwind_relperm = upwind_cell < 0 ? bdy_.relperm : fluid_data_.rel_perm[upwind_cell];
            PhaseVec upwind_viscosity = upwind_cell < 0 ? bdy_.viscosity : fluid_data_.viscosity[upwind_cell];
            PhaseVec upwind_ff = upwind_cell < 0 ? bdy_.fractional_flow : fluid_data_.frac_flow[upwind_cell];
            PhaseVec phase_flux(upwind_ff);
            phase_flux *= face_flux[face];
            CompVec change(0.0);

            // Estimate max derivative of ff.
            double face_max_ff_deriv = 0.0;
            if (downwind_cell >= 0) { // Only contribution on inflow and internal faces.
                // Evaluating all functions at upwind viscosity.
                PhaseVec downwind_mob(0.0);
                double downwind_totmob = 0.0;
                for (int phase = 0; phase < numPhases; ++phase) {
                    downwind_mob[phase] = fluid_data_.rel_perm[downwind_cell][phase]/upwind_viscosity[phase];
                    downwind_totmob += downwind_mob[phase];
                }
                PhaseVec downwind_ff = downwind_mob;
                downwind_ff /= downwind_totmob;
                PhaseVec ff_diff = upwind_ff;
                ff_diff -= downwind_ff;
                for (int phase = 0; phase < numPhases; ++phase) {
                    if (std::fabs(ff_diff[phase]) > 1e-10) {
                        double ff_deriv = ff_diff[phase]/(upwind_sat[phase] - fluid_data_.saturation[downwind_cell][phase]);
                        ASSERT(ff_deriv >= 0.0);
                        face_max_ff_deriv = std::max(face_max_ff_deriv, ff_deriv);
                    }
                }
            }

            // Compute z change.
            for (int phase = 0; phase < numPhases; ++phase) {
                CompVec z_in_phase = bdy_.phase_to_comp[phase];
                if (upwind_cell >= 0) {
                    for (int comp = 0; comp < numComponents; ++comp) {
                        z_in_phase[comp] = fluid_data_.cellA[numPhases*numComponents*upwind_cell + numComponents*phase + comp];
                    }
                }
                z_in_phase *= phase_flux[phase];
                change += z_in_phase;
            }

            // Update output variables.
            if (upwind_cell >= 0) {
                cell_outflux[upwind_cell] += std::fabs(face_flux[face]);
            }
            if (c0 >= 0) {
                comp_change[c0] -= change;
                cell_max_ff_deriv[c0] = std::max(cell_max_ff_deriv[c0], face_max_ff_deriv);
            }
            if (c1 >= 0) {
                comp_change[c1] += change;
                cell_max_ff_deriv[c1] = std::max(cell_max_ff_deriv[c1], face_max_ff_deriv);
            }
        }

        // Done with all faces, now deal with well perforations.
        int num_perf = perf_cells_.size();
        for (int perf = 0; perf < num_perf; ++perf) {
            int cell = perf_cells_[perf];
            double flow = perf_flow_[perf];
            ASSERT(flow != 0.0);
            const TransportFluidData& fl = perf_props_[perf];
            // For injection, phase volumes depend on fractional
            // flow of injection perforation, for production we
            // use the fractional flow of the producing cell.
            PhaseVec phase_flux = flow > 0.0 ? fl.fractional_flow : fluid_data_.frac_flow[cell];
            phase_flux *= flow;
            // Conversion to mass flux is given at perforation state
            // if injector, cell state if producer (this is ensured by
            // updateFluidProperties()).
            CompVec change(0.0);
            for (int phase = 0; phase < numPhases; ++phase) {
                CompVec z_in_phase = fl.phase_to_comp[phase];
                z_in_phase *= phase_flux[phase];
                change += z_in_phase;
            }
            comp_change[cell] += change;
        }
    }


};


} // namespace Opm


#endif // OPM_COMPONENTTRANSPORT_HEADER_INCLUDED
