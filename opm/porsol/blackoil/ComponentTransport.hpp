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
#include <opm/porsol/blackoil/fluid/BlackoilDefs.hpp>
#include <opm/porsol/blackoil/BlackoilFluid.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>

#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <limits.h>
#include <assert.h>
#include <stdio.h>

namespace Opm
{


template <class Grid, class Rock, class Fluid, class Wells>
class ExplicitCompositionalTransport : public BlackoilDefs
{
public:
    /// @brief
    ///    Default constructor. Does nothing.
    ExplicitCompositionalTransport()
        : pgrid_(0), prock_(0), pfluid_(0), pwells_(0), ptrans_(0),
          min_surfvol_threshold_(0.0),
          single_step_only_(false),
          min_vtime_(0.0)
    {
    }

    void init(const Opm::parameter::ParameterGroup& param)
    {
        min_surfvol_threshold_ = param.getDefault("min_surfvol_threshold", min_surfvol_threshold_);
        single_step_only_ = param.getDefault("single_step_only", single_step_only_);
        min_vtime_ = param.getDefault("min_vtime",  min_vtime_);
    }

    void setup(const Grid& grid,
               const Rock& rock,
               const Fluid& fluid,
               const Wells& wells,
               const std::vector<double>& face_trans,
               const typename Grid::Vector& gravity)
    {
        pgrid_ = &grid;
        prock_ = &rock;
        pfluid_ = &fluid;
        pwells_ = &wells;
        ptrans_ = &face_trans;
        gravity_ = gravity;
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
        std::vector<double> cell_gravflux;
        std::vector<double> cell_max_mob_deriv;
        double cur_time = 0.0;
        updateFluidProperties(cell_pressure, face_pressure, cell_z,
                              external_pressure, external_composition);
//         if (!volumeDiscrepancyAcceptable(voldisclimit)) {
//             return 0.0;
//         }
        std::vector<CompVec> cell_z_start;
        std::cout << "Transport solver target time: " << dt << std::endl;
        std::cout << "   Step               Stepsize           Remaining time\n";
        int count = 0;
        while (cur_time < dt) {
            cell_z_start = cell_z;
            computeChange(face_flux, comp_change, cell_outflux, cell_max_ff_deriv);
            double min_time = 1e100;
            for (int cell = 0; cell < num_cells; ++cell) {
                double pvol = prock_->porosity(cell)*pgrid_->cellVolume(cell);
                double vtime = pvol/(cell_outflux[cell]*cell_max_ff_deriv[cell]);
                double gtime = 1e100; // No working CFL for gravity yet.
                double max_nonzero_time = 1e100;
                for (int comp = 0; comp < numComponents; ++comp) {
                    if (comp_change[cell][comp] < 0.0) {
                        if (cell_z[cell][comp] > min_surfvol_threshold_) {
                            max_nonzero_time = std::min(max_nonzero_time,
                                                        -cell_z[cell][comp]*pvol/comp_change[cell][comp]);
                        } else {
                            comp_change[cell][comp] = 0.0;
                            cell_z[cell][comp] = 0.0;
                        }
                    }
                }
                vtime = std::max(vtime,min_vtime_);
                double time = std::min(std::min(vtime, gtime), max_nonzero_time);
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
            // Check if remaining number of steps is excessive.
            if ((dt - cur_time)/step_time > 10000) {
                std::cout << "Collapsing transport stepsize detected." << std::endl;
                return cur_time;
            }
            std::cout.precision(10);
            std::cout << std::setw(6) << count++
                      << std::setw(24) << step_time
                      << std::setw(24) << dt - cur_time << std::endl;
            std::cout.precision(16);
            for (int cell = 0; cell < num_cells; ++cell) {
                double pv = pgrid_->cellVolume(cell)*prock_->porosity(cell);
                comp_change[cell] *= (step_time/pv);
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
            if (single_step_only_) {
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
    const std::vector<double>* ptrans_;
    typename Grid::Vector gravity_;

    struct AllTransportFluidData : public Opm::AllFluidData
    {
        std::vector<double> total_mobility;
        std::vector<PhaseVec> fractional_flow;

        template <class G, class R>
        void computeNew(const G& grid,
                        const R& rock,
                        const BlackoilFluid& fluid,
                        const typename Grid::Vector gravity,
                        const std::vector<PhaseVec>& cell_pressure,
                        const std::vector<PhaseVec>& face_pressure,
                        const std::vector<CompVec>& cell_z,
                        const CompVec& bdy_z,
                        const double dt)
        {
            Opm::AllFluidData::computeNew(grid, rock, fluid, gravity,
                                          cell_pressure, face_pressure,
                                          cell_z, bdy_z, dt);
            int num = grid.numCells();
            total_mobility.resize(num);
            fractional_flow.resize(num);
#pragma omp parallel for
            for (int i = 0; i < num; ++i) {
                total_mobility[i] = 0.0;
                for (int phase = 0; phase < numPhases; ++phase) {
                    total_mobility[i] += cell_data.mobility[i][phase];
                }
                fractional_flow[i] = cell_data.mobility[i];
                fractional_flow[i] *= (1.0/total_mobility[i]);
            }
        }
    };
    AllTransportFluidData fluid_data_;
    struct TransportFluidData
    {
        PhaseVec saturation;
        PhaseVec mobility;
        PhaseVec fractional_flow;
        PhaseToCompMatrix phase_to_comp;
        PhaseVec relperm;
        PhaseVec viscosity;
    };
    TransportFluidData bdy_;
    std::vector<int> perf_cells_;
    std::vector<double> perf_flow_;
    std::vector<TransportFluidData> perf_props_;
    double min_surfvol_threshold_;
    bool single_step_only_;
    double min_vtime_;

private: // Methods

    TransportFluidData computeProps(const PhaseVec& pressure,
                                    const CompVec& composition)
    {
        BlackoilFluid::FluidState state = pfluid_->computeState(pressure, composition);
        TransportFluidData data;
        data.saturation = state.saturation_;
        data.mobility = state.mobility_;
        double total_mobility = 0.0;
        for (int phase = 0; phase < numPhases; ++phase) {
            total_mobility += state.mobility_[phase];
        }
        data.fractional_flow = state.mobility_;
        data.fractional_flow /= total_mobility;
        data.phase_to_comp = state.phase_to_comp_;
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
        fluid_data_.computeNew(*pgrid_, *prock_, *pfluid_, gravity_,
                               cell_pressure, face_pressure, cell_z, external_composition, dummy_dt);

        // Properties on boundary. \TODO no need to ever recompute this.
        bdy_ = computeProps(external_pressure, external_composition);

        // Properties at well perforations.
        // \TODO only need to recompute this once per pressure update.
        // No, that is false, at production perforations the cell z is
        // used, which may change every step.
        perf_cells_.clear();
        perf_flow_.clear();
        perf_props_.clear();
        Wells::WellReport::report()->clearAll();
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
                Wells::WellReport::report()->perfPressure.push_back(pwells_->perforationPressure(cell));
                Wells::WellReport::report()->cellPressure.push_back(cell_pressure[cell][0]);
                Wells::WellReport::report()->cellId.push_back(cell);
            }
        }
    }




    void cellData(int cell, TransportFluidData& tfd) const
    {
        tfd.saturation = fluid_data_.cell_data.saturation[cell];
        tfd.mobility = fluid_data_.cell_data.mobility[cell];
        tfd.fractional_flow = fluid_data_.fractional_flow[cell];
        tfd.phase_to_comp = fluid_data_.cell_data.state_matrix[cell];
        tfd.relperm = fluid_data_.cell_data.relperm[cell];
        tfd.viscosity = fluid_data_.cell_data.viscosity[cell];
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
        const int num_cells = pgrid_->numCells();
        comp_change.clear();
        CompVec zero(0.0);
        comp_change.resize(num_cells, zero);
        cell_outflux.clear();
        cell_outflux.resize(num_cells, 0.0);
        cell_max_ff_deriv.clear();
        cell_max_ff_deriv.resize(num_cells, 0.0);
        CompVec surf_dens = pfluid_->surfaceDensities();
	const int num_faces = pgrid_->numFaces();
	std::vector<CompVec> face_change(num_faces);
	std::vector<double> faces_max_ff_deriv(num_faces);
#pragma omp parallel for
        for (int face = 0; face < num_faces; ++face) {
            // Compute phase densities on face.
            PhaseVec phase_dens(0.0);
            for (int phase = 0; phase < numPhases; ++phase) {
                const double* At = &fluid_data_.face_data.state_matrix[face][0][0]; // Already transposed since in Fortran order...
                for (int comp = 0; comp < numPhases; ++comp) {
                    phase_dens[phase] += At[numPhases*phase + comp]*surf_dens[comp];
                }
            }
            // Collect data from adjacent cells (or boundary).
            int c[2];
            TransportFluidData d[2];
            for (int ix = 0; ix < 2; ++ix) {
                c[ix] = pgrid_->faceCell(face, ix);
                if (c[ix] >= 0) {
                    cellData(c[ix], d[ix]);
                } else {
                    d[ix] = bdy_;
                }
            }
            // Compute upwind directions.
            int upwind_dir[numPhases] = { 0, 0, 0 };
            PhaseVec vstar(face_flux[face]);
            // double gravity_flux = gravity_*pgrid_->faceNormal(face)*pgrid_->faceArea(face);

            typename Grid::Vector centroid_diff = c[0] >= 0 ? pgrid_->cellCentroid(c[0]) : pgrid_->faceCentroid(face);
            centroid_diff -= c[1] >= 0 ? pgrid_->cellCentroid(c[1]) : pgrid_->faceCentroid(face);
            double gravity_flux = gravity_*centroid_diff*ptrans_->operator[](face);
            PhaseVec rho_star(phase_dens);
            process_face(&(d[0].mobility[0]), &(d[1].mobility[0]),
                         &vstar[0], gravity_flux, numPhases, &rho_star[0], upwind_dir);

            // Compute phase fluxes.
            PhaseVec phase_mob;
            double tot_mob = 0.0;
            for (int phase = 0; phase < numPhases; ++phase) {
                phase_mob[phase] = d[upwind_dir[phase] - 1].mobility[phase];
                tot_mob += phase_mob[phase];
            }
            PhaseVec ff = phase_mob;
            ff /= tot_mob;
            PhaseVec phase_flux = ff;
            phase_flux *= face_flux[face];
            // Until we have proper bcs for transport, assume no gravity flow across bdys.
            if (gravity_flux != 0.0 && c[0] >= 0 && c[1] >= 0) {
                // Gravity contribution.
                double omega = ff*phase_dens;
                for (int phase = 0; phase < numPhases; ++phase) {
                    double gf = (phase_dens[phase] - omega)*gravity_flux;
                    phase_flux[phase] -= phase_mob[phase]*gf;
                }
            }

            // Estimate max derivative of ff.
            double face_max_ff_deriv = 0.0;
            // Only using total flux upwinding for this purpose.
            // Aim is to reproduce old results first. \TODO fix, include gravity.
            int downwind_cell = c[upwind_dir[0]%2]; // Keep in mind that upwind_dir[] \in {1, 2}
            if (downwind_cell >= 0) { // Only contribution on inflow and internal faces.
                // Evaluating all functions at upwind viscosity.
                // Added for this version.
                PhaseVec upwind_viscosity = d[upwind_dir[0] - 1].viscosity;
                PhaseVec upwind_sat = d[upwind_dir[0] - 1].saturation;
                PhaseVec upwind_relperm = d[upwind_dir[0] - 1].relperm;
                PhaseVec downwind_mob(0.0);
                PhaseVec upwind_mob(0.0);
                double downwind_totmob = 0.0;
                double upwind_totmob = 0.0;
                for (int phase = 0; phase < numPhases; ++phase) {
                    downwind_mob[phase] = fluid_data_.cell_data.relperm[downwind_cell][phase]/upwind_viscosity[phase];
                    downwind_totmob += downwind_mob[phase];
                    upwind_mob[phase] = upwind_relperm[phase]/upwind_viscosity[phase];
                    upwind_totmob += upwind_mob[phase];
                }
                PhaseVec downwind_ff = downwind_mob;
                downwind_ff /= downwind_totmob;
                PhaseVec upwind_ff = upwind_mob;
                upwind_ff /= upwind_totmob;
                PhaseVec ff_diff = upwind_ff;
                ff_diff -= downwind_ff;
                for (int phase = 0; phase < numPhases; ++phase) {
                    if (std::fabs(ff_diff[phase]) > 1e-10) {
                        if (face_flux[face] != 0.0) {
                            double ff_deriv = ff_diff[phase]/(upwind_sat[phase] - fluid_data_.cell_data.saturation[downwind_cell][phase]);
                            // ASSERT(ff_deriv >= 0.0);
                            face_max_ff_deriv = std::max(face_max_ff_deriv, std::fabs(ff_deriv));
                        }
                    }
                }
            }
	    faces_max_ff_deriv[face] = face_max_ff_deriv;

            // Compute z change.
            CompVec change(0.0);
            for (int phase = 0; phase < numPhases; ++phase) {
                int upwind_ix = upwind_dir[phase] - 1; // Since process_face returns 1 or 2.
                CompVec z_in_phase = d[upwind_ix].phase_to_comp[phase];
                z_in_phase *= phase_flux[phase];
                change += z_in_phase;
            }
	    face_change[face] = change;
	}

	// Update output variables
#pragma omp parallel for
	for (int cell = 0; cell < num_cells; ++cell) {
	    const int num_local_faces = pgrid_->numCellFaces(cell);
	    for (int local = 0; local < num_local_faces; ++local) {
		int face = pgrid_->cellFace(cell,local);
		if (cell == pgrid_->faceCell(face, 0)) {
		    comp_change[cell] -= face_change[face];
		    if (face_flux[face] >= 0.0) {
			cell_outflux[cell] += std::fabs(face_flux[face]);
		    }
		} else if (cell == pgrid_->faceCell(face, 1)) {
		    comp_change[cell] += face_change[face];
		    if (face_flux[face] < 0.0) {
			cell_outflux[cell] += std::fabs(face_flux[face]);
		    }
		}
		cell_max_ff_deriv[cell] = std::max(cell_max_ff_deriv[cell],
						   faces_max_ff_deriv[face]);
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
            PhaseVec phase_flux = flow > 0.0 ? fl.fractional_flow : fluid_data_.fractional_flow[cell];
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
            Wells::WellReport::report()->massRate.push_back(change);
        }
    }


    // Arguments are:
    //  [in] cmob1, cmob2: cell mobilities for adjacent cells
    //  [in, destroyed] vstar: total velocity (equal for all phases) on input, modified on output
    //  [in] gf: gravity flux
    //  [in] np: number of phases
    //  [in, destroyed] rho: density on face
    //  [out] ix: upwind cells for each phase (1 or 2)
    static void 
    process_face(double *cmob1, double *cmob2, double *vstar, double gf, 
                 int np, double *rho, int *ix)
    {
        int i,j,k,a,b,c;
        double v, r, m;
        for (i=0; i<np; ++i) {
            
            /* ============= */
            /* same sign     */

            /* r = max(rho)*/
            j = 0; r = DBL_MIN;
            for(k=0; k<np; ++k) { if (rho[k] > r) { r = rho[k]; j = k; } }
            a = !(vstar[j]<0) && !(gf<0);
            b = !(vstar[j]>0) && !(gf>0);
            c = a || b;

            if ( !c ) {
                rho[j] = NAN;        
                v = vstar[j] - gf;
                if (v < 0) { ix[j] = 2; m = cmob2[j]; }
                else       { ix[j] = 1; m = cmob1[j]; }
                for (k=0; k<np; ++k) { vstar[k] -= m *(rho[k]-r)*gf;}
                continue;
            }


            /* ============= */
            /* opposite sign */

            /* r = min(rho)*/
            j = 0; r = DBL_MAX;
            for(k=0; k<np; ++k) { if (rho[k] < r) { r = rho[k]; j = k; } }
        
            if ( c ) {            
                rho[j] = NAN;       
                v      = vstar[j] + gf;
                if (v < 0) { ix[j] = 2; m = cmob2[j]; }
                else       { ix[j] = 1; m = cmob1[j]; }
                for (k=0; k<np; ++k) { vstar[k] -= m *(rho[k]-r)*gf;}
                continue;
            }
        }
    }




    static void
    phase_upwind_directions(int number_of_faces, int *face_cells, double *dflux, double *gflux, int np, 
                            double *cmobility, double *fdensity, int *ix)
    {
        int k,f;
        int c1, c2;
        double *cmob1, *cmob2;
        double *fden  = (double*)malloc(np * sizeof *fden);
        double *vstar = (double*)malloc(np * sizeof *vstar);
        double gf;
 
        for (f=0; f<number_of_faces; ++f) {
            gf = gflux[f];

            for (k=0; k<np; ++k) {
                vstar[k] = dflux[f];
                fden [k] = fdensity[np*f+k];
            }
            c1 = face_cells[2*f + 0];
            c2 = face_cells[2*f + 1];
            if (c1 != -1) {cmob1 = cmobility+np*c1;}
            else          {cmob1 = cmobility+np*c2;}
            if (c2 != -1) {cmob2 = cmobility+np*c2;}
            else          {cmob2 = cmobility+np*c1;}
            process_face(cmob1, cmob2, vstar, gf, np, fden, ix+np*f);
        
        }
        free(fden);
        free(vstar);
    }

};


} // namespace Opm


#endif // OPM_COMPONENTTRANSPORT_HEADER_INCLUDED
