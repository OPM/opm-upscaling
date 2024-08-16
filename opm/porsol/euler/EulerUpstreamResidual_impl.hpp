//===========================================================================
//
// File: EulerUpstreamResidual_impl.hpp
//
// Created: Thu May  6 11:22:04 2010
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//            Jostein R Natvig    <jostein.r.natvig@sintef.no>
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2010 Statoil ASA.

  This file is part of The Open Reservoir Simulator Project (OpenRS).

  OpenRS is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OpenRS is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OpenRS.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPENRS_EULERUPSTREAMRESIDUAL_IMPL_HEADER
#define OPENRS_EULERUPSTREAMRESIDUAL_IMPL_HEADER



#include <opm/core/utility/Average.hpp>
#include <opm/porsol/common/Matrix.hpp>

#ifdef USE_TBB
#include <tbb/parallel_for.h>
#endif

#include <iostream>


namespace Opm
{




    namespace EulerUpstreamResidualDetails
    {
	template <typename T, template <typename> class StoragePolicy, class OrderingPolicy>
	FullMatrix<T, OwnData, OrderingPolicy>
        arithAver(const FullMatrix<T, StoragePolicy, OrderingPolicy>& m1,
                  const FullMatrix<T, StoragePolicy, OrderingPolicy>& m2)
	{
	    return Opm::utils::arithmeticAverage<FullMatrix<T, StoragePolicy, OrderingPolicy>,
		FullMatrix<T, OwnData, OrderingPolicy> >(m1, m2);
	}

        template <class UpstreamSolver, class PressureSolution>
        struct UpdateForCell
        {
            typedef typename UpstreamSolver::Vector Vector;
            typedef typename UpstreamSolver::FIt FIt;
            typedef typename UpstreamSolver::RP::PermTensor PermTensor;
            typedef typename UpstreamSolver::RP::MutablePermTensor MutablePermTensor;

            const UpstreamSolver& s;
            const std::vector<double>& saturation;
            const Vector& gravity;
            const PressureSolution& pressure_sol;
            std::vector<double>& residual;

            UpdateForCell(const UpstreamSolver& solver,
                          const std::vector<double>& sat,
                          const Vector& grav,
                          const PressureSolution& psol,
                          std::vector<double>& res)
                : s(solver), saturation(sat), gravity(grav), pressure_sol(psol), residual(res)
            {
            }

            template <class CIt>
            void operator()(const CIt& c) const
            {
                // This is constant for the whole run.
                const double delta_rho = s.preservoir_properties_->densityDifference();
                int cell[2];
                double cell_sat[2];
                cell[0] = c->index();
                cell_sat[0] = saturation[cell[0]];

                // Loop over all cell faces.
                for (FIt f = c->facebegin(); f != c->faceend(); ++f) {
                    // Neighbour face, will be changed if on a periodic boundary.
                    FIt nbface = f;
                    double dS = 0.0;
                    // Compute cell[1], cell_sat[1]
                    if (f->boundary()) {
                        if (s.pboundary_->satCond(*f).isPeriodic()) {
                            nbface = s.bid_to_face_[s.pboundary_->getPeriodicPartner(f->boundaryId())];
                            assert(nbface != f);
                            cell[1] = nbface->cellIndex();
                            assert(cell[0] != cell[1]);
                            // Periodic faces will be visited twice, but only once
                            // should they contribute. We make sure that we skip the
                            // periodic faces half the time.
                            if (cell[0] > cell[1]) {
                                // We skip this face.
                                continue;
                            }
                            cell_sat[1] = saturation[cell[1]];
                        } else {
                            assert(s.pboundary_->satCond(*f).isDirichlet());
                            cell[1] = cell[0];
                            cell_sat[1] = s.pboundary_->satCond(*f).saturation();
                        }
                    } else {
                        cell[1] = f->neighbourCellIndex();
                        assert(cell[0] != cell[1]);
                        if (cell[0] > cell[1]) {
                            // We skip this face.
                            continue;
                        }
                        cell_sat[1] = saturation[cell[1]];
                    }

                    // Get some local properties.
                    const double loc_area = f->area();
                    const double loc_flux = pressure_sol.outflux(f);
                    const Vector loc_normal = f->normal();

                    // We will now try to establish the upstream directions for each
                    // phase. They may be the same, or different (due to gravity).
                    // Recall the equation for v_w (water phase velocity):
                    //   v_w  = lambda_w * (lambda_o + lambda_w)^{-1}
                    //          * (v + lambda_o * K * grad p_{cow} + lambda_o * K * (rho_w - rho_o) * g)
                    //             ^   ^^^^^^^^^^^^^^^^^^^^^^^^^^^   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                    //     viscous term       capillary term                    gravity term
                    //
                    // For the purpose of upstream weighting, we only consider the viscous and gravity terms.
                    // The question is, in which direction does v_w and v_o point? That is, what is the sign
                    // of v_w*loc_normal and v_o*loc_normal?
                    //
                    // For the case when the mobilities are scalar, the following analysis applies:
                    // The viscous contribution to v_w is loc_area*loc_normal*f_w*v == f_w*loc_flux.
                    // Then the phase fluxes become
                    //     flux_w = f_w*(loc_flux + loc_area*loc_normal*lambda_o*K*(rho_w - rho_o)*g)
                    //                              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                    //                                           := lambda_o*G (only scalar case)
                    //     flux_o = f_o*(loc_flux - lambda_w*G)
                    // In the above, we must decide where to evaluate K, and for this purpose (deciding
                    // upstream directions) we use a K averaged between the two cells.
                    // Since all mobilities and fractional flow functions are positive, the sign
                    // of one of these cases is trivial. If G >= 0, flux_w is in the same direction as
                    // loc_flux, if G <= 0, flux_o is in the same direction as loc_flux.
                    // The phase k for which flux_k and loc_flux are of same sign, is called the trivial
                    // phase in the code below.
                    //
                    // Assuming for the moment that G >=0, we know the direction of the water flux
                    // (same as loc_flux) and evaluate lambda_w in the upstream cell. Then we may use
                    // that lambda_w to evaluate flux_o using the above formula. Knowing flux_o, we know
                    // the direction of the oil flux, and can evaluate lambda_o in the corresponding
                    // upstream cell. Finally, we can use the equation for flux_w to compute that flux.
                    // The opposite case is similar.
                    //
                    // What about tensorial mobilities? In the following code, we make the assumption
                    // that the directions of vectors are not so changed by the multiplication with
                    // mobility tensors that upstream directions change. In other words, we let all
                    // the upstream logic stand as it is. This assumption may need to be revisited.
                    // A worse problem is that
                    // 1) we do not have v, just loc_area*loc_normal*v,
                    // 2) we cannot define G, since the lambdas do not commute with the dot product.

                    typedef typename UpstreamSolver::RP::Mobility Mob;
                    using Opm::utils::arithmeticAverage;
                    // Doing arithmetic averages. Should we consider harmonic or geometric instead?
                    const MutablePermTensor aver_perm
                        = arithAver(s.preservoir_properties_->permeability(cell[0]),
                                    s.preservoir_properties_->permeability(cell[1]));
                    // Computing the raw gravity influence vector = (rho_w - rho_o)Kg
                    Vector grav_influence = prod(aver_perm, gravity);
                    grav_influence *= delta_rho;
                    // Computing G. Note that we do not multiply with the mobility,
                    // so this G is wrong in case of anisotropic relperm.
                    const double G = s.method_gravity_ ?
                        loc_area*inner(loc_normal, grav_influence) 
                        : 0.0;
                    const int triv_phase = G >= 0.0 ? 0 : 1;
                    const int ups_cell = loc_flux >= 0.0 ? 0 : 1;
                    // Compute mobility of the trivial phase.
                    Mob m_ups[2];
                    s.preservoir_properties_->phaseMobility(triv_phase, cell[ups_cell],
                                                          cell_sat[ups_cell], m_ups[triv_phase].mob);
                    // Compute gravity flow of the nontrivial phase.
                    const double sign_G[2] = { -1.0, 1.0 };
                    double grav_flux_nontriv = sign_G[triv_phase]*loc_area
                        *inner(loc_normal, m_ups[triv_phase].multiply(grav_influence));
                    // Find flow direction of nontrivial phase.
                    const int ups_cell_nontriv = (loc_flux + grav_flux_nontriv >= 0.0) ? 0 : 1;
                    const int nontriv_phase = (triv_phase + 1) % 2;
                    s.preservoir_properties_->phaseMobility(nontriv_phase, cell[ups_cell_nontriv],
                                                          cell_sat[ups_cell_nontriv], m_ups[nontriv_phase].mob);
                    // Now we have the upstream phase mobilities in m_ups[].
                    Mob m_tot;
                    m_tot.setToSum(m_ups[0], m_ups[1]);
                    Mob m_totinv;
                    m_totinv.setToInverse(m_tot);


                    const double aver_sat
                        = Opm::utils::arithmeticAverage<double, double>(cell_sat[0], cell_sat[1]);

                    Mob m1c0, m1c1, m2c0, m2c1;
                    s.preservoir_properties_->phaseMobility(0, cell[0], aver_sat, m1c0.mob);
                    s.preservoir_properties_->phaseMobility(0, cell[1], aver_sat, m1c1.mob);
                    s.preservoir_properties_->phaseMobility(1, cell[0], aver_sat, m2c0.mob);
                    s.preservoir_properties_->phaseMobility(1, cell[1], aver_sat, m2c1.mob);
                    Mob m_aver[2];
                    m_aver[0].setToAverage(m1c0, m1c1);
                    m_aver[1].setToAverage(m2c0, m2c1);
                    Mob m_aver_tot;
                    m_aver_tot.setToSum(m_aver[0], m_aver[1]);
                    Mob m_aver_totinv;
                    m_aver_totinv.setToInverse(m_aver_tot);

                    // Viscous (pressure driven) term.
                    if (s.method_viscous_) {
                        // v is not correct for anisotropic relperm.
                        Vector v(loc_normal);
                        v *= loc_flux;
                        const double visc_change = inner(loc_normal, m_ups[0].multiply(m_totinv.multiply(v)));
                        // 		    const double visc_change = (m_ups[0].mob/(m_ups[1].mob + m_ups[0].mob))*loc_flux;
                        // 		    std::cout << "New: " << visc_change_2 << "   old: " << visc_change << '\n';
                        dS += visc_change;
                    }

                    // Gravity term.
                    if (s.method_gravity_) {
                        if (cell[0] != cell[1]) {
                            // We only add gravity flux on internal or periodic faces.
                            const double grav_change = loc_area
                                *inner(loc_normal, m_ups[0].multiply(m_totinv.multiply(m_ups[1].multiply(grav_influence))));
                            // const double grav_change = (lambda_one*lambda_two/(lambda_two+lambda_one))*G;
                            // const double grav_change = (lambda_one*lambda_two/(lambda_two+lambda_one))*loc_gravity_flux;
                            dS += grav_change;
                        }
                    }

                    // Capillary term.
                    if (s.method_capillary_) {
                        // J(s_w) = \frac{p_c(s_w)\sqrt{k/\phi}}{\sigma \cos\theta}
                        // p_c = \frac{J \sigma \cos\theta}{\sqrt{k/\phi}}
                        Vector cap_influence = prod(aver_perm, s.estimateCapPressureGradient(f, nbface));
                        const double cap_change = loc_area
			    *inner(loc_normal, m_aver[0].multiply(m_aver_totinv.multiply(m_aver[1].multiply(cap_influence))));
                        dS += cap_change;
                    }

                    // Modify saturation.
                    if (cell[0] != cell[1]){
                        residual[cell[0]] -= dS;
                        residual[cell[1]] += dS;
                    } else {
                        assert(cell[0] == cell[1]);
                        residual[cell[0]] -= dS;
                    }
                }
                // Source term.
                double rate = s.pinjection_rates_->element(cell[0]);
                if (rate < 0.0) {
                    // For anisotropic relperm, fractionalFlow does not really make sense
                    // as a scalar
                    rate *= s.preservoir_properties_->fractionalFlow(cell[0], cell_sat[0]);
                }
                residual[cell[0]] += rate;
            }
        };

        template <typename Iter>
        struct IndirectRange
        {
            typedef Iter Iterator;
            explicit IndirectRange(const std::vector<Iter>& iters)
                : iters_(iters), beg_(0), end_(iters_.size() - 1)
            {
                assert(iters_.size() >= 2);
            }
#ifdef USE_TBB
            IndirectRange(IndirectRange& r, tbb::split)
                : iters_(r.iters_)
            {
                int m = (r.beg_ + r.end_)/2;
                beg_ = m;
                end_ = r.end_;
                r.end_ = m;
            }
#endif
            bool empty() const
            {
                return beg_ == end_;
            }
            bool is_divisible() const
            {
                return end_ - beg_ > 1;
            }
            Iter begin() const
            {
                return iters_[beg_];
            }
            Iter end() const
            {
                return iters_[end_];
            }
        private:
            const std::vector<Iter>& iters_;
            int beg_;
            int end_;
        };

        template <class Updater>
        struct UpdateLoopBody
        {
            explicit UpdateLoopBody(const Updater& upd)
                : updater(upd)
            {
            }
            const Updater& updater;
            template <class Range>
            void operator()(const Range& r) const
            {
                typename Range::Iterator c = r.begin();
                typename Range::Iterator cend = r.end();
                for (; c != cend; ++c) {
                    updater(c);
                }
            }
        };

    } // namespace EulerUpstreamResidualDetails


    // --------- Member functions -----------



    template <class GI, class RP, class BC>
    inline EulerUpstreamResidual<GI, RP, BC>::EulerUpstreamResidual()
	: pgrid_(0),
	  preservoir_properties_(0),
	  pboundary_(0)
    {
    }


    template <class GI, class RP, class BC>
    inline EulerUpstreamResidual<GI, RP, BC>::EulerUpstreamResidual(const GI& g, const RP& r, const BC& b)
	: pgrid_(&g),
	  preservoir_properties_(&r),
	  pboundary_(&b)
    {
        initFinal();
    }



    template <class GI, class RP, class BC>
    inline void EulerUpstreamResidual<GI, RP, BC>::initObj(const GI& g, const RP& r, const BC& b)
    {
	pgrid_ = &g;
        preservoir_properties_ = &r;
        pboundary_ = &b;
        initFinal();
    }




    template <class GI, class RP, class BC>
    inline void EulerUpstreamResidual<GI, RP, BC>::initFinal()
    {
	// Build bid_to_face_ mapping for handling periodic conditions.
	int maxbid = 0;
	for (typename GI::CellIterator c = pgrid_->cellbegin(); c != pgrid_->cellend(); ++c) {
	    for (typename GI::CellIterator::FaceIterator f = c->facebegin(); f != c->faceend(); ++f) {
		int bid = f->boundaryId();
		maxbid = std::max(maxbid, bid);
	    }
	}
	bid_to_face_.clear();
	bid_to_face_.resize(maxbid + 1);
	for (typename GI::CellIterator c = pgrid_->cellbegin(); c != pgrid_->cellend(); ++c) {
	    for (typename GI::CellIterator::FaceIterator f = c->facebegin(); f != c->faceend(); ++f) {
		if (f->boundary() && pboundary_->satCond(*f).isPeriodic()) {
		    bid_to_face_[f->boundaryId()] = f;
		}
	    }
	}

        // Build cell_iters_.
        const int num_cells_per_iter = std::min(50, pgrid_->numberOfCells());
        int counter = 0;
	for (typename GI::CellIterator c = pgrid_->cellbegin(); c != pgrid_->cellend(); ++c, ++counter) {
            if (counter % num_cells_per_iter == 0) {
                cell_iters_.push_back(c);
            }
        }
        cell_iters_.push_back(pgrid_->cellend());
    }



    template <class GI, class RP, class BC>
    inline const GI& EulerUpstreamResidual<GI, RP, BC>::grid() const
    {
        return *pgrid_;
    }



    template <class GI, class RP, class BC>
    inline const RP& EulerUpstreamResidual<GI, RP, BC>::reservoirProperties() const
    {
        return *preservoir_properties_;
    }


    template <class GI, class RP, class BC>
    inline const BC& EulerUpstreamResidual<GI, RP, BC>::boundaryConditions() const
    {
        return *pboundary_;
    }



    template <class GI, class RP, class BC>
    inline void EulerUpstreamResidual<GI, RP, BC>::computeCapPressures(const std::vector<double>& saturation) const
    {
	int num_cells = saturation.size();
	cap_pressures_.resize(num_cells);
	for (int cell = 0; cell < num_cells; ++cell) {
	    cap_pressures_[cell] = preservoir_properties_->capillaryPressure(cell, saturation[cell]);
	}
    }




    template <class GI, class RP, class BC>
    template <class PressureSolution>
    inline void EulerUpstreamResidual<GI, RP, BC>::
    computeResidual(const std::vector<double>& saturation,
                    const typename GI::Vector& gravity,
                    const PressureSolution& pressure_sol,
                    const Opm::SparseVector<double>& injection_rates,
                    const bool method_viscous,
                    const bool method_gravity,
                    const bool method_capillary,
                    std::vector<double>& residual) const
    {
	// Make sure sat_change is zero, and has the right size.
	residual.clear();
	residual.resize(saturation.size(), 0.0);

        pinjection_rates_ = &injection_rates;
        method_viscous_ = method_viscous;
        method_gravity_ = method_gravity;
        method_capillary_ = method_capillary;

	// For every face, we will modify residual for adjacent cells.
	// We loop over every cell and intersection, and modify only if
	// this cell has lower index than the neighbour, or we are on the boundary.
        typedef EulerUpstreamResidualDetails::UpdateForCell<EulerUpstreamResidual<GI,RP,BC>, PressureSolution> CellUpdater;
        CellUpdater update_cell(*this, saturation, gravity, pressure_sol, residual);
        EulerUpstreamResidualDetails::UpdateLoopBody<CellUpdater> body(update_cell);
        EulerUpstreamResidualDetails::IndirectRange<CIt> r(cell_iters_);
#ifdef USE_TBB
        tbb::parallel_for(r, body);
#else
        body(r);
#endif
    }




    template <class GI, class RP, class BC>
    inline typename GI::Vector
    EulerUpstreamResidual<GI, RP, BC>::
    estimateCapPressureGradient(const FIt& f, const FIt& nbf) const
    {
	// At nonperiodic boundaries, we return a zero gradient.
	// That is (sort of) a trivial Neumann (noflow) condition for the capillary pressure.
	if (f->boundary() && !pboundary_->satCond(*f).isPeriodic()) {
	    return Vector(0.0);
	}
	// Find neighbouring cell and face: nbc and nbf.
	// If we are not on a periodic boundary, nbf is of course equal to f.
	auto c = f->cell();
	auto nb = f->boundary() ? (f == nbf ? c : nbf->cell()) : f->neighbourCell();

	// Estimate the gradient like a finite difference between
	// cell centers, except that in order to handle periodic
	// conditions we pass through the face centroid(s).
	auto cell_c = c.centroid();
	auto nb_c = nb.centroid();
	auto f_c = f->centroid();
	auto nbf_c = nbf->centroid();
	double d0 = (cell_c - f_c).two_norm();
	double d1 = (nb_c - nbf_c).two_norm();
	int cell = c.index();
	int nbcell = nb.index();
	double cp0 = cap_pressures_[cell];
	double cp1 = cap_pressures_[nbcell];
	double val = (cp1 - cp0)/(d0 + d1);
	auto res = nb_c - nbf_c + f_c - cell_c;
	res /= res.two_norm();
	res *= val;
	return res;
    }


} // namespace Opm


#endif // OPENRS_EULERUPSTREAMRESIDUAL_IMPL_HEADER
