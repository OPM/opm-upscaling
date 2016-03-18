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

#ifndef OPM_FLUIDMATRIXINTERACTIONBLACKOIL_HEADER_INCLUDED
#define OPM_FLUIDMATRIXINTERACTIONBLACKOIL_HEADER_INCLUDED

#include <opm/core/utility/UniformTableLinear.hpp>
#include <opm/core/utility/buildUniformMonotoneTable.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/Parser/ParseContext.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/SwofTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/SgofTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/TableManager.hpp>
#include "BlackoilDefs.hpp"
#include <iostream>
#include <stdexcept>

namespace Opm
{

// Forward declaration needed by associated parameters class.
template <class ScalarT, class ParamsT>
class FluidMatrixInteractionBlackoil;

template <class ScalarT>
class FluidMatrixInteractionBlackoilParams
{
public:
    typedef ScalarT Scalar;
    void init(Opm::DeckConstPtr deck)
    {
        ParseContext parseContext;
        EclipseState eclipseState(deck , parseContext);
        // Extract input data.
        const auto& tables    = eclipseState.getTableManager();
        const auto& swofTable = tables->getSwofTables().getTable<SwofTable>(0);
        const auto& sgofTable = tables->getSgofTables().getTable<SgofTable>(0);

        std::vector<double> sw = swofTable.getSwColumn().vectorCopy();
        std::vector<double> krw = swofTable.getKrwColumn().vectorCopy();
        std::vector<double> krow = swofTable.getKrowColumn().vectorCopy();
        std::vector<double> pcow = swofTable.getPcowColumn().vectorCopy();
        std::vector<double> sg = sgofTable.getSgColumn().vectorCopy();
        std::vector<double> krg = sgofTable.getKrgColumn().vectorCopy();
        std::vector<double> krog = sgofTable.getKrogColumn().vectorCopy();
        std::vector<double> pcog = sgofTable.getPcogColumn().vectorCopy();

        // Create tables for krw, krow, krg and krog.
        int samples = 200;
        buildUniformMonotoneTable(sw, krw,  samples, krw_);
        buildUniformMonotoneTable(sw, krow, samples, krow_);
        buildUniformMonotoneTable(sg, krg,  samples, krg_);
        buildUniformMonotoneTable(sg, krog, samples, krog_);
        krocw_ = krow[0]; // At connate water -> ecl. SWOF

        // Create tables for pcow and pcog.
        buildUniformMonotoneTable(sw, pcow, samples, pcow_);
        buildUniformMonotoneTable(sg, pcog, samples, pcog_);
    }

private:
    template <class S, class P>
    friend class FluidMatrixInteractionBlackoil;

    Opm::utils::UniformTableLinear<Scalar> krw_;
    Opm::utils::UniformTableLinear<Scalar> krow_;
    Opm::utils::UniformTableLinear<Scalar> pcow_;
    Opm::utils::UniformTableLinear<Scalar> krg_;
    Opm::utils::UniformTableLinear<Scalar> krog_;
    Opm::utils::UniformTableLinear<Scalar> pcog_;
    Scalar krocw_; // = krow_(s_wc)
};


/*!
 * \ingroup material
 *
 * \brief Capillary pressures and relative permeabilities for a black oil system.
 */
template <class ScalarT, class ParamsT = FluidMatrixInteractionBlackoilParams<ScalarT> >
class FluidMatrixInteractionBlackoil : public BlackoilDefs
{
public:
    typedef ParamsT Params;
    typedef typename Params::Scalar Scalar;

    /*!
     * \brief The linear capillary pressure-saturation curve.
     *
     * This material law is linear:
     * \f[
     p_C = (1 - \overline{S}_w) (p_{C,max} - p_{C,entry}) + p_{C,entry}
     \f]
     *
     * \param Swe Effective saturation of of the wetting phase \f$\overline{S}_w\f$
     */
    template <class pcContainerT, class SatContainerT>
    static void pC(pcContainerT &pc,
                   const Params &params, 
                   const SatContainerT &saturations,
                   Scalar /*temperature*/)
    {
        Scalar sw = saturations[Aqua];
        Scalar sg = saturations[Vapour];
        pc[Liquid] = 0.0;
        pc[Aqua] = params.pcow_(sw);
        pc[Vapour] = params.pcog_(sg);
    }

    /*!
     * \brief The saturation-capillary pressure curve.
     *
     * This is the inverse of the capillary pressure-saturation curve:
     * \f[
     S_w = 1 - \frac{p_C - p_{C,entry}}{p_{C,max} - p_{C,entry}}
     \f]
     *
     * \param pC Capillary pressure \f$\p_C\f$
     * \return The effective saturaion of the wetting phase \f$\overline{S}_w\f$
     */
    template <class SatContainerT, class pcContainerT>
    static void S(SatContainerT &saturations,
                  const Params &params, 
                  const pcContainerT &pc,
                  Scalar /*temperature*/)
    {
        std::cerr << "FluidMatrixInteractionBlackoil::S() is not implemented yet\n";
        throw std::logic_error("Not implemented");
    }


    /*!
     * \brief The relative permeability of all phases.
     */
    template <class krContainerT, class SatContainerT>
    static void kr(krContainerT &kr,
                   const Params &params, 
                   const SatContainerT &saturations,
                   Scalar /*temperature*/)
    {
        // Stone-II relative permeability model.
        Scalar sw = saturations[Aqua];
        Scalar sg = saturations[Vapour];
        Scalar krw = params.krw_(sw);
        Scalar krg = params.krg_(sg);
        Scalar krow = params.krow_(sw);
        Scalar krog = params.krog_(sg);
        Scalar krocw = params.krocw_;
        kr[Aqua] = krw;
        kr[Vapour] = krg;
        kr[Liquid] = krocw*((krow/krocw + krw)*(krog/krocw + krg) - krw - krg);
        if (kr[Liquid] < 0.0) {
            kr[Liquid] = 0.0;
        }
    }


    /*!
     * \brief The saturation derivatives of relative permeability of all phases.
     */
    template <class krContainerT, class SatContainerT>
    static void dkr(krContainerT &dkr,
                    const Params &params, 
                    const SatContainerT &saturations,
                    Scalar /*temperature*/)
    {
        for (int p1 = 0; p1 < numPhases; ++p1) {
            for (int p2 = 0; p2 < numPhases; ++p2) {
                dkr[p1][p2] = 0.0;
            }
        }
        // Stone-II relative permeability model.
        Scalar sw = saturations[Aqua];
        Scalar sg = saturations[Vapour];
        Scalar krw = params.krw_(sw);
        Scalar dkrww = params.krw_.derivative(sw);
        Scalar krg = params.krg_(sg);
        Scalar dkrgg = params.krg_.derivative(sg);
        Scalar krow = params.krow_(sw);
        Scalar dkrow = params.krow_.derivative(sw);
        Scalar krog = params.krog_(sg);
        Scalar dkrog = params.krog_.derivative(sg);
        Scalar krocw = params.krocw_;
        dkr[Aqua][Aqua] = dkrww;
        dkr[Vapour][Vapour] = dkrgg;
        dkr[Liquid][Aqua] = krocw*((dkrow/krocw + dkrww)*(krog/krocw + krg) - dkrww);
        dkr[Liquid][Vapour] = krocw*((krow/krocw + krw)*(dkrog/krocw + dkrgg) - dkrgg);
    }
};

} // namespace Opm




#endif // OPM_FLUIDMATRIXINTERACTIONBLACKOIL_HEADER_INCLUDED
