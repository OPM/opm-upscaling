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

#ifndef OPM_TABULATEDPHASE_HEADER_INCLUDED
#define OPM_TABULATEDPHASE_HEADER_INCLUDED


#include <dune/porsol/common/UniformTableLinear.hpp>

namespace Opm
{


    class TabulatedPhase;

    class TabulatedPhaseParams
    {
    public:
        typedef double Scalar;
    private:
        friend class TabulatedPhase;
        UniformTableLinear<double> kr_s;       // kr(s)
        UniformTableLinear<double> dkr_ds;     // dkr/ds(s)
        UniformTableLinear<double> pc_s;       // pc(s)
        UniformTableLinear<double> dpc_ds;     // dpc/ds(s)
        UniformTableLinear<double> s_pc;       // s(pc)
        UniformTableLinear<double> ds_dpc;     // ds/dpc(pc)
    };


    /*!\ingroup material
     *
     * \brief Implementation of tabulated phase behaviour.
     */
    class TabulatedPhase
    {
    public:
        typedef TabulatedPhaseParams Params;
        typedef Params::Scalar Scalar;

        /*!
         * \brief The relative permeability function.
         *
         * \param s The mobile saturation of the phase under consideration.
         */
        static Scalar kr(const Params& params, Scalar s)
        {
            return params.kr(s);
        }

        /*!
         * \brief The derivative of the relative permeability.
         *
         * \param s The mobile saturation of the phase under consideration.
         */
        static Scalar dkr_ds(const Params& params, Scalar s)
        {
            return params.dkr_ds(s);
        }

        /*!
         * \brief The capillary pressure-saturation curve.
         *
         * \param s   Effective saturation of the phase under consideration.
         */
        static Scalar pc(const Params& params, Scalar s)
        {
            return params.pc(s);
        }

        /*!
         * \brief Returns the partial derivative of the capillary
         *        pressure to the effective saturation.
         *
         * This is equivalent to
         * \f[
         \frac{\partial p_c}{\partial \overline{s}} =
         -\frac{p_e}{\alpha} \overline{s}^{-1/\alpha - 1}
         \f]
        */
        static Scalar dpc_ds(const Params& params, Scalar s)
        {
            return params.dpc_ds(s);
        }

        /*!
         * \brief The saturation-capillary pressure curve.
         *
         * This is the inverse of the capillary pressure-saturation curve:
         * \f[
         \overline{s} = (\frac{p_c}{p_e})^{-\alpha}
         \f]
         *
         * \param pc Capillary pressure \f$p_c\f$
         * \return The effective saturaion of the wetting phase \f$\overline{s}\f$
         */
        static Scalar s(const Params& params, Scalar pc)
        {
            return params.s(pc);
        }

        /*!
         * \brief Returns the partial derivative of the effective
         *        saturation to the capillary pressure.
         */
        static Scalar ds_dpc(const Params& params, Scalar pc)
        {
            return params.ds_dpc(pc);
        }

    };

} // namespace Opm

#endif // OPM_TABULATEDPHASE_HEADER_INCLUDED
