/*****************************************************************************
 *   Copyright (C) 2008-2010 by Melanie Darcis                               *
 *   Copyright (C) 2009-2010 by Andreas Lauser                               *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief A fluid system with water and gas as phases and brine and CO2
 *        as components.
 */
#ifndef OPM_BRINE_CO2_SYSTEM_HH
#define OPM_BRINE_CO2_SYSTEM_HH

#include <dune/porsol/blackoil/co2fluid/opm/material/idealgas.hh>

#include <dune/porsol/blackoil/co2fluid/opm/material/components/brine.hh>
#include <dune/porsol/blackoil/co2fluid/opm/material/components/co2.hh>
#include <dune/porsol/blackoil/co2fluid/opm/material/components/simpleco2.hh>
#include <dune/porsol/blackoil/co2fluid/opm/material/components/simpleh2o.hh>
#include <dune/porsol/blackoil/co2fluid/opm/material/components/tabulatedcomponent.hh>
#include <dune/porsol/blackoil/co2fluid/opm/material/binarycoefficients/brine_co2.hh>

#include <dune/porsol/blackoil/co2fluid/opm/material/settablephase.hh>

namespace Opm
{

/*!
 * \brief A compositional fluid with water and molecular nitrogen as
 *        components in both, the liquid and the gas phase.
 */
template <class CO2Tables, bool verbose=true>
class Brine_CO2_System
{
    typedef Brine_CO2_System<CO2Tables> ThisType;
    typedef double Scalar;

    typedef Opm::IdealGas<Scalar> IdealGas;
    typedef Opm::BinaryCoeff::Brine_CO2<Scalar, CO2Tables, verbose> Brine_CO2;
    typedef Opm::SettablePhase<Scalar, ThisType> SettablePhase;

    typedef Opm::H2O<Scalar> H2O_IAPWS;
    typedef Opm::Brine<Scalar, H2O_IAPWS> Brine_IAPWS;
    typedef Opm::TabulatedComponent<Scalar, H2O_IAPWS, verbose> H2O_Tabulated;
    typedef Opm::TabulatedComponent<Scalar, Brine_IAPWS, verbose> Brine_Tabulated;

public:
    typedef H2O_Tabulated H2O;
    typedef Brine_Tabulated Brine;
    //typedef H2O_IAPWS H2O;
    //typedef Opm::SimpleH2O<Scalar>                   H2O;
    typedef Opm::CO2<Scalar, CO2Tables>              CO2;

    static const int numComponents = 2;
    static const int numPhases = 2;

    static const int lPhaseIdx = 0; // index of the liquid phase
    static const int gPhaseIdx = 1; // index of the gas phase

    static const int wPhaseIdx = lPhaseIdx; // index of the wetting phase
    static const int nPhaseIdx = gPhaseIdx; // index of the non-wetting phase

    static const int BrineIdx = 0;
    static const int CO2Idx = 1;

    static void init()
    {
        std::cout << "Initializing tables for the H2O fluid properties.\n";
        H2O_Tabulated::init(273.15, 623.15, 100,
                            -10, 40e6, 200);
        std::cout << "Initializing tables for the brine fluid properties.\n";
        // set the salinity of brine to the one used by the CO2 tables
        Brine_IAPWS::salinity = CO2Tables::brineSalinity;
        Brine_Tabulated::init(273.15, 623.15, 100,
                              -10, 40e6, 200);
    }

    static void init(Scalar startTemp, Scalar endTemp, int tempSteps,
    		Scalar startPressure, Scalar endPressure, int pressureSteps)
    {
    	std::cout << "Using problem-specific bounds for tables";
        std::cout << "Initializing tables for the H2O fluid properties.\n";
        H2O_Tabulated::init(startTemp,endTemp, tempSteps,
        					startPressure, endPressure, pressureSteps);
        std::cout << "Initializing tables for the brine fluid properties.\n";
        // set the salinity of brine to the one used by the CO2 tables
        Brine_IAPWS::salinity = CO2Tables::brineSalinity;
        Brine_Tabulated::init(startTemp,endTemp, tempSteps,
        					startPressure, endPressure, pressureSteps);
    }

    /*!
     * \brief Return the human readable name of a component
     */
    static const char *componentName(int compIdx)
    {
        switch (compIdx) {
        case BrineIdx: return Brine::name();
        case CO2Idx: return CO2::name();
        };
        OPM_THROW(Opm::InvalidStateException, "Invalid component index " << compIdx);
    }

    /*!
     * \brief Return the molar mass of a component in [kg/mol].
     */
    static Scalar molarMass(int compIdx)
    {
        switch (compIdx) {
        case BrineIdx: return Brine::molarMass();
        case CO2Idx: return CO2::molarMass();
        };
        OPM_THROW(Opm::InvalidStateException, "Invalid component index " << compIdx);
    }

    /*!
     * \brief Given all mole fractions in a phase, return the phase
     *        density [kg/m^3].
     */
    template <class FluidState>
    static Scalar phaseDensity(int phaseIdx,
                               Scalar temperature,
                               Scalar pressure,
                               const FluidState &fluidState)
    {
        switch (phaseIdx) {
        case lPhaseIdx:
        {
            return liquidDensity_(temperature,
                                  pressure,
                                  fluidState.moleFraction(lPhaseIdx, BrineIdx),
                                  fluidState.moleFraction(lPhaseIdx, CO2Idx));
        }
        case gPhaseIdx:
        {
            return gasDensity_(temperature,
                               pressure,
                               fluidState.moleFraction(gPhaseIdx, BrineIdx),
                               fluidState.moleFraction(gPhaseIdx, CO2Idx));
        };
        }
        OPM_THROW(Opm::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    /*!
     * \brief Return the viscosity of a phase.
     */
    template <class FluidState>
    static Scalar phaseViscosity(int phaseIdx,
                                 Scalar temperature,
                                 Scalar pressure,
                                 const FluidState &fluidState)
    {
        if (phaseIdx == lPhaseIdx) {
            // assume pure brine for the liquid phase. TODO: viscosity
            // of mixture
            return Brine::liquidViscosity(temperature, pressure);
        }
        else {
            return CO2::gasViscosity(temperature, pressure);
        }
    }

    /*!
     * \brief Returns the solubility of a component in a
     *        phase.
     *
     * beware, the function name "activityCoeff" is misleading, as
     * this only determines how much is soluted in the phase!
     */
    template <class FluidState>
    static Scalar activityCoeff(int phaseIdx,
                                int compIdx,
                                Scalar temperature,
                                Scalar pressure,
                                const FluidState &state)
    {
        if (phaseIdx == gPhaseIdx)
        {
            if (compIdx == BrineIdx)
                return Brine_CO2::fugacityCoefficientH2O(temperature, pressure) * pressure;
            else if (compIdx == CO2Idx)
                return Brine_CO2::fugacityCoefficientCO2(temperature, pressure) * pressure;

            OPM_THROW(Opm::InvalidStateException, "Invalid component index " << compIdx);
        }
        Scalar xlCO2(0.), ygH2O(0.);
        Brine_CO2::calculateMoleFractions(temperature, pressure, Brine_IAPWS::salinity, -1, xlCO2, ygH2O);

        switch (compIdx) {
        case BrineIdx:
            {
                return (ygH2O / (1-0.0329216 - xlCO2)) * pressure;
            }
        case CO2Idx:
            {
                // sol = xgCO2 / xlCO2 * p_n is the activity if one assumes:
                //
                Scalar sol = ((1-0.0329216 - ygH2O) / xlCO2);
                return sol*pressure;
            }
        };
        OPM_THROW(Opm::InvalidStateException, "Invalid component index " << compIdx);
    }

    /*!
     * \brief Assuming the composition of a single phase and the
     *        pressure of all phases is known or that all phases are
     *        present, compute the thermodynamic equilibrium from the
     *        temperature and phase pressures. If the known phase
     *        index
     *
     */
    template <class FluidState>
    static void computeEquilibrium(FluidState &fluidState,
                                   int knownPhaseIdx = -1)
    {
        const Scalar T = fluidState.temperature();
        const Scalar pg = fluidState.pressure(gPhaseIdx);
        const Scalar pl = fluidState.pressure(lPhaseIdx);
        Scalar xlCO2, ygH2O, xlBrine, ygCO2;
        Scalar xlCO2_max, ygH2O_max;

        if (knownPhaseIdx < 0)
        {
            // we only have all phase pressures and temperature and
            // know that all phases are present
            SettablePhase liquid;
            xlCO2 = 0.0;
            ygH2O = 0.0;
            Brine_CO2::calculateMoleFractions(T, pg, Brine_IAPWS::salinity, knownPhaseIdx, xlCO2, ygH2O);
            Scalar xlBrine = 1 - xlCO2;
            Scalar rhol = liquidDensity_(T, pl, xlBrine, xlCO2);
            liquid.moleFrac_[BrineIdx] = xlBrine;
            liquid.moleFrac_[CO2Idx] = xlCO2;
            liquid.pressure_ = pl;
            liquid.density_ = rhol;
            liquid.xToX(); // update mass fractions from mole fractions
            fluidState.assignPhase(lPhaseIdx, liquid);

            SettablePhase gas;
            ygCO2 = 1 - ygH2O;
            Scalar rhog = gasDensity_(T, pg, ygH2O, ygCO2);
            gas.moleFrac_[BrineIdx] = ygH2O;
            gas.moleFrac_[CO2Idx] = ygCO2;
            gas.pressure_ = pg;
            gas.density_ = rhog;
            gas.xToX(); // update mass fractions from mole fractions
            fluidState.assignPhase(gPhaseIdx, gas);
        }
        else if (knownPhaseIdx == lPhaseIdx) {
            // calculate composition of liquid phase
            SettablePhase gas;
            xlCO2 = fluidState.moleFraction(lPhaseIdx, CO2Idx) ;
            ygH2O = 0.0;
            Brine_CO2::calculateMoleFractions(T, pg, Brine_IAPWS::salinity, -1, xlCO2_max, ygH2O_max);
            Brine_CO2::calculateMoleFractions(T, pg, Brine_IAPWS::salinity, knownPhaseIdx, xlCO2, ygH2O);
            Scalar ygCO2 = 1 - ygH2O;
            Scalar ygH2O = ygH2O_max;
            Scalar rhog = gasDensity_(T, pg, ygH2O, ygCO2);
            gas.moleFrac_[BrineIdx] = ygH2O;
            gas.moleFrac_[CO2Idx] = ygCO2;
            gas.pressure_ = pg;
            gas.density_ = rhog;
            gas.xToX(); // update mass fractions from mole fractions
            fluidState.assignPhase(gPhaseIdx, gas);
            // the composition of the liquid phase is given
        }

        else if (knownPhaseIdx == gPhaseIdx) {
            // calculate composition of liquid phase
            SettablePhase liquid;
            ygH2O = fluidState.moleFraction(gPhaseIdx, BrineIdx);
            xlCO2 = 0.0;
            Brine_CO2::calculateMoleFractions(T, pg, Brine_IAPWS::salinity, -1, xlCO2_max, ygH2O_max);
            Brine_CO2::calculateMoleFractions(T, pg, Brine_IAPWS::salinity, knownPhaseIdx, xlCO2, ygH2O);
            xlBrine = 1 - xlCO2;
            xlCO2 = xlCO2_max;
            Scalar rhol = liquidDensity_(T, pl, xlBrine, xlCO2);
            liquid.moleFrac_[BrineIdx] = xlBrine;
            liquid.moleFrac_[CO2Idx] = xlCO2;
            liquid.pressure_ = pl;
            liquid.density_ = rhol;
            liquid.xToX(); // update mass fractions from mole fractions
            fluidState.assignPhase(lPhaseIdx, liquid);
        }
    }

    /*!
     * \brief Given the phase compositions, return the binary
     *        diffusion coefficent of two components in a phase.
     */
    template <class FluidState>
    static Scalar diffCoeff(int phaseIdx,
                            int compIIdx,
                            int compJIdx,
                            Scalar temperature,
                            Scalar pressure,
                            const FluidState &fluidState)
    {
        if (compIIdx > compJIdx)
            std::swap(compIIdx, compJIdx);

        switch (phaseIdx) {
        case lPhaseIdx:
            switch (compIIdx) {
            case BrineIdx:
                switch (compJIdx) {
                case CO2Idx: return Brine_CO2::liquidDiffCoeff(temperature, pressure);
                default:
                    OPM_THROW(Opm::InvalidStateException,
                               "Binary diffusion coefficients of trace "
                               "substances in liquid phase is undefined!\n");
                }
            }
        case gPhaseIdx:
            switch (compIIdx) {
            case BrineIdx:
                switch (compJIdx) {
                case CO2Idx: return Brine_CO2::gasDiffCoeff(temperature, pressure);
                }
            }
        }

        OPM_THROW(Opm::InvalidStateException,
                   "Binary diffusion coefficient of components "
                   << compIIdx << " and " << compJIdx
                   << " in phase " << phaseIdx << " is undefined!\n");
    };

    /*!
     * \brief Given the phase composition, return the specific
     *        phase enthalpy [J/kg].
     */
    template <class FluidState>
    static Scalar phaseEnthalpy(int phaseIdx,
                           Scalar temperature,
                           Scalar pressure,
                           const FluidState &fluidState)
    {
        if (phaseIdx == lPhaseIdx) {
            return liquidEnthalpyBrineCO2_(temperature,
                                           pressure,
                                           Brine_IAPWS::salinity,
                                           fluidState.massFrac(phaseIdx, CO2Idx));
        }
        else {
            Scalar result = 0;
            result +=
                Brine::gasEnthalpy(temperature, pressure) *
                fluidState.massFrac(gPhaseIdx, BrineIdx);
            result +=
                CO2::gasEnthalpy(temperature, pressure) *
                fluidState.massFrac(gPhaseIdx, CO2Idx);
            return result;
        }
    }

    /*!
     * \brief Given the phase composition, return the phase's internal
     *        energy [J/kg].
     */
    template <class FluidState>
    static Scalar phaseInternalEnergy(int phaseIdx,
                                 Scalar temperature,
                                 Scalar pressure,
                                 const FluidState &fluidState)
    {
        Scalar h = phaseEnthalpy(phaseIdx, temperature, pressure, fluidState);
        Scalar rho = phaseDensity(phaseIdx, temperature, pressure, fluidState);
        return h - pressure/rho;
    }

private:
    static Scalar gasDensity_(Scalar T,
                              Scalar pg,
                              Scalar xgH2O,
                              Scalar xgCO2)
    {
        Scalar gasDensity = CO2::gasDensity(T, pg);
        return gasDensity;
    }

    /***********************************************************************/
    /*                                                                     */
    /* Total brine density with dissolved CO2                              */
    /* rho_{b,CO2} = rho_w + contribution(salt) + contribution(CO2)        */
    /*                                                                     */
    /***********************************************************************/
    static Scalar liquidDensity_(Scalar T,
                                 Scalar pl,
                                 Scalar xlH2O,
                                 Scalar xlCO2)
    {
        if(T < 273.15) {
            OPM_THROW(NumericalProblem,
                       "Liquid density for Brine and CO2 is only "
                       "defined above 273.15K (is" << T << ")");
        }
        if(pl >= 2.5e8) {
            std::cout << "*** brine_co2_system:  outside range " << pl << " ****" << std::endl;
            pl = std::min(pl,2.5e8);
            /*
            OPM_THROW(NumericalProblem,
                       "Liquid density for Brine and CO2 is only "
                       "defined below 250MPa (is" << pl << ")");
            */
        }

        Scalar rho_brine = Brine::liquidDensity(T, pl);
        Scalar rho_pure = H2O::liquidDensity(T, pl);
        Scalar rho_lCO2 = liquidDensityWaterCO2_(T, pl, xlH2O, xlCO2);
        Scalar contribCO2 = rho_lCO2 - rho_pure;

        return rho_brine + contribCO2;
    }

    static Scalar liquidDensityWaterCO2_(Scalar temperature,
                                         Scalar pl,
                                         Scalar xlH2O,
                                         Scalar xlCO2)
    {
        const Scalar M_CO2 = CO2::molarMass();
        const Scalar M_H2O = H2O::molarMass();

        const Scalar tempC = temperature - 273.15;        /* tempC : temperature in °C */
        const Scalar rho_pure = H2O::liquidDensity(temperature, pl);
        xlH2O = 1.0 - xlCO2; // xlH2O is available, but in case of a pure gas phase
                             // the value of M_T for the virtual liquid phase can become very large
        const Scalar M_T = M_H2O * xlH2O + M_CO2 * xlCO2;
        const Scalar V_phi =
            (37.51 +
             tempC*(-9.585e-2 +
                    tempC*(8.74e-4 -
                           tempC*5.044e-7))) / 1.0e6;
        return 1/ (xlCO2 * V_phi/M_T + M_H2O * xlH2O / (rho_pure * M_T));
    }

    static Scalar liquidEnthalpyBrineCO2_(Scalar T,
                                          Scalar p,
                                          Scalar S,
                                          Scalar X_CO2_w)
    {
        /* X_CO2_w : mass fraction of CO2 in brine */

        /* same function as enthalpy_brine, only extended by CO2 content */

        /*Numerical coefficents from PALLISER*/
        static const Scalar f[] = {
            2.63500E-1, 7.48368E-6, 1.44611E-6, -3.80860E-10
        };

        /*Numerical coefficents from MICHAELIDES for the enthalpy of brine*/
        static const Scalar a[4][3] = {
            { 9633.6, -4080.0, +286.49 },
            { +166.58, +68.577, -4.6856 },
            { -0.90963, -0.36524, +0.249667E-1 },
            { +0.17965E-2, +0.71924E-3, -0.4900E-4 }
        };

        Scalar theta, h_NaCl;
        Scalar m, h_ls, h_ls1, d_h;
        Scalar S_lSAT, delta_h;
        int i, j;
        Scalar delta_hCO2, hg, hw;

        theta = T - 273.15;

        S_lSAT = f[0] + f[1]*theta + f[2]*theta*theta + f[3]*theta*theta*theta;
        /*Regularization*/
        if (S>S_lSAT) {
            S = S_lSAT;
        }

        hw = H2O::liquidEnthalpy(T, p) /1E3; /* kJ/kg */

        /*DAUBERT and DANNER*/
        /*U=*/h_NaCl = (3.6710E4*T + 0.5*(6.2770E1)*T*T - ((6.6670E-2)/3)*T*T*T
                        +((2.8000E-5)/4)*(T*T*T*T))/(58.44E3)- 2.045698e+02; /* kJ/kg */

        m = (1E3/58.44)*(S/(1-S));
        i = 0;
        j = 0;
        d_h = 0;

        for (i = 0; i<=3; i++) {
            for (j=0; j<=2; j++) {
                d_h = d_h + a[i][j] * pow(theta, i) * pow(m, j);
            }
        }
	/* heat of dissolution for halite according to Michaelides 1971 */
        delta_h = (4.184/(1E3 + (58.44 * m)))*d_h;

        /* Enthalpy of brine without CO2 */
        h_ls1 =(1-S)*hw + S*h_NaCl + S*delta_h; /* kJ/kg */

        /* heat of dissolution for CO2 according to Fig. 6 in Duan and Sun 2003. (kJ/kg)
	   In the relevant temperature ranges CO2 dissolution is
	   exothermal */
        delta_hCO2 = (-57.4375 + T * 0.1325) * 1000/44;

	/* enthalpy contribution of CO2 (kJ/kg) */
        hg = CO2::liquidEnthalpy(T, p)/1E3 + delta_hCO2;

        /* Enthalpy of brine with dissolved CO2 */
        h_ls = (h_ls1 - X_CO2_w*hw + hg*X_CO2_w)*1E3; /*J/kg*/

        return (h_ls);
    };
};

} // end namepace

#endif
