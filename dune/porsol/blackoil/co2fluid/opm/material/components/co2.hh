/*****************************************************************************
 *   Copyright (C) 2009-2010 by Melanie Darcis                               *
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
 * \brief A class for the CO2 fluid properties
 */
#ifndef OPM_CO2_HH
#define OPM_CO2_HH

#include <dune/porsol/blackoil/co2fluid/opm/common/exceptions.hh>
#include <dune/porsol/blackoil/co2fluid/opm/material/components/component.hh>

#include <cmath>
#include <iostream>

namespace Opm
{
/*!
 * \brief A class for the CO2 fluid properties
 */
template <class Scalar, class CO2Tables>
class CO2 : public Component<Scalar, CO2<Scalar, CO2Tables> >
{
public:
    /*!
     * \brief A human readable name for the CO2.
     */
    static const char *name()
    { return "CO2"; }

    /*!
     * \brief The mass in [kg] of one mole of CO2.
     */
    static Scalar molarMass()
    { return 44e-3; }

    /*!
     * \brief Returns the critical temperature [K] of CO2
     */
    static Scalar criticalTemperature()
    { return 273.15 + 30.95; /* [K] */ }

    /*!
     * \brief Returns the critical pressure [Pa] of CO2
     */
    static Scalar criticalPressure()
    { return 73.8e5; /* [N/m^2] */ }

    /*!
     * \brief Returns the temperature [K]at CO2's triple point.
     */
    static Scalar tripleTemperature()
    { return 273.15 - 56.35; /* [K] */ }

    /*!
     * \brief Returns the pressure [Pa] at CO2's triple point.
     */
    static Scalar triplePressure()
    { return 5.11e5; /* [N/m^2] */ }

    /*!
     * \brief The vapor pressure in [N/m^2] of pure CO2
     *        at a given temperature.
     *
     * See:
     *
     * TODO
     */
    static Scalar vaporPressure(Scalar T)
    { OPM_THROW(Opm::NotImplemented, "vaporPressure of CO2"); }

    /*!
     * \brief Specific enthalpy of gaseous CO2 [J/kg].
     */
    static Scalar gasEnthalpy(Scalar temperature,
                              Scalar pressure)
    {
        // for critical CO2 there's no difference between liquid and
        // gas
        return
            CO2Tables::tabulatedEnthalpy.at(temperature, pressure);
    }

    /*!
     * \brief Specific enthalpy of liquid CO2 [J/kg].
     */
    static Scalar liquidEnthalpy(Scalar temperature,
                                 Scalar pressure)
    {
        // for critical CO2 there's no difference between liquid and
        // gas
        return gasEnthalpy(temperature, pressure);
    }

    /*!
     * \brief Specific internal energy of CO2 [J/kg].
     */
    static Scalar gasInternalEnergy(Scalar temperature,
                                    Scalar pressure)
    {
        // for critical CO2 there's no difference between liquid and
        // gas
        Scalar h = gasEnthalpy(temperature, pressure);

        if (pressure < CO2Tables::tabulatedDensity.yMin()) {
            Scalar p0 = CO2Tables::tabulatedDensity.yMin();
            Scalar rho0 = CO2Tables::tabulatedDensity(temperature, p0);
            return h - p0/rho0;
        }

        Scalar rho = gasDensity(temperature, pressure);
        return h - (pressure / rho);
    }

    /*!
     * \brief Specific internal energy of liquid CO2 [J/kg].
     */
    static Scalar liquidInternalEnergy(Scalar temperature,
                                       Scalar pressure)
    {
        // for critical CO2 there's no difference between liquid and
        // gas
        return gasInternalEnergy(temperature, pressure);
    }

    /*!
     * \brief The density of CO2 at a given pressure and temperature [kg/m^3].
    */
    static Scalar gasDensity(Scalar temperature, Scalar pressure)
    {
        if (pressure < CO2Tables::tabulatedDensity.yMin()) {
            Scalar p0 = CO2Tables::tabulatedDensity.yMin();
            Scalar rho0 = CO2Tables::tabulatedDensity.at(temperature, p0);
            return pressure/p0 * rho0;
        }
        return CO2Tables::tabulatedDensity.at(temperature, pressure);
    }

    /*!
     * \brief The density of pure CO2 at a given pressure and temperature [kg/m^3].
     */
    static Scalar liquidDensity(Scalar temperature, Scalar pressure)
    {
        // we assume critical CO2!
        return gasDensity(temperature, pressure);
    }

    /*!
     * \brief The dynamic viscosity [N/m^3*s] of CO2.
     * Equations given in: - Vesovic et al., 1990
     *                        - Fenhour etl al., 1998
     */
    static Scalar gasViscosity(Scalar temperature, Scalar pressure)
    {
        static const double a0 = 0.235156;
        static const double a1 = -0.491266;
        static const double a2 = 5.211155E-2;
        static const double a3 = 5.347906E-2;
        static const double a4 = -1.537102E-2;

        static const double d11 = 0.4071119E-2;
        static const double d21 = 0.7198037E-4;
        static const double d64 = 0.2411697E-16;
        static const double d81 = 0.2971072E-22;
        static const double d82 = -0.1627888E-22;

        static const double ESP = 251.196;

        double mu0, SigmaStar, TStar;
        double dmu, rho;
        double visco_CO2;

        if(temperature < 275.) // regularisation
        {
            // temperature = 275;
            OPM_THROW(Opm::NumericalProblem,
                       "ConstrelCO2: Temperature " << temperature << " out of range at " << __FILE__ << ":" << __LINE__);
        }


        TStar = temperature/ESP;

        /* mu0: viscosity in zero-density limit */
        SigmaStar = exp(a0 + a1*log(TStar)
                        + a2*log(TStar)*log(TStar)
                        + a3*log(TStar)*log(TStar)*log(TStar)
                        + a4*log(TStar)*log(TStar)*log(TStar)*log(TStar) );

        mu0 = 1.00697*sqrt(temperature) / SigmaStar;

        /* dmu : excess viscosity at elevated density */
        rho = gasDensity(temperature, pressure); /* CO2 mass density [kg/m^3] */

        dmu = d11*rho + d21*rho*rho + d64*pow(rho,6)/(TStar*TStar*TStar)
            + d81*pow(rho,8) + d82*pow(rho,8)/TStar;

        /* dmucrit : viscosity increase near the critical point */

        // False (Lybke 2July2007)
        //e1 = 5.5930E-3;
        //e2 = 6.1757E-5;
        //e4 = 2.6430E-11;
        //dmucrit = e1*rho + e2*rho*rho + e4*rho*rho*rho;
        //visco_CO2 = (mu0 + dmu + dmucrit)/1.0E6;   /* conversion to [Pa s] */

        visco_CO2 = (mu0 + dmu)/1.0E6;   /* conversion to [Pa s] */

        return visco_CO2;
    };

    /*!
     * \brief The dynamic viscosity [N/m^3*s] of pure CO2.
     */
    static Scalar liquidViscosity(Scalar temperature, Scalar pressure)
    {
      // no difference for supercritical CO2
      return gasViscosity(temperature, pressure);;
    };
};

} // end namepace

#endif
