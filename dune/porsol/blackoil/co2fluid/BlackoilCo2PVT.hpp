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

#ifndef OPM_BLACKOILCO2PVT_HEADER_INCLUDED
#define OPM_BLACKOILCO2PVT_HEADER_INCLUDED

#define OPM_DEPRECATED __attribute__((deprecated))
#define OPM_DEPRECATED_MSG(msg) __attribute__((deprecated))

#include <opm/core/eclipse/EclipseGridParser.hpp>

#include <dune/porsol/blackoil/co2fluid/opm/common/exceptions.hh>
#include <dune/porsol/blackoil/co2fluid/opm/material/fluidsystems/brine_co2_system.hh>
#include <dune/porsol/blackoil/co2fluid/opm/material/fluidstates/compositionalfluidstate.hh>
#include <dune/porsol/blackoil/co2fluid/benchmark3co2tables.hh>

#include <dune/porsol/blackoil/fluid/MiscibilityProps.hpp>
#include <dune/porsol/blackoil/fluid/BlackoilDefs.hpp>
#include <string>


namespace Opm
{
    class BlackoilCo2PVT : public BlackoilDefs
    {
    public:        
        
        typedef Opm::Brine_CO2_System<Opm::Benchmark3::CO2Tables, false> FluidSystem;
        typedef Opm::CompositionalFluidState<double, FluidSystem> CompositionalFluidState;
	
	void init(const Opm::EclipseGridParser& ep);

        double getViscosity(double press,
                            const CompVec& surfvol,
			    PhaseIndex phase) const;
	CompVec surfaceDensities() const;
	double getSaturation(double press, 
	                     const CompVec& surfvol, 
	                     PhaseIndex phase) const;
        double B   (double press,
                    const CompVec& surfvol,
		    PhaseIndex phase) const;
        double dBdp(double press,
                    const CompVec& surfvol,
		    PhaseIndex phase) const;
        double R   (double press,
                    const CompVec& surfvol,
		    PhaseIndex phase) const;
        double dRdp(double press,
                    const CompVec& surfvol,
		    PhaseIndex phase) const;

        void getViscosity(const std::vector<PhaseVec>& pressures,
                          const std::vector<CompVec>& surfvol,
                          std::vector<PhaseVec>& output) const;
        void B(const std::vector<PhaseVec>& pressures,
               const std::vector<CompVec>& surfvol,
               std::vector<PhaseVec>& output) const;
        void dBdp(const std::vector<PhaseVec>& pressures,
                  const std::vector<CompVec>& surfvol,
                  std::vector<PhaseVec>& output_B,
                  std::vector<PhaseVec>& output_dBdp) const;
        void R(const std::vector<PhaseVec>& pressures,
               const std::vector<CompVec>& surfvol,
               std::vector<PhaseVec>& output) const;
        void dRdp(const std::vector<PhaseVec>& pressures,
                  const std::vector<CompVec>& surfvol,
                  std::vector<PhaseVec>& output_R,
                  std::vector<PhaseVec>& output_dRdp) const;

    private:
	CompVec surfaceDensities_;
	FluidSystem brineCo2_;
        struct SubState
        {               
            Opm::FieldVector<double,2> density;
            Opm::FieldMatrix<double,2,2> massfrac;
            Opm::FieldVector<double,2> phaseVolume;
            Opm::FieldVector<double,2> phaseViscosity;
            double saturation;
        };
        void computeState(SubState& ss, double zBrine, double zCO2, double pressure, double temperature = 300.0) const;
        enum {
            wPhase = FluidSystem::wPhaseIdx,
            nPhase = FluidSystem::nPhaseIdx,
            
            wComp = FluidSystem::BrineIdx,
            nComp = FluidSystem::CO2Idx    
        };
   
    }; // class BlackoilCo2PVT

    // ------------ Method implementations --------------

    void BlackoilCo2PVT::init(const Opm::EclipseGridParser& ep)
    {
	surfaceDensities_[Aqua]   = 1000.;
	surfaceDensities_[Vapour] = 2.0;
	surfaceDensities_[Liquid] = 1000.;

        brineCo2_.init();
    }

    BlackoilCo2PVT::CompVec BlackoilCo2PVT::surfaceDensities() const
    {
        return surfaceDensities_;
    }

    double BlackoilCo2PVT::getViscosity(double press, const CompVec& surfvol, PhaseIndex phase) const
    {
        SubState ss;
        computeState(ss, surfvol[Liquid], surfvol[Vapour], press);
        switch(phase) {
        case Aqua: return 1.0e-10;
        case Liquid: return ss.phaseViscosity[wPhase];
        case Vapour: return ss.phaseViscosity[nPhase];
        };
        return 1.0e-10;
    }


    double BlackoilCo2PVT::getSaturation(double press, const CompVec& surfvol, PhaseIndex phase) const
    {
        SubState ss;
        computeState(ss, surfvol[Liquid], surfvol[Vapour], press);
        switch(phase) {
        case Aqua: return 0.0;
        case Liquid: return ss.saturation;
        case Vapour: return 1.0 - ss.saturation;
        };
        return 0.0;
    }
    
    
    double BlackoilCo2PVT::B(double press, const CompVec& surfvol, PhaseIndex phase) const
    {
        SubState ss;
        computeState(ss, surfvol[Liquid], surfvol[Vapour], press);
        switch(phase) {
        case Aqua: return 1.0;
        case Liquid: return surfaceDensities_[Liquid]/(ss.massfrac[wPhase][wComp]*ss.density[wPhase]+1.0e-10);
        case Vapour: return surfaceDensities_[Vapour]/(ss.massfrac[nPhase][nComp]*ss.density[nPhase]+1.0e-10);
        };
        return 1.0;
    }

    double BlackoilCo2PVT::dBdp(double press, const CompVec& surfvol, PhaseIndex phase) const
    {
        const double dp = 100.;
        return (B(press+dp,surfvol,phase)-B(press,surfvol,phase))/dp;
    }

    double BlackoilCo2PVT::R(double press, const CompVec& surfvol, PhaseIndex phase) const
    {
        SubState ss;
        computeState(ss, surfvol[Oil], surfvol[Gas], press);
        switch(phase) {
        case Aqua: return 0.0;
        case Liquid: return (ss.massfrac[wPhase][nComp]*surfaceDensities_[Liquid])/(ss.massfrac[wPhase][wComp]*surfaceDensities_[Vapour]+1.0e-10);
        case Vapour: return (ss.massfrac[nPhase][wComp]*surfaceDensities_[Vapour])/(ss.massfrac[nPhase][nComp]*surfaceDensities_[Liquid]+1.0e-10);
        };
        return 0.0;
    }

    double BlackoilCo2PVT::dRdp(double press, const CompVec& surfvol, PhaseIndex phase) const
    {
        const double dp = 100.;
        return (R(press+dp,surfvol,phase)-R(press,surfvol,phase))/dp;
    }

    void BlackoilCo2PVT::getViscosity(const std::vector<PhaseVec>& pressures,
                                   const std::vector<CompVec>& surfvol,
                                   std::vector<PhaseVec>& output) const
    {
        int num = pressures.size();
        output.resize(num);
        for (int i = 0; i < num; ++i) {
            output[i][Aqua] = getViscosity(pressures[i][Aqua],surfvol[i],Aqua);
            output[i][Liquid] = getViscosity(pressures[i][Liquid],surfvol[i],Liquid);
            output[i][Vapour] = getViscosity(pressures[i][Vapour],surfvol[i],Vapour);
        }
    }

    void BlackoilCo2PVT::B(const std::vector<PhaseVec>& pressures,
                        const std::vector<CompVec>& surfvol,
                        std::vector<PhaseVec>& output) const
    {
        int num = pressures.size();
        output.resize(num);
        for (int i = 0; i < num; ++i) {
            output[i][Aqua] = B(pressures[i][Aqua],surfvol[i],Aqua);
            output[i][Liquid] = B(pressures[i][Liquid],surfvol[i],Liquid);
            output[i][Vapour] = B(pressures[i][Vapour],surfvol[i],Vapour);
        }
    }

    void BlackoilCo2PVT::dBdp(const std::vector<PhaseVec>& pressures,
                           const std::vector<CompVec>& surfvol,
                           std::vector<PhaseVec>& output_B,
                           std::vector<PhaseVec>& output_dBdp) const
    {
        int num = pressures.size();
        output_B.resize(num);
        output_dBdp.resize(num);
        for (int i = 0; i < num; ++i) {
            output_B[i][Aqua] = B(pressures[i][Aqua],surfvol[i],Aqua);
            output_B[i][Liquid] = B(pressures[i][Liquid],surfvol[i],Liquid);
            output_B[i][Vapour] = B(pressures[i][Vapour],surfvol[i],Vapour);
            output_dBdp[i][Aqua] = dBdp(pressures[i][Aqua],surfvol[i],Aqua);
            output_dBdp[i][Liquid] = dBdp(pressures[i][Liquid],surfvol[i],Liquid);
            output_dBdp[i][Vapour] = dBdp(pressures[i][Vapour],surfvol[i],Vapour);
        }
    }

    void BlackoilCo2PVT::R(const std::vector<PhaseVec>& pressures,
                        const std::vector<CompVec>& surfvol,
                        std::vector<PhaseVec>& output) const
    {
        int num = pressures.size();
        output.resize(num);
        for (int i = 0; i < num; ++i) {
            output[i][Aqua] = R(pressures[i][Aqua],surfvol[i],Aqua);
            output[i][Liquid] = R(pressures[i][Liquid],surfvol[i],Liquid);
            output[i][Vapour] = R(pressures[i][Vapour],surfvol[i],Vapour);
        }
    }

    void BlackoilCo2PVT::dRdp(const std::vector<PhaseVec>& pressures,
                           const std::vector<CompVec>& surfvol,
                           std::vector<PhaseVec>& output_R,
                           std::vector<PhaseVec>& output_dRdp) const
    {
        int num = pressures.size();
        output_R.resize(num);
        output_dRdp.resize(num);
        for (int i = 0; i < num; ++i) {
            output_R[i][Aqua] = R(pressures[i][Aqua],surfvol[i],Aqua);
            output_R[i][Liquid] = R(pressures[i][Liquid],surfvol[i],Liquid);
            output_R[i][Vapour] = R(pressures[i][Vapour],surfvol[i],Vapour);
            output_dRdp[i][Aqua] = dRdp(pressures[i][Aqua],surfvol[i],Aqua);
            output_dRdp[i][Liquid] = dRdp(pressures[i][Liquid],surfvol[i],Liquid);
            output_dRdp[i][Vapour] = dRdp(pressures[i][Vapour],surfvol[i],Vapour);
        }
    }
    
    void BlackoilCo2PVT::computeState(BlackoilCo2PVT::SubState& ss, double zBrine, double zCO2, double pressure, double temperature) const
    {
               
        CompositionalFluidState state;     
        state.setTemperature(temperature);
        state.setPressure(wPhase, pressure);
        state.setPressure(nPhase, pressure);
        
        double massH20 = surfaceDensities_[Liquid]*zBrine;
        double massCO2 = surfaceDensities_[Vapour]*zCO2;
        
        // A priori, assume presence of both phases
        brineCo2_.computeEquilibrium(state); 
        ss.density[wPhase] = state.density(wPhase);
        ss.density[nPhase] = state.density(nPhase);
        ss.massfrac[wPhase][nComp] = state.massFraction(wPhase, nComp);
        ss.massfrac[nPhase][wComp] = state.massFraction(nPhase, wComp);
        ss.massfrac[wPhase][wComp] = 1.0 - ss.massfrac[wPhase][nComp];
        ss.massfrac[nPhase][nComp] = 1.0 - ss.massfrac[nPhase][wComp];

        double detX = ss.massfrac[wPhase][wComp]*ss.massfrac[nPhase][nComp]-ss.massfrac[wPhase][nComp]*ss.massfrac[nPhase][wComp];
        ss.phaseVolume[wPhase] = (massH20*ss.massfrac[nPhase][nComp] - massCO2*ss.massfrac[nPhase][wComp])/(ss.density[wPhase]*detX);
        ss.phaseVolume[nPhase] = (massCO2*ss.massfrac[wPhase][wComp] - massH20*ss.massfrac[wPhase][nComp])/(ss.density[nPhase]*detX);
        
        // Determine number of phase
        if (ss.phaseVolume[wPhase] > 0.0 && ss.phaseVolume[nPhase] > 0.0) { // Both phases
            ss.saturation = ss.phaseVolume[wPhase]/(ss.phaseVolume[wPhase]+ss.phaseVolume[nPhase]);
            state.setSaturation(wPhase, ss.saturation);
            state.setSaturation(nPhase, 1.0 - ss.saturation);
        }
        else if (ss.phaseVolume[wPhase] <= 0.0) { // Wetting phase only
            ss.saturation = 0.0;
            // Gas phase:
            ss.massfrac[nPhase][nComp] = massCO2/(massCO2+massH20);
            ss.massfrac[nPhase][wComp] = 1.0 - ss.massfrac[nPhase][nComp];
            double M1 = FluidSystem::molarMass(wComp);
            double M2 = FluidSystem::molarMass(nComp);
            double avgMolarMass = M1*M2/(M2 + ss.massfrac[nPhase][nComp]*(M1 - M2));
            state.setMoleFraction(nPhase, nComp, ss.massfrac[nPhase][nComp]*avgMolarMass/M2);
            state.setMoleFraction(nPhase, wComp, ss.massfrac[nPhase][wComp]*avgMolarMass/M1);
            ss.density[nPhase] = brineCo2_.phaseDensity(nPhase, state.temperature(nPhase), state.pressure(nPhase), state);
            state.setDensity(nPhase, ss.density[nPhase]);
            ss.phaseVolume[nPhase] = (massH20+massCO2)/ss.density[nPhase];
            state.setSaturation(nPhase, 1.0 - ss.saturation);
            // Virtual properties of non-existing liquid phase:
            brineCo2_.computeEquilibrium(state, nPhase);
            ss.massfrac[wPhase][wComp] = state.massFraction(wPhase, wComp);
            ss.massfrac[wPhase][nComp] = state.massFraction(wPhase, nComp);
            ss.density[wPhase] = state.density(wPhase);
            ss.phaseVolume[wPhase] = 0.0;
            state.setSaturation(wPhase, ss.saturation);
        }
        else if (ss.phaseVolume[nPhase] <= 0.0) { // Non-wetting phase only
            ss.saturation = 1.0;
            // Liquid phase:
            ss.massfrac[wPhase][wComp] = massH20/(massCO2+massH20);
            ss.massfrac[wPhase][nComp] = 1.0 - ss.massfrac[wPhase][wComp];
            double M1 = FluidSystem::molarMass(wComp);
            double M2 = FluidSystem::molarMass(nComp);
            double avgMolarMass = M1*M2/(M2 + ss.massfrac[wPhase][nComp]*(M1 - M2));
            state.setMoleFraction(wPhase, nComp, ss.massfrac[wPhase][nComp]*avgMolarMass/M2);
            state.setMoleFraction(wPhase, wComp, ss.massfrac[wPhase][wComp]*avgMolarMass/M1);
            ss.density[wPhase] = brineCo2_.phaseDensity(wPhase, state.temperature(wPhase), state.pressure(wPhase), state);
            state.setDensity(wPhase, ss.density[wPhase]);
            ss.phaseVolume[wPhase] = (massH20+massCO2)/ss.density[wPhase];
            state.setSaturation(wPhase, ss.saturation);
            // Virtual properties of non-existing gas phase:
            brineCo2_.computeEquilibrium(state, wPhase);
            ss.massfrac[nPhase][nComp] = state.massFraction(nPhase, nComp);
            ss.massfrac[nPhase][wComp] = state.massFraction(nPhase, wComp);
            ss.density[nPhase] = state.density(nPhase);
            ss.phaseVolume[nPhase] = 0.0;
            state.setSaturation(nPhase, 1.0 - ss.saturation);
        } 
        
        ss.phaseViscosity[wPhase] = brineCo2_.phaseViscosity(wPhase, temperature, pressure, state);
        ss.phaseViscosity[nPhase] = brineCo2_.phaseViscosity(nPhase, temperature, pressure, state);        

    }
    
    
    

} // Opm


#endif // OPM_BLACKOILCO2PVT_HEADER_INCLUDED
