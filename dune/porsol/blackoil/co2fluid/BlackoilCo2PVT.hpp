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

        void generateBlackOilTables(double temperature);

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
        double temperature_;
        struct SubState
        {               
            Opm::FieldVector<double,2> density;
            Opm::FieldMatrix<double,2,2> massfrac;
            Opm::FieldVector<double,2> phaseVolume;
            Opm::FieldVector<double,2> phaseViscosity;
            double saturation;
        };
        void computeState(SubState& ss, double zBrine, double zCO2, double pressure) const;
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
	surfaceDensities_[Water]   = 1000.;
	surfaceDensities_[Gas] = 2.0;
	surfaceDensities_[Oil] = 1000.;

        temperature_ = 300.;

        brineCo2_.init();
    }

    void BlackoilCo2PVT::generateBlackOilTables(double temperature)
    {
        std::cout << "\n Generating pvt tables for the eclipse black oil formulation\n using the oil component as brine and the gas component as co_2." << std::endl;
        if (std::fabs(temperature-400.) > 100.0) {
            std::cout << "T=" << temperature << " is outside the allowed range [300,500] Kelvin" << std::endl;
            exit(-1);
        }

        temperature_ = temperature;

        CompVec z;
        z[Water] = 0.0;
        z[Oil] = 1.0;
        z[Gas] = 1.0e6;

        std::ofstream black("blackoil_pvt");
        black.precision(8);
        black << std::fixed << std::showpoint;

        const double pMin=150.e5;
        const double pMax=400.e5;
        const unsigned int np=11;
        std::vector<double> pValue(np+1,0.0);
        std::vector<double> rs(np+1,0.0);

        pValue[0] = 101330.0;
        rs[0] = 0.0;

        // Buble points 
        z[Gas] = 1.0e4;
        for (unsigned int i=0; i<np; ++i) {
            pValue[i+1] = pMin + i*(pMax-pMin)/(np-1);
            rs[i+1] = R(pValue[i+1], z, Liquid);       
        }

        const unsigned int np_fill_in = 10;
        const double dr = (rs[1] - rs[0])/(np_fill_in+1);
        const double dp = (pValue[1] - pValue[0])/(np_fill_in+1);
        rs.insert(rs.begin()+1,np_fill_in,0.0);
        pValue.insert(pValue.begin()+1,np_fill_in,0.0);
        for (unsigned int i=1; i<=np_fill_in; ++i) {
           rs[i] = rs[i-1] + dr;
           pValue[i] = pValue[i-1] + dp;
        }

        // Brine with dissolved co_2 ("live oil")
        black << "PVTO\n";
        black << "--     Rs       Pbub        Bo          Vo\n";
        black << "--  sm3/sm3     barsa       rm3/sm3     cP\n";

        // Undersaturated
        for (unsigned int i=0; i<np+np_fill_in; ++i) {
            z[Gas] = rs[i];
            black << std::setw(14) << rs[i];
            for (unsigned int j=i; j<np+1+np_fill_in; ++j) {
                if (j<=np_fill_in) {
                   if (j==i) black << std::setw(14) << pValue[j]*1.e-5 << std::setw(14) << 1.0-j*0.001 << std::setw(14) << 1.06499;
                   continue; 
                }
                if (j>i) black << std::endl << std::setw(14) << ' ';
                black << std::setw(14) << pValue[j]*1.e-5
                      << std::setw(14) << B(pValue[j], z, Liquid)
                      << std::setw(14) << getViscosity(pValue[j], z, Liquid)*1.e3;     
            }
            black << " /" <<  std::endl;
        }
        black << "/ " <<  std::endl;

        // We provide tables for co_2 both with and without dissolved water:

        // Co_2 neglecting dissolved water ("dry gas")
        black << "\nPVDG\n";
        black << "--    Pg          Bg            Vg\n";
        black << "--   barsa        rm3/sm3       cP\n";
 
        for (unsigned int i=0; i<np; ++i) {
            z[Oil] = 0.0;
            z[Gas] = 1.0;
            black << std::setw(14) << pValue[i+np_fill_in+1]*1.e-5
                  << std::setw(14) << B(pValue[i+np_fill_in+1], z, Vapour)
                  << std::setw(14) << getViscosity(pValue[i+np_fill_in+1], z, Vapour)*1.e3
                  << std::endl;     
        }
        black << "/ " <<  std::endl;

        // Co_2 with dissolved water ("wet gas")
        black << "\nPVTG\n";
        black << "--    Pg          Rv            Bg            Vg\n";
        black << "--   barsa        sm3/sm3       rm3/sm3       cP\n";
        for (unsigned int i=0; i<np; ++i) {
            z[Oil] = 1000.0;
            z[Gas] = 1.0;
            black << std::setw(14) << pValue[i+np_fill_in+1]*1.e-5
                  << std::setw(14) << R(pValue[i+np_fill_in+1], z, Vapour)
                  << std::setw(14) << B(pValue[i+np_fill_in+1], z, Vapour)
                  << std::setw(14) << getViscosity(pValue[i+np_fill_in+1], z, Vapour)*1.e3
                  << std::endl;
            z[Oil] = 0.0;
            black << std::setw(14) << ' '
                  << std::setw(14) << R(pValue[i+np_fill_in+1], z, Vapour)
                  << std::setw(14) << B(pValue[i+np_fill_in+1], z, Vapour)
                  << std::setw(14) << getViscosity(pValue[i+np_fill_in+1], z, Vapour)*1.e3
                  << " /" << std::endl;          
        }
        black << "/ " <<  std::endl;
        black << std::endl;
        std::cout << " Pvt tables for temperature=" << temperature << " Kelvin is written to file blackoil_pvt. " << std::endl;
        std::cout << " NOTE that the file contains tables for both PVDG (dry gas) and PVTG (wet gas)." << std::endl;
    }

    BlackoilCo2PVT::CompVec BlackoilCo2PVT::surfaceDensities() const
    {
        return surfaceDensities_;
    }

    double BlackoilCo2PVT::getViscosity(double press, const CompVec& surfvol, PhaseIndex phase) const
    {
        SubState ss;
        computeState(ss, surfvol[Oil], surfvol[Gas], press);
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
        computeState(ss, surfvol[Oil], surfvol[Gas], press);
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
        computeState(ss, surfvol[Oil], surfvol[Gas], press);
        switch(phase) {
        case Aqua: return 1.0;
        case Liquid: return surfaceDensities_[Oil]/(ss.massfrac[wPhase][wComp]*ss.density[wPhase]+1.0e-10);
        case Vapour: return surfaceDensities_[Gas]/(ss.massfrac[nPhase][nComp]*ss.density[nPhase]+1.0e-10);
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
        case Liquid: return (ss.massfrac[wPhase][nComp]*surfaceDensities_[Oil])/(ss.massfrac[wPhase][wComp]*surfaceDensities_[Gas]+1.0e-10);
        case Vapour: return (ss.massfrac[nPhase][wComp]*surfaceDensities_[Gas])/(ss.massfrac[nPhase][nComp]*surfaceDensities_[Oil]+1.0e-10);
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
        SubState ss;
        for (int i = 0; i < num; ++i) {
            computeState(ss, surfvol[i][Oil], surfvol[i][Gas], pressures[i][Liquid]);
            output[i][Aqua] = 1.0e-10;
            output[i][Liquid] = ss.phaseViscosity[wPhase];
            output[i][Vapour] = ss.phaseViscosity[nPhase];
        }
    }

    void BlackoilCo2PVT::B(const std::vector<PhaseVec>& pressures,
                        const std::vector<CompVec>& surfvol,
                        std::vector<PhaseVec>& output) const
    {
        int num = pressures.size();
        output.resize(num);
        SubState ss;
        for (int i = 0; i < num; ++i) {
            computeState(ss, surfvol[i][Oil], surfvol[i][Gas], pressures[i][Liquid]);
            output[i][Aqua] = 1.0;
            output[i][Liquid] = surfaceDensities_[Oil]/(ss.massfrac[wPhase][wComp]*ss.density[wPhase]+1.0e-10);
            output[i][Vapour] = surfaceDensities_[Gas]/(ss.massfrac[nPhase][nComp]*ss.density[nPhase]+1.0e-10);
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
        SubState ss;
        const double dp = 100.;
        for (int i = 0; i < num; ++i) {
            computeState(ss, surfvol[i][Oil], surfvol[i][Gas], pressures[i][Liquid]);
            output_B[i][Aqua] = 1.0;
            output_B[i][Liquid] = surfaceDensities_[Oil]/(ss.massfrac[wPhase][wComp]*ss.density[wPhase]+1.0e-10);
            output_B[i][Vapour] = surfaceDensities_[Gas]/(ss.massfrac[nPhase][nComp]*ss.density[nPhase]+1.0e-10);
            computeState(ss, surfvol[i][Oil], surfvol[i][Gas], pressures[i][Liquid]+dp);
            output_dBdp[i][Aqua] = 0.0;
            output_dBdp[i][Liquid] = (surfaceDensities_[Oil]/(ss.massfrac[wPhase][wComp]*ss.density[wPhase]+1.0e-10) - output_B[i][Liquid])/dp;
            output_dBdp[i][Vapour] = (surfaceDensities_[Gas]/(ss.massfrac[nPhase][nComp]*ss.density[nPhase]+1.0e-10) - output_B[i][Vapour])/dp;
        }
    }

    void BlackoilCo2PVT::R(const std::vector<PhaseVec>& pressures,
                        const std::vector<CompVec>& surfvol,
                        std::vector<PhaseVec>& output) const
    {
        int num = pressures.size();
        output.resize(num);
        SubState ss;
        for (int i = 0; i < num; ++i) {
            computeState(ss, surfvol[i][Oil], surfvol[i][Gas], pressures[i][Liquid]);
            output[i][Aqua] = 0.0;
            output[i][Liquid] = (ss.massfrac[wPhase][nComp]*surfaceDensities_[Oil])/(ss.massfrac[wPhase][wComp]*surfaceDensities_[Gas]+1.0e-10);
            output[i][Vapour] = (ss.massfrac[nPhase][wComp]*surfaceDensities_[Gas])/(ss.massfrac[nPhase][nComp]*surfaceDensities_[Oil]+1.0e-10);
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
        SubState ss;
        const double dp = 100.;
        for (int i = 0; i < num; ++i) {
            computeState(ss, surfvol[i][Oil], surfvol[i][Gas], pressures[i][Liquid]);
            output_R[i][Aqua] = 0.0;
            output_R[i][Liquid] = (ss.massfrac[wPhase][nComp]*surfaceDensities_[Oil])/(ss.massfrac[wPhase][wComp]*surfaceDensities_[Gas]+1.0e-10);
            output_R[i][Vapour] = (ss.massfrac[nPhase][wComp]*surfaceDensities_[Gas])/(ss.massfrac[nPhase][nComp]*surfaceDensities_[Oil]+1.0e-10);
            computeState(ss, surfvol[i][Oil], surfvol[i][Gas], pressures[i][Liquid]+dp);
            output_dRdp[i][Aqua] = 0.0;
            output_dRdp[i][Liquid] = ((ss.massfrac[wPhase][nComp]*surfaceDensities_[Oil])/(ss.massfrac[wPhase][wComp]*surfaceDensities_[Gas]+1.0e-10) - output_R[i][Liquid])/dp;
            output_dRdp[i][Vapour] = ((ss.massfrac[nPhase][wComp]*surfaceDensities_[Gas])/(ss.massfrac[nPhase][nComp]*surfaceDensities_[Oil]+1.0e-10) - output_R[i][Vapour])/dp;
        }
    }
    
    void BlackoilCo2PVT::computeState(BlackoilCo2PVT::SubState& ss, double zBrine, double zCO2, double pressure) const
    {
               
        CompositionalFluidState state;     
        state.setTemperature(temperature_);
        state.setPressure(wPhase, pressure);
        state.setPressure(nPhase, pressure);
        
        double massH20 = surfaceDensities_[Oil]*zBrine;
        double massCO2 = surfaceDensities_[Gas]*zCO2;
        
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
        
        ss.phaseViscosity[wPhase] = brineCo2_.phaseViscosity(wPhase, temperature_, pressure, state);
        ss.phaseViscosity[nPhase] = brineCo2_.phaseViscosity(nPhase, temperature_, pressure, state);        

    }
    
    
    

} // Opm


#endif // OPM_BLACKOILCO2PVT_HEADER_INCLUDED
