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

#include "config.h"

#include "BlackoilPVT.hpp"
#include <opm/core/utility/Units.hpp>
#include "MiscibilityDead.hpp"
#include "MiscibilityLiveOil.hpp"
#include "MiscibilityLiveGas.hpp"
#include "MiscibilityWater.hpp"
#include <opm/common/ErrorMacros.hpp>
#include <opm/core/utility/linearInterpolation.hpp>

#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/PvdoTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/PvdgTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/TableManager.hpp>
#include <opm/parser/eclipse/Parser/ParseContext.hpp>

using namespace Opm;

namespace Opm
{


    void BlackoilPVT::init(Opm::DeckConstPtr deck)
    {
        Opm::ParseContext parseContext;
        const auto eclipseState = std::make_shared<EclipseState>(deck , parseContext);
	region_number_ = 0;

	// Surface densities. Accounting for different orders in eclipse and our code.
	if (deck->hasKeyword("DENSITY")) {
        const auto& densityRecord =
            deck->getKeyword("DENSITY").getRecord(/*regionIdx=*/0);
	    densities_[Aqua]   = densityRecord.getItem("WATER").getSIDouble(0);
	    densities_[Vapour] = densityRecord.getItem("GAS").getSIDouble(0);
	    densities_[Liquid] = densityRecord.getItem("OIL").getSIDouble(0);
	} else {
	    OPM_THROW(std::runtime_error, "Input is missing DENSITY\n");
	}

        // Water PVT
        if (deck->hasKeyword("PVTW")) {
            water_props_.reset(new MiscibilityWater(deck->getKeyword("PVTW")));
        } else {
            water_props_.reset(new MiscibilityWater(0.5*Opm::prefix::centi*Opm::unit::Poise)); // Eclipse 100 default 
        }

        // Oil PVT
        const auto& tables     = eclipseState->getTableManager();
        const auto& pvdoTables = tables->getPvdoTables();
        const auto& pvtoTables = tables->getPvtoTables();
        if (!pvdoTables.empty()) {
            const auto& pvdoTable = pvdoTables.getTable<PvdoTable>(0);
            oil_props_.reset(new MiscibilityDead(pvdoTable));
        } else if (pvtoTables.empty()) {
            // PVTOTables is a std::vector<>
            const auto& pvtoTable = pvtoTables[0];
            oil_props_.reset(new MiscibilityLiveOil(pvtoTable));
        } else if (deck->hasKeyword("PVCDO")) {
            auto *misc_water = new MiscibilityWater(0);
            misc_water->initFromPvcdo(deck->getKeyword("PVCDO"));
            oil_props_.reset(misc_water);
        } else {
            OPM_THROW(std::runtime_error, "Input is missing PVDO and PVTO\n");
        }

	// Gas PVT
        const auto& pvdgTables = tables->getPvdgTables();
        const auto& pvtgTables = tables->getPvtgTables();
        if (!pvdgTables.empty()) {
            const auto& pvdgTable = pvdgTables.getTable<PvdgTable>(0);
            gas_props_.reset(new MiscibilityDead(pvdgTable));
        } else if (pvtgTables.empty()) {
            gas_props_.reset(new MiscibilityLiveGas(pvtgTables[0]));
        } else {
	    OPM_THROW(std::runtime_error, "Input is missing PVDG and PVTG\n");
        }
    }

    BlackoilPVT::CompVec BlackoilPVT::surfaceDensities() const
    {
        return densities_;
    }

    double BlackoilPVT::getViscosity(double press, const CompVec& surfvol, PhaseIndex phase) const
    {
        return propsForPhase(phase).getViscosity(region_number_, press, surfvol);
    }

    double BlackoilPVT::B(double press, const CompVec& surfvol, PhaseIndex phase) const
    {
        return propsForPhase(phase).B(region_number_, press, surfvol);
    }

    double BlackoilPVT::dBdp(double press, const CompVec& surfvol, PhaseIndex phase) const
    {
        return propsForPhase(phase).dBdp(region_number_, press, surfvol);
    }

    double BlackoilPVT::R(double press, const CompVec& surfvol, PhaseIndex phase) const
    {
        return propsForPhase(phase).R(region_number_, press, surfvol);
    }

    double BlackoilPVT::dRdp(double press, const CompVec& surfvol, PhaseIndex phase) const
    {
        return propsForPhase(phase).dRdp(region_number_, press, surfvol);
    }

    const MiscibilityProps& BlackoilPVT::propsForPhase(PhaseIndex phase) const
    {
        switch (phase) {
        case Aqua:
            return *water_props_;
        case Liquid:
            return *oil_props_;
        case Vapour:
            return *gas_props_;
        default:
            OPM_THROW(std::runtime_error, "Unknown phase accessed: " << phase);
        }
    }

    void BlackoilPVT::getViscosity(const std::vector<PhaseVec>& pressures,
                                   const std::vector<CompVec>& surfvol,
                                   std::vector<PhaseVec>& output) const
    {
        int num = pressures.size();
        output.resize(num);
        for (int phase = 0; phase < numPhases; ++phase) {
            propsForPhase(PhaseIndex(phase)).getViscosity(pressures, surfvol, phase, data1_);
            for (int i = 0; i < num; ++i) {
                output[i][phase] = data1_[i];
            }
        }
    }

    void BlackoilPVT::B(const std::vector<PhaseVec>& pressures,
                        const std::vector<CompVec>& surfvol,
                        std::vector<PhaseVec>& output) const
    {
        int num = pressures.size();
        output.resize(num);
        for (int phase = 0; phase < numPhases; ++phase) {
            propsForPhase(PhaseIndex(phase)).B(pressures, surfvol, phase, data1_);
            for (int i = 0; i < num; ++i) {
                output[i][phase] = data1_[i];
            }
        }
    }

    void BlackoilPVT::dBdp(const std::vector<PhaseVec>& pressures,
                           const std::vector<CompVec>& surfvol,
                           std::vector<PhaseVec>& output_B,
                           std::vector<PhaseVec>& output_dBdp) const
    {
        int num = pressures.size();
        output_B.resize(num);
        output_dBdp.resize(num);
        for (int phase = 0; phase < numPhases; ++phase) {
            propsForPhase(PhaseIndex(phase)).dBdp(pressures, surfvol, phase, data1_, data2_);
            for (int i = 0; i < num; ++i) {
                output_B[i][phase] = data1_[i];
                output_dBdp[i][phase] = data2_[i];
            }
        }
    }

    void BlackoilPVT::R(const std::vector<PhaseVec>& pressures,
                        const std::vector<CompVec>& surfvol,
                        std::vector<PhaseVec>& output) const
    {
        int num = pressures.size();
        output.resize(num);
        for (int phase = 0; phase < numPhases; ++phase) {
            propsForPhase(PhaseIndex(phase)).R(pressures, surfvol, phase, data1_);
            for (int i = 0; i < num; ++i) {
                output[i][phase] = data1_[i];
            }
        }
    }

    void BlackoilPVT::dRdp(const std::vector<PhaseVec>& pressures,
                           const std::vector<CompVec>& surfvol,
                           std::vector<PhaseVec>& output_R,
                           std::vector<PhaseVec>& output_dRdp) const
    {
        int num = pressures.size();
        output_R.resize(num);
        output_dRdp.resize(num);
        for (int phase = 0; phase < numPhases; ++phase) {
            propsForPhase(PhaseIndex(phase)).dRdp(pressures, surfvol, phase, data1_, data2_);
            for (int i = 0; i < num; ++i) {
                output_R[i][phase] = data1_[i];
                output_dRdp[i][phase] = data2_[i];
            }
        }
    }

} // namespace Opm
