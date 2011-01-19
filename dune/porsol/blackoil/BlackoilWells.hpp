/*
  Copyright 2011 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_BLACKOILWELLS_HEADER_INCLUDED
#define OPM_BLACKOILWELLS_HEADER_INCLUDED

#include <dune/porsol/blackoil/fluid/BlackoilDefs.hpp>
#include <dune/common/ErrorMacros.hpp>
#include <dune/common/SparseTable.hpp>
#include <dune/common/fvector.hh>
#include <vector>

// Forward declaration.
namespace Dune
{
    class EclipseGridParser;
}

namespace
{

int prod_control_mode(const std::string& control){
    const int num_prod_control_modes = 7;
    static std::string prod_control_modes[num_prod_control_modes] =
	{std::string("ORAT"), std::string("WRAT"), std::string("GRAT"),
	 std::string("LRAT"), std::string("RESV"), std::string("BHP"),
	 std::string("THP") };   //  @bsp missing GRUP
    int m = -1;
    for (int i=0; i<num_prod_control_modes; ++i) {
	if (control[0] == prod_control_modes[i][0]) {
	    m = i;
	    break;
	}
    }
    if (m >= 0) {
	return m;
    } else {
	THROW("Unknown well control mode = " << control << " in input file");
    }
}

int inje_control_mode(const std::string& control)
{
    const int num_inje_control_modes = 4;
    static std::string inje_control_modes[num_inje_control_modes] =
	{std::string("RATE"), std::string("BHP"), std::string("THP"),
	 std::string("GRUP") };   //  @bsp missing"RESV" ??? 
    int m = -1;
    for (int i=0; i<num_inje_control_modes; ++i) {
	if (control[0] == inje_control_modes[i][0]) {
	    m = i;
	    break;
	}
    }
 
    if (m >= 0) {
	return m;
    } else {
	THROW("Unknown well control mode = " << control << " in input file");
    }
}
} // anon namespace

namespace Opm
{



    /// A class designed to encapsulate a set of rate- or
    /// pressure-controlled wells in the black-oil setting.
    class BlackoilWells : public BlackoilDefs
    {
    public:
        void init(const Dune::EclipseGridParser& parser,
		  const Dune::CpGrid& grid);

        // Well-centric interface. Mostly used by pressure solver.
        int numWells() const;
        enum WellType { Injector, Producer };
        WellType type(int wellnum) const;
        enum WellControl { Rate, Pressure };
        WellControl control(int wellnum) const;
        double target(int wellnum) const;
        int numPerforations(int wellnum) const;
        int wellCell(int wellnum, int perfnum) const;
        double wellIndex(int wellnum, int perfnum) const;
        double pressureDelta(int wellnum, int perfnum) const;

        // Updating rates and pressures after pressure solve.
        void update(int num_cells,
                    const std::vector<double>& well_pressures,
                    const std::vector<double>& well_fluxes);

        // Cell-centric interface. Mostly used by transport solver.
        double perforationPressure(int cell) const;
        double wellToReservoirFlux(int cell) const;
        Dune::FieldVector<double, 3> injectionMixture(int cell) const;

    private:
        struct WellData { WellType type; WellControl control; double target; };
        std::vector<WellData> well_data_;
        struct PerfData { int cell; double well_index; double pdelta; };
        Dune::SparseTable<PerfData> perf_data_;
        std::vector<double> well_flux_;
        std::vector<double> well_pressure_;
        Dune::FieldVector<double, 3> injection_mixture_;
	std::vector<std::string> well_names_;
    };


    // ------------ Method implementations --------------

    // inline void BlackoilWells::init(const Dune::EclipseGridParser& parser)
    // {
    //     // Temporary test hack.
    //     WellData w1 = { Injector, Pressure, 2e7 };
    //     WellData w2 = { Producer, Pressure, 1e7 };
    //     well_data_.push_back(w1);
    //     well_data_.push_back(w2);
    //     PerfData p1 = { 0, 1e-11, 0.0 };
    //     PerfData p2 = { 99, 1e-11, 0.0 };
    //     perf_data_.appendRow(&p1, &p1 + 1);
    //     perf_data_.appendRow(&p2, &p2 + 1);
    //     injection_mixture_ = 0.0;
    //     injection_mixture_[Gas] = 1.0;
    // }

    inline void BlackoilWells::init(const Dune::EclipseGridParser& parser,
				    const Dune::CpGrid& grid)
    {
	std::vector<std::string> keywords;
	keywords.push_back("WELSPECS");
	keywords.push_back("COMPDAT");
	keywords.push_back("WELTARG");
	if (!parser.hasFields(keywords)) {
	    THROW("Needed field is missing in file");
	}
	if (!(parser.hasField("WCONINJE") || parser.hasField("WCONPROD")) ) {
	    THROW("Needed field is missing in file");
	}
        using namespace Dune;
	const WELSPECS& welspecs = parser.getWELSPECS();
	const int num_welspecs   = welspecs.welspecs.size();

	const COMPDAT& compdats = parser.getCOMPDAT();
	const int num_compdats  = compdats.compdat.size();

	const WELTARG& weltargs = parser.getWELTARG();
	const int num_weltargs  = weltargs.weltarg.size();	

	const WCONINJE& wconinjes = parser.getWCONINJE();
	const int num_wconinjes   = wconinjes.wconinje.size();

	const WCONPROD& wconprods = parser.getWCONPROD();
	const int num_wconprods   = wconprods.wconprod.size();

	// Get WELLSPECS data
	well_names_.reserve(num_welspecs);
	well_data_.reserve(num_welspecs);
	for (int i=0; i<num_welspecs; ++i) {
	    well_names_.push_back(welspecs.welspecs[i].name_);
	    WellData wd;
	    well_data_.push_back(wd);
	}

	// Get COMPDAT data   
	for (int kw=0; kw<num_compdats; ++kw) {
	    std::string name = compdats.compdat[kw].well_;
	    std::string::size_type len = name.find('*');
	    if (len != std::string::npos) {
		name = name.substr(0, len);
	    }

	    // global_cell is a map from compressed cells to Cartesian grid cells 
	    const std::vector<int>& global_cell = grid.globalCell();
	    const boost::array<int, 3>& cpgdim = grid.logicalCartesianSize();
	    std::map<int,int> cartesian_to_compressed;
	    for (int i=0; i<int(global_cell.size()); ++i) {
		cartesian_to_compressed.insert(std::make_pair(global_cell[i], i));
	    }
	    bool found = false;
	    for (int wix=0; wix<num_welspecs; ++wix) {
		if (well_names_[wix].compare(0,len, name) == 0) { //equal
		    int ix = compdats.compdat[kw].grid_ind_[0] - 1;
		    int jy = compdats.compdat[kw].grid_ind_[1] - 1;
		    int kz = compdats.compdat[kw].grid_ind_[2] - 1;
		    int cart_grid_indx = ix + cpgdim[0]*(jy + cpgdim[1]*kz);
		    std::map<int, int>::const_iterator cgit = 
			cartesian_to_compressed.find(cart_grid_indx);
		    if (cgit != cartesian_to_compressed.end()) {
			PerfData pd;
			pd.cell       = cgit->second;
			pd.well_index = 1.e20;
			pd.pdelta     = 0.0;
			perf_data_.appendRow(&pd, &pd + 1);
		    } else {
			THROW("Cell with i,j,k indices " << ix << ' ' << jy << ' '
			      << kz << " not found!");
		    }
		    found = true;
		    break;
		}
	    }
	    if (!found) {
		THROW("Undefined well name: " << compdats.compdat[kw].well_
		      << " in COMPDAT");
	    }
	}	

	// Get WCONINJE data
	for (int kw=0; kw<num_wconinjes; ++kw) {
	    std::string name = wconinjes.wconinje[kw].well_;
	    std::string::size_type len = name.find('*');
	    if (len != std::string::npos) {
		name = name.substr(0, len);
	    }

	    bool well_found = false;
	    for (int wix=0; wix<num_welspecs; ++wix) {
		if (well_names_[wix].compare(0,len, name) == 0) { //equal
		    well_found = true;
		    well_data_[wix].type = Injector;
		    int m = inje_control_mode(wconinjes.wconinje[kw].control_mode_);
		    switch(m) {
		    case 0:  // RATE
			well_data_[wix].control = Rate;
			well_data_[wix].target = 
			    wconinjes.wconinje[kw].surface_flow_max_rate_;
			break;
		    case 1:  // BHP
			well_data_[wix].control = Pressure;
			well_data_[wix].target = wconinjes.wconinje[kw].BHP_limit_;
			break;
		    case 2:  // THP
			well_data_[wix].control = Pressure;
			well_data_[wix].target = wconinjes.wconinje[kw].THP_limit_;
			break;
		    case 3:  // GRUP 
			well_data_[wix].control = Rate;   // @bsp ???
			break;
		    default:
			THROW("Unknown well control mode; WCONIJE  = "
			      << wconinjes.wconinje[kw].control_mode_
			      << " in input file");
		    }
		}
	    }
	    if (!well_found) {
		THROW("Undefined well name: " << wconinjes.wconinje[kw].well_
		      << " in WCONINJE");
	    }
	}	

	// Get WCONPROD data   
	for (int kw=0; kw<num_wconprods; ++kw) {
	    std::string name = wconprods.wconprod[kw].well_;
	    std::string::size_type len = name.find('*');
	    if (len != std::string::npos) {
		name = name.substr(0, len);
	    }

	    bool well_found = false;
	    for (int wix=0; wix<num_welspecs; ++wix) {
		if (well_names_[wix].compare(0,len, name) == 0) { //equal
		    well_found = true;
		    well_data_[wix].type = Producer;
		    int m = prod_control_mode(wconprods.wconprod[kw].control_mode_);
		    switch(m) {
		    case 0:  // ORAT
			well_data_[wix].control = Rate;
			well_data_[wix].target = wconprods.wconprod[kw].oil_max_rate_;
			break;
		    case 1:  // WRAT
			well_data_[wix].control = Rate;
			well_data_[wix].target = wconprods.wconprod[kw].water_max_rate_;
			break;
		    case 2:  // GRAT
			well_data_[wix].control = Rate;
			well_data_[wix].target = wconprods.wconprod[kw].gas_max_rate_;
			break;
		    case 3:  // LRAT
			well_data_[wix].control = Rate;
			well_data_[wix].target = wconprods.wconprod[kw].liquid_max_rate_;
			break;
		    case 4:  // RESV 
			well_data_[wix].control = Rate;
			well_data_[wix].target =
			    wconprods.wconprod[kw].fluid_volume_max_rate_;
			break;
		    case 5:  // BHP
			well_data_[wix].control = Pressure; 
			well_data_[wix].target = wconprods.wconprod[kw].BHP_limit_;
			break;
		    case 6:  // THP 
			well_data_[wix].control = Pressure;
			well_data_[wix].target = wconprods.wconprod[kw].THP_limit_;
			break;
		    default:
			THROW("Unknown well control mode; WCONPROD  = "
			      << wconprods.wconprod[kw].control_mode_
			      << " in input file");
		    }
		}
	    }
	    if (!well_found) {
		THROW("Undefined well name: " << wconprods.wconprod[kw].well_
		      << " in WCONPROD");
	    }
	}	

	// Get WELTARG data   
	for (int kw=0; kw<num_weltargs; ++kw) {
	    std::string name = weltargs.weltarg[kw].well_;
	    std::string::size_type len = name.find('*');
	    if (len != std::string::npos) {
		name = name.substr(0, len);
	    }

	    bool well_found = false;
	    for (int wix=0; wix<num_welspecs; ++wix) {
		if (well_names_[wix].compare(0,len, name) == 0) { //equal
		    well_found = true;
		    well_data_[wix].target = weltargs.weltarg[kw].new_value_;
		    break;
		}
	    }
	    if (!well_found) {
		THROW("Undefined well name: " << weltargs.weltarg[kw].well_
		      << " in WELTARG");
	    }
	}

	std::cout << "\t WELL DATA" << std::endl;
	for(int i=0; i< int(well_data_.size()); ++i) {
	    std::cout << i << ": " << well_data_[i].type << "  "
		      << well_data_[i].control << "  " << well_data_[i].target
		      << std::endl;
	}

	std::cout << "\n\t PERF DATA" << std::endl;
	for(int i=0; i< int(perf_data_.size()); ++i) {
	    for(int j=0; j< int(perf_data_[i].size()); ++j) {
		std::cout << i << ": " << perf_data_[i][j].cell << "  "
			  << perf_data_[i][j].well_index << "  "
			  << perf_data_[i][j].pdelta << std::endl;
	    }
	}
    }

    inline int BlackoilWells::numWells() const
    {
        return well_data_.size();
    }

    inline BlackoilWells::WellType BlackoilWells::type(int wellnum) const
    {
        return well_data_[wellnum].type;
    }

    inline BlackoilWells::WellControl BlackoilWells::control(int wellnum) const
    {
        return well_data_[wellnum].control;
    }

    inline double BlackoilWells::target(int wellnum) const
    {
        return well_data_[wellnum].target;
    }

    inline int BlackoilWells::numPerforations(int wellnum) const
    {
        return perf_data_[wellnum].size();
    }

    inline int BlackoilWells::wellCell(int wellnum, int perfnum) const
    {
        return perf_data_[wellnum][perfnum].cell;
    }

    inline double BlackoilWells::wellIndex(int wellnum, int perfnum) const
    {
        return perf_data_[wellnum][perfnum].well_index;
    }

    inline double BlackoilWells::pressureDelta(int wellnum, int perfnum) const
    {
        return perf_data_[wellnum][perfnum].pdelta;
    }

    inline void BlackoilWells::update(int num_cells,
                                      const std::vector<double>& well_pressures,
                                      const std::vector<double>& well_fluxes)
    {
        // Input is per perforation, data members store for all cells.
        ASSERT(perf_data_.dataSize() == int(well_pressures.size()));
        well_pressure_.resize(num_cells, -1e100);
        well_flux_.resize(num_cells, 0.0);
        int pcount = 0;
        for (int w = 0; w < numWells(); ++w) {
            for (int perf = 0; perf < numPerforations(w); ++perf) {
                int cell = wellCell(w, perf);
                well_pressure_[cell] = well_pressures[pcount];
                well_flux_[cell] = well_fluxes[pcount];
                ++pcount;
            }
        }
        ASSERT(pcount == perf_data_.dataSize());
    }

    inline double BlackoilWells::wellToReservoirFlux(int cell) const
    {
        return well_flux_[cell];
    }

    inline double BlackoilWells::perforationPressure(int cell) const
    {
        return well_pressure_[cell];
    }

    inline Dune::FieldVector<double, 3> BlackoilWells::injectionMixture(int cell) const
    {
        return injection_mixture_;
    }

} // namespace Opm


#endif // OPM_BLACKOILWELLS_HEADER_INCLUDED
