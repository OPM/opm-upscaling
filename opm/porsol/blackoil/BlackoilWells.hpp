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

#include <opm/porsol/blackoil/fluid/BlackoilDefs.hpp>
#include <opm/core/io/eclipse/EclipseGridParser.hpp>
#include <dune/grid/CpGrid.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/SparseTable.hpp>
#include <opm/porsol/common/Rock.hpp>
#include <dune/common/fvector.hh>
#include <vector>


namespace Opm
{
    /// A class designed to encapsulate a set of rate- or
    /// pressure-controlled wells in the black-oil setting.
    class BlackoilWells : public BlackoilDefs
    {
    public:
        void init(const Opm::EclipseGridParser& parser,
		  const Dune::CpGrid& grid,
		  const Opm::Rock<3>& rock);

        // Well-centric interface. Mostly used by pressure solver.
        int numWells() const;
        enum WellType { Injector, Producer };
        WellType type(int wellnum) const;
        enum WellControl { Rate, Pressure };
        WellControl control(int wellnum) const;
        double target(int wellnum) const;
        double referenceDepth(int wellnum) const;
        int numPerforations(int wellnum) const;
        int wellCell(int wellnum, int perfnum) const;
        double wellIndex(int wellnum, int perfnum) const;

        // Updating rates and pressures after pressure solve.
        void update(int num_cells,
                    const std::vector<double>& well_perf_pressures,
                    const std::vector<double>& well_perf_fluxes);

        // Cell-centric interface. Mostly used by transport solver.
        double perforationPressure(int cell) const;
        double wellToReservoirFlux(int cell) const;
        CompVec injectionMixture(int cell) const;
   
        // Simple well report ...
        class WellReport
        {
        public:
            static WellReport* report()
            {
                if (well_report_ == 0x0)
                  well_report_ = new WellReport;
                return well_report_;
            }
            void clearAll()
            {
                perfPressure.clear();
                cellPressure.clear();
                massRate.clear();
                cellId.clear();
            }       
            std::vector<double> perfPressure;
            std::vector<double> cellPressure;
            std::vector<BlackoilDefs::CompVec>  massRate; // (surface volumes)
            std::vector<int> cellId;
        protected:
            WellReport() {}
        private:
            static WellReport* well_report_;
        };

    private:
	// Use the Peaceman well model to compute well indices
	double computeWellIndex(double radius, const Dune::FieldVector<double, 3>& cubical,
				const Opm::Rock<3>::PermTensor& permeability, double skin_factor) const;

	struct WellData { WellType type; WellControl control; double target; double reference_bhp_depth; };
        std::vector<WellData> well_data_;
        struct PerfData { int cell; double well_index; };
        Opm::SparseTable<PerfData> perf_data_;
        std::vector<double> well_cell_flux_;
        std::vector<double> well_cell_pressure_;
        Dune::FieldVector<double, 3> injection_mixture_;
	std::vector<std::string> well_names_;
    };

    BlackoilWells::WellReport* BlackoilWells::WellReport::well_report_ = 0x0;

    // ------------ Method implementations --------------


    // Forward declaration of some helper functions.
    namespace
    {
        int prod_control_mode(const std::string& control);
        int inje_control_mode(const std::string& control);
        template<class grid_t>
        const Dune::FieldVector<double,3> getCubeDim(const grid_t& grid, int cell);

    } // anon namespace

    inline void BlackoilWells::init(const Opm::EclipseGridParser& parser,
				    const Dune::CpGrid& grid,
				    const Opm::Rock<3>& rock) 
    {
	std::vector<std::string> keywords;
	keywords.push_back("WELSPECS");
	keywords.push_back("COMPDAT");
// 	keywords.push_back("WELTARG");
	if (!parser.hasFields(keywords)) {
	    MESSAGE("Missing well keywords in deck, initializing no wells.");
            return;
	}
	if (!(parser.hasField("WCONINJE") || parser.hasField("WCONPROD")) ) {
	    THROW("Needed field is missing in file");
	}
        using namespace Opm;

	// Get WELLSPECS data
	const WELSPECS& welspecs = parser.getWELSPECS();
	const int num_welspecs   = welspecs.welspecs.size();
	well_names_.reserve(num_welspecs);
	well_data_.reserve(num_welspecs);
	for (int i=0; i<num_welspecs; ++i) {
	    well_names_.push_back(welspecs.welspecs[i].name_);
	    WellData wd;
	    well_data_.push_back(wd);
            well_data_.back().reference_bhp_depth = welspecs.welspecs[i].datum_depth_BHP_;
            if (welspecs.welspecs[i].datum_depth_BHP_ < 0.0) {
                // Set refdepth to a marker value, will be changed
                // after getting perforation data to the centroid of
                // the cell of the top well perforation.
                well_data_.back().reference_bhp_depth = -1e100;
            }
        }

	// Get COMPDAT data   
	const COMPDAT& compdats = parser.getCOMPDAT();
	const int num_compdats  = compdats.compdat.size();
        std::vector<std::vector<PerfData> > wellperf_data(num_welspecs);
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
		    int kz1 = compdats.compdat[kw].grid_ind_[2] - 1;
                    int kz2 = compdats.compdat[kw].grid_ind_[3] - 1;
                    for (int kz = kz1; kz <= kz2; ++kz) {
                        int cart_grid_indx = ix + cpgdim[0]*(jy + cpgdim[1]*kz);
                        std::map<int, int>::const_iterator cgit = 
                            cartesian_to_compressed.find(cart_grid_indx);
                        if (cgit == cartesian_to_compressed.end()) {
                            THROW("Cell with i,j,k indices " << ix << ' ' << jy << ' '
                                  << kz << " not found!");
                        }
                        int cell = cgit->second;
                        PerfData pd;
                        pd.cell = cell;
                        if (compdats.compdat[kw].connect_trans_fac_ > 0.0) {
                            pd.well_index = compdats.compdat[kw].connect_trans_fac_;
                        } else {
                            double radius = 0.5*compdats.compdat[kw].diameter_;
                            if (radius <= 0.0) {
                                radius = 0.5*unit::feet;
                                MESSAGE("Warning: Well bore internal radius set to " << radius);
                            }
                            Dune::FieldVector<double, 3> cubical = getCubeDim(grid, cell);
                            const Rock<3>::PermTensor permeability = rock.permeability(cell);  
                            pd.well_index = computeWellIndex(radius, cubical, permeability,
                                                             compdats.compdat[kw].skin_factor_);
                        }
                        wellperf_data[wix].push_back(pd);
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
        for (int w = 0; w < num_welspecs; ++w) {
            perf_data_.appendRow(wellperf_data[w].begin(), wellperf_data[w].end());
            if (well_data_[w].reference_bhp_depth == -1e100) {
                // It was defaulted. Set reference depth to minimum perforation depth.
                double min_depth = 1e100;
                int num_wperfs = wellperf_data[w].size();
                for (int perf = 0; perf < num_wperfs; ++perf) {
                    double depth = grid.cellCentroid(wellperf_data[w][perf].cell)[2];
                    min_depth = std::min(min_depth, depth);
                }
                well_data_[w].reference_bhp_depth = min_depth;
            }
        }
 
	// Get WCONINJE data
        injection_mixture_ = 0.0;
        if (parser.hasField("WCONINJE")) {
            const WCONINJE& wconinjes = parser.getWCONINJE();
            const int num_wconinjes = wconinjes.wconinje.size();
            int injector_component = -1;
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
                            well_data_[wix].target = wconinjes.wconinje[kw].surface_flow_max_rate_;
                            break;
                        case 1:  // RESV
                            well_data_[wix].control = Rate;
                            well_data_[wix].target = wconinjes.wconinje[kw].reservoir_flow_max_rate_;
                            break;
                        case 2:  // BHP
                            well_data_[wix].control = Pressure;
                            well_data_[wix].target = wconinjes.wconinje[kw].BHP_limit_;
                            break;
                        case 3:  // THP
                            well_data_[wix].control = Pressure;
                            well_data_[wix].target = wconinjes.wconinje[kw].THP_limit_;
                            break;
                        default:
                            THROW("Unknown well control mode; WCONIJE  = "
                                  << wconinjes.wconinje[kw].control_mode_
                                  << " in input file");
                        }
                        int itp = -1;
                        if (wconinjes.wconinje[kw].injector_type_ == "WATER") {
                            itp = Water;
                        } else if (wconinjes.wconinje[kw].injector_type_ == "OIL") {
                            itp = Oil;
                        } else if (wconinjes.wconinje[kw].injector_type_ == "GAS") {
                            itp = Gas;
                        }
                        if (itp == -1 || (injector_component != -1 && itp != injector_component)) {
                            if (itp == -1) {
                                THROW("Error in injector specification, found no known fluid type.");
                            } else {
                                THROW("Error in injector specification, we can only handle a single injection fluid.");
                            }
                        } else {
                            injector_component = itp;
                        }
                    }
                }
                if (!well_found) {
                    THROW("Undefined well name: " << wconinjes.wconinje[kw].well_
                          << " in WCONINJE");
                }
            }
            if (injector_component != -1) {
                injection_mixture_[injector_component] = 1.0;
            }
        } else {
            // No WCONINJE.
            // This default is only invoked if production wells
            // start injecting, (and only if there are no injection wells).
            // In other words, only if we have a primary production type scenario.
            // In this case we expect the flow to be very small, and not affect
            // the simulation in any significant way.
            injection_mixture_[Oil] = 1.0;
        }

	// Get WCONPROD data
        if (parser.hasField("WCONPROD")) {
            const WCONPROD& wconprods = parser.getWCONPROD();
            const int num_wconprods   = wconprods.wconprod.size();
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
                            well_data_[wix].target = wconprods.wconprod[kw].reservoir_flow_max_rate_;
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
        }

	// Get WELTARG data
        if (parser.hasField("WELTARG")) {
            const WELTARG& weltargs = parser.getWELTARG();
            const int num_weltargs  = weltargs.weltarg.size();	
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
        }

        // Debug output.
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
			  << perf_data_[i][j].well_index << std::endl;
	    }
	}

        // Ensuring that they have the right size.
        well_cell_pressure_.resize(grid.numCells(), -1e100);
        well_cell_flux_.resize(grid.numCells(), 0.0);
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

    inline double BlackoilWells::referenceDepth(int wellnum) const
    {
        return well_data_[wellnum].reference_bhp_depth;
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

    inline void BlackoilWells::update(int num_cells,
                                      const std::vector<double>& well_perf_pressures,
                                      const std::vector<double>& well_perf_fluxes)
    {
        // Input is per perforation, data members store for all cells.
        ASSERT(perf_data_.dataSize() == int(well_perf_pressures.size()));
        well_cell_pressure_.resize(num_cells, -1e100);
        well_cell_flux_.resize(num_cells, 0.0);
        int pcount = 0;
        for (int w = 0; w < numWells(); ++w) {
            for (int perf = 0; perf < numPerforations(w); ++perf) {
                int cell = wellCell(w, perf);
                well_cell_pressure_[cell] = well_perf_pressures[pcount];
                well_cell_flux_[cell] = well_perf_fluxes[pcount];
                ++pcount;
            }
        }
        ASSERT(pcount == perf_data_.dataSize());
    }

    inline double BlackoilWells::wellToReservoirFlux(int cell) const
    {
        return well_cell_flux_[cell];
    }

    inline double BlackoilWells::perforationPressure(int cell) const
    {
        return well_cell_pressure_[cell];
    }

    inline Dune::FieldVector<double, 3> BlackoilWells::injectionMixture(int cell) const
    {
        return injection_mixture_;
    }

    inline double BlackoilWells::computeWellIndex(double radius,
						  const Dune::FieldVector<double, 3>& cubical,
						  const Opm::Rock<3>::PermTensor& permeability,
						  double skin_factor) const
    {
	// Use the Peaceman well model to compute well indices.
	// radius is the radius of the well.
	// cubical contains [dx, dy, dz] of the cell.
	// (Note that the well model asumes that each cell is a cuboid).
	// permeability is the permeability of the given cell.
	// returns the well index of the cell.

	// sse: Using the Peaceman modell.
	// NOTE: The formula is valid for cartesian grids, so the result can be a bit
	// (in worst case: there is no upper bound for the error) off the mark.
	double effective_perm = sqrt(permeability(0,0) * permeability(1,1));
	// sse: The formula for r_0 can be found on page 39 of
	// "Well Models for Mimetic Finite Differerence Methods and Improved Representation
	//  of Wells in Multiscale Methods" by Ingeborg Skjelkvåle Ligaarden.
	assert(permeability(0,0) > 0.0);
	assert(permeability(1,1) > 0.0);
	double kxoy = permeability(0,0) / permeability(1,1);
	double kyox = permeability(1,1) / permeability(0,0);
	double r0_denominator = pow(kyox, 0.25) + pow(kxoy, 0.25);
	double r0_numerator = sqrt((sqrt(kyox)*cubical[0]*cubical[0]) +
				   (sqrt(kxoy)*cubical[1]*cubical[1]));
	assert(r0_denominator > 0.0);
	double r0 = 0.28 * r0_numerator / r0_denominator;
	assert(radius > 0.0);
	assert(r0 > 0.0);
	if (r0 < radius) {
	    std::cout << "ERROR: Too big well radius detected.";
	    std::cout << "Specified well radius is " << radius
		      << " while r0 is " << r0 << ".\n";
	}

        const double two_pi = 6.2831853071795864769252867665590057683943387987502116419498;

	double wi_denominator = log(r0 / radius) + skin_factor;
	double wi_numerator = two_pi * cubical[2];
	assert(wi_denominator > 0.0);
	double wi = effective_perm * wi_numerator / wi_denominator;
	assert(wi > 0.0);
	return wi;
    }



    // ------------- Helper functions for init() --------------

    namespace
    {

        int prod_control_mode(const std::string& control){
            const int num_prod_control_modes = 8;
            static std::string prod_control_modes[num_prod_control_modes] =
                {std::string("ORAT"), std::string("WRAT"), std::string("GRAT"),
                 std::string("LRAT"), std::string("RESV"), std::string("BHP"),
                 std::string("THP"), std::string("GRUP") };
            int m = -1;
            for (int i=0; i<num_prod_control_modes; ++i) {
                if (control == prod_control_modes[i]) {
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
            const int num_inje_control_modes = 5;
            static std::string inje_control_modes[num_inje_control_modes] =
                {std::string("RATE"), std::string("RESV"), std::string("BHP"),
                 std::string("THP"), std::string("GRUP") };
            int m = -1;
            for (int i=0; i<num_inje_control_modes; ++i) {
                if (control == inje_control_modes[i]) {
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


        template<class grid_t>
        const Dune::FieldVector<double,3> getCubeDim(const grid_t& grid, int cell)
        {
            Dune::FieldVector<double, 3> cube;
            int num_local_faces = grid.numCellFaces(cell);
            std::vector<double> x(num_local_faces);
            std::vector<double> y(num_local_faces);
            std::vector<double> z(num_local_faces);
            for (int lf=0; lf<num_local_faces; ++ lf) {
                int face = grid.cellFace(cell,lf);
                const Dune::FieldVector<double,3>& centroid = 
                    grid.faceCentroid(face);
                x[lf] = centroid[0];
                y[lf] = centroid[1];
                z[lf] = centroid[2];
            }
            cube[0] = *max_element(x.begin(), x.end()) -
                *min_element(x.begin(), x.end());
            cube[1] = *max_element(y.begin(), y.end()) -
                *min_element(y.begin(), y.end());
            cube[2] = *max_element(z.begin(), z.end()) -
                *min_element(z.begin(), z.end());
            return cube;
        }

    } // anon namespace

} // namespace Opm


#endif // OPM_BLACKOILWELLS_HEADER_INCLUDED
