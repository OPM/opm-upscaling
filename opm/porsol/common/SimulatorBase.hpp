//===========================================================================
//
// File: SimulatorBase.hpp
//
// Created: Tue Aug 11 15:01:48 2009
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//            Bård Skaflestad     <bard.skaflestad@sintef.no>
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2009, 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2009, 2010 Statoil ASA.

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

#ifndef OPENRS_SIMULATORBASE_HEADER
#define OPENRS_SIMULATORBASE_HEADER


#include <opm/common/utility/parameters/ParameterGroup.hpp>

#include <opm/common/utility/numeric/SparseVector.hpp>
#include <opm/grid/utility/SparseTable.hpp>
#include <opm/input/eclipse/Units/Units.hpp>

#include <opm/grid/common/Volumes.hpp>
#include <opm/grid/CpGrid.hpp>

#include <opm/porsol/common/GridInterfaceEuler.hpp>
#include <opm/porsol/common/ReservoirPropertyCapillary.hpp>
#include <opm/porsol/common/BoundaryConditions.hpp>
#include <opm/porsol/common/setupGridAndProps.hpp>
#include <opm/porsol/common/setupBoundaryConditions.hpp>
#include <opm/porsol/common/SimulatorUtilities.hpp>

#include <opm/porsol/euler/EulerUpstream.hpp>
#include <opm/porsol/euler/ImplicitCapillarity.hpp>

#include <opm/porsol/mimetic/MimeticIPEvaluator.hpp>
#include <opm/porsol/mimetic/IncompFlowSolverHybrid.hpp>


#include <opm/common/utility/platform_dependent/disable_warnings.h>

#include <dune/grid/yaspgrid.hh>

#include <opm/common/utility/platform_dependent/reenable_warnings.h>


#include <fstream>
#include <iterator>
#include <iostream>


namespace Opm
{




    /// @brief
    /// @todo Doc me!
    /// @tparam
    template <class SimTraits>
    class SimulatorBase
    {
    public:

	/// @brief
	/// @todo Doc me!
	SimulatorBase()
	    : simulation_steps_(1),
	      stepsize_(1.0),   // init() expects units of days! Yes, but now the meaning of stepsize_ changes
	                        // from days (here) to seconds (after init()). Solution to that?
              residual_tolerance_(1e-8),
              linsolver_verbosity_(1),
              linsolver_type_(1)
	{
	}

	/// @brief Initialization from parameters.
	/// @param param a parameter object
	void init(const Opm::ParameterGroup& param)
	{
	    initControl(param);
	    initGridAndProps(param);
	    initInitialConditions(param);
	    initBoundaryConditions(param);
            initSources(param);
	    initSolvers(param);

	    // Write any unused parameters.
	    std::cout << "====================   Unused parameters:   ====================\n";
	    param.displayUsage();
	    std::cout << "================================================================\n";
	}

    protected:
	typedef Dune::CpGrid                                         GridType;
 	enum { Dimension = GridType::dimension };
	typedef Dune::FieldVector<double, Dimension>                 Vector;
 	typedef typename SimTraits::template ResProp<Dimension>::Type ResProp;
	typedef GridInterfaceEuler<GridType>                   GridInterface;
	typedef GridInterface::CellIterator                    CellIter;
	typedef CellIter::FaceIterator                         FaceIter;
	typedef BasicBoundaryConditions<true, true>                 BCs;
	typedef typename SimTraits::template FlowSolver<GridInterface, BCs>::Type FlowSolver;
        typedef typename SimTraits::template TransportSolver<GridInterface, BCs>::Type TransportSolver;

	int simulation_steps_;
	double stepsize_;
        std::vector<double> init_saturation_;
        Vector gravity_;
	double residual_tolerance_;
	int linsolver_verbosity_;
	int linsolver_type_;

	GridType grid_;
	GridInterface ginterf_;
	ResProp res_prop_;
	BCs bcond_;
	Opm::SparseVector<double> injection_rates_;
        std::vector<double> injection_rates_psolver_;	// Should modify psolver to take SparseVector
        FlowSolver flow_solver_;
	TransportSolver transport_solver_;


	virtual void initControl(const Opm::ParameterGroup& param)
	{
	    simulation_steps_ = param.getDefault("simulation_steps", simulation_steps_);
	    stepsize_ = Opm::unit::convert::from(param.getDefault("stepsize", stepsize_),
                                                  Opm::unit::day);
	}

	virtual void initGridAndProps(const Opm::ParameterGroup& param)
	{
	    setupGridAndProps(param, grid_, res_prop_);
	    ginterf_.init(grid_);

            gravity_[0] = param.getDefault("gx", 0.0);
            gravity_[1] = param.getDefault("gy", 0.0);
            gravity_[2] = param.getDefault("gz", 0.0); //Dune::unit::gravity);
	}

	virtual void initInitialConditions(const Opm::ParameterGroup& param)
	{
            if (param.getDefault("init_saturation_from_file", false)) {
                std::string filename = param.get<std::string>("init_saturation_filename");
                std::ifstream satfile(filename.c_str());
                if (!satfile) {
                    OPM_THROW(std::runtime_error, "Could not open initial saturation file: " << filename);
                }
                int num_sats;
                satfile >> num_sats;
                if (num_sats != ginterf_.numberOfCells()) {
                    OPM_THROW(std::runtime_error, "Number of saturation values claimed different from number of grid cells: "
                          << num_sats << " vs. " << ginterf_.numberOfCells());
                }
                std::istream_iterator<double> beg(satfile);
                std::istream_iterator<double> end;
                init_saturation_.assign(beg, end);
                if (int(init_saturation_.size()) != num_sats) {
                    OPM_THROW(std::runtime_error, "Number of saturation values claimed different from actual file content: "
                          << num_sats << " vs. " << init_saturation_.size());
                }
            } else {
                double init_s = param.getDefault("init_saturation", 0.0);
                init_saturation_.clear();
                init_saturation_.resize(ginterf_.numberOfCells(), init_s);
            }
	}

	virtual void initBoundaryConditions(const Opm::ParameterGroup& param)
	{
	    setupBoundaryConditions(param, ginterf_, bcond_);
	}

        virtual void initSources(const Opm::ParameterGroup& /* param */)
        {
            int nc = ginterf_.numberOfCells();
	    injection_rates_ = Opm::SparseVector<double>(nc);
	    injection_rates_psolver_.resize(nc, 0.0);
//             injection_rates_.addElement(1.0, 0);
//             injection_rates_.addElement(-1.0, nc - 1);
//             injection_rates_psolver_[0] = 1.0;
//             injection_rates_psolver_[nc - 1] = -1.0;
        }

	virtual void initSolvers(const Opm::ParameterGroup& param)
	{
	    // Initialize flow solver.
	    flow_solver_.init(ginterf_, res_prop_, gravity_, bcond_);
            residual_tolerance_ = param.getDefault("residual_tolerance", residual_tolerance_);
            linsolver_verbosity_ = param.getDefault("linsolver_verbosity", linsolver_verbosity_);
            linsolver_type_ = param.getDefault("linsolver_type", linsolver_type_);
	    //flow_solver_.assembleStatic(ginterf_, res_prop_);
	    // Initialize transport solver.
	    transport_solver_.init(param, ginterf_, res_prop_, bcond_);
	}


    };



} // namespace Opm



#endif // OPENRS_SIMULATORBASE_HEADER
