//==============================================================================
//!
//! \file upscale_elasticity.cpp
//!
//! \date Nov 9 2011
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Elasticity upscaling on cornerpoint grids
//!
//==============================================================================
#ifdef HAVE_CONFIG_H
# include "config.h"     
#endif

#include <iostream>
#include <unistd.h>
#include <cstring>
#include <dune/common/exceptions.hh> // We use exceptions
#include <dune/common/version.hh>
#include <opm/core/utility/StopWatch.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#if HAVE_OPENMP
#include <omp.h>
#endif

#include <dune/common/version.hh>
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 3)
#include <dune/common/parallel/mpihelper.hh>
#else
#include <dune/common/mpihelper.hh>
#endif

#include <opm/elasticity/elasticity_upscale.hpp>
#include <opm/elasticity/matrixops.hpp>


//! \brief Display the available command line parameters
void syntax(char** argv)
{
  std::cerr << "Usage: " << argv[0] << " gridfilename=filename.grdecl [method=]" << std::endl
            << "\t[xmax=] [ymax=] [zmax=] [xmin=] [ymin=] [zmin=] [linsolver_type=]" << std::endl
            <<" \t[Emin=] [ctol=] [ltol=] [rock_list=] [vtufilename=] [output=] [linsolver_verbose=]" << std::endl << std::endl
            << "\t gridfilename             - the grid file. can be 'uniform'" << std::endl
            << "\t vtufilename              - save results to vtu file" << std::endl
            << "\t rock_list                - use material information from given rocklist" << std::endl
            << "\t output                   - output results in text report format" << std::endl
            << "\t method                   - can be MPC or Mortar" << std::endl
            << "\t\t if not specified, Mortar couplings are used" << std::endl
            << "\t (x|y|z)max/(x|y|z)min - maximum/minimum grid coordinates." << std::endl
            << "\t\t if not specified, coordinates are found by grid traversal" << std::endl
            << "\t cells(x|y|z)             - number of cells if gridfilename=uniform." << std::endl
            << "\t lambda(x|y)              - order of lagrangian multipliers in first/second direction" << std::endl
            << "\t Emin                     - minimum E modulus (to avoid numerical issues)" << std::endl
            << "\t ctol                     - collapse tolerance in grid parsing" << std::endl
            << "\t ltol                     - tolerance in iterative linear solvers" << std::endl
            << "\t linsolver_type=iterative - use a suitable iterative method (cg or gmres)" << std::endl
            << "\t linsolver_type=direct    - use the SuperLU sparse direct solver" << std::endl
            << "\t verbose                  - set to true to get verbose output" << std::endl
            << "\t fastamg                  - use the fast AMG" << std::endl
            << "\t linsolver_restart        - number of iterations before gmres is restarted" << std::endl
            << "\t linsolver_presteps       - number of pre-smooth steps in the AMG" << std::endl
            << "\t linsolver_poststeps      - number of post-smooth steps in the AMG" << std::endl
            << "\t\t affects memory usage" << std::endl
            << "\t linsolver_symmetric      - use symmetric linear solver. Defaults to true" << std::endl
            << "\t mortar_precond           - preconditioner for mortar block. Defaults to schur-amg" << std::endl;
}


enum UpscaleMethod {
  UPSCALE_NONE   = 0,
  UPSCALE_MPC    = 1,
  UPSCALE_MORTAR = 2
};

//! \brief Structure holding parameters configurable from command line
struct Params {
  //! \brief The eclipse grid file
  std::string file;
  //! \brief Rocklist overriding material params in .grdecl
  std::string rocklist;
  //! \brief VTU output file
  std::string vtufile;
  //! \brief Text output file
  std::string output;
  //! \brief Method
  UpscaleMethod method;
  //! \brief A scaling parameter for the E-modulus to avoid numerical issues
  double Emin;
  //! \brief The tolerance for collapsing nodes in the Z direction
  double ctol;
  //! \brief Minimum grid coordinates
  double min[3];
  //! \brief Maximum grid coordinates
  double max[3];
  //! \brief Polynomial order of multipliers in first / second direction
  int lambda[2];
  //! \brief Number of elements on interface grid in the x direction
  int n1;
  //! \brief Number of elements on interface grid in the y direction
  int n2;
  Opm::Elasticity::LinSolParams linsolver;

  // Debugging options

  //! \brief cell in x
  int cellsx;
  //! \brief cell in y
  int cellsy;
  //! \brief cell in z
  int cellsz;
  //! \brief verbose output
  bool verbose;
};

//! \brief Parse the command line arguments
void parseCommandLine(int argc, char** argv, Params& p)
{
  Opm::parameter::ParameterGroup param(argc, argv);
  p.max[0]    = param.getDefault("xmax",-1);
  p.max[1]    = param.getDefault("ymax",-1);
  p.max[2]    = param.getDefault("zmax",-1);
  p.min[0]    = param.getDefault("xmin",-1);
  p.min[1]    = param.getDefault("ymin",-1);
  p.min[2]    = param.getDefault("zmin",-1);
  p.lambda[0] = param.getDefault("lambdax", 2);
  p.lambda[1] = param.getDefault("lambday", 2);
  std::string method = param.getDefault<std::string>("method","mortar");
  if (method == "mpc")
    p.method = UPSCALE_MPC;
  if (method == "mortar")
    p.method = UPSCALE_MORTAR;
  if (method == "none")
    p.method = UPSCALE_NONE;
  p.Emin     = param.getDefault<double>("Emin",0.0);
  p.ctol     = param.getDefault<double>("ctol",1.e-6);
  p.file     = param.get<std::string>("gridfilename");
  p.rocklist = param.getDefault<std::string>("rock_list","");
  p.vtufile  = param.getDefault<std::string>("vtufilename","");
  p.output   = param.getDefault<std::string>("output","");
  p.verbose  = param.getDefault<bool>("verbose",false);
  size_t i;
  if ((i=p.vtufile.find(".vtu")) != std::string::npos)
    p.vtufile = p.vtufile.substr(0,i);

  if (p.file == "uniform") {
    p.cellsx   = param.getDefault("cellsx",3);
    p.cellsy   = param.getDefault("cellsy",3);
    p.cellsz   = param.getDefault("cellsz",3);
  }
  p.linsolver.parse(param);
  p.n1       = -1;
  p.n2       = -1;
}

//! \brief Write a log of the simulation to a text file
void writeOutput(const Params& p, Opm::time::StopWatch& watch, int cells,
                 const std::vector<double>& volume, 
                 const Dune::FieldMatrix<double,6,6>& C)
{
  // get current time
  time_t rawtime;
  struct tm* timeinfo;
  time(&rawtime);
  timeinfo = localtime(&rawtime);

  // get hostname
  char hostname[1024];
  gethostname(hostname,1024);

  std::string method = "mortar";
  if (p.method == UPSCALE_MPC)
    method = "mpc";
  if (p.method == UPSCALE_NONE)
    method = "none";

  // write log
  std::ofstream f;
  f.open(p.output.c_str());
  f << "######################################################################" << std::endl
    << "# Results from upscaling elastic moduli." << std::endl
    << "#" << std::endl
    << "# Finished: " << asctime(timeinfo)
    << "# Hostname: " << hostname << std::endl
    << "#" << std::endl
    << "# Upscaling time: " << watch.secsSinceStart() << " secs" << std::endl
    << "#" << std::endl;
  if (p.file == "uniform") {
    f << "# Uniform grid used" << std::endl
      << "#\t cells: " << p.cellsx*p.cellsy*p.cellsz << std::endl;
  }
  else {
    f  << "# Eclipse file: " << p.file << std::endl
       << "#\t cells: " << cells << std::endl;
  }
  f << "#" << std::endl;
  if (!p.rocklist.empty()) {
    f << "# Rock list: " << p.rocklist << std::endl
      << "#" << std::endl;
  }
  f << "# Options used:" << std::endl
    << "#\t         method: " << method << std::endl
    << "#\t linsolver_type: " << (p.linsolver.type==Opm::Elasticity::DIRECT?"direct":"iterative")
                              << std::endl;
  if (p.linsolver.type == Opm::Elasticity::ITERATIVE)
    f << "#\t           ltol: " << p.linsolver.tol << std::endl;
  if (p.file == "uniform") {
    f << "#\t          cellsx: " << p.cellsx << std::endl
      << "#\t          cellsy: " << p.cellsy << std::endl
      << "#\t          cellsz: " << p.cellsz << std::endl;
  }
  f << "#" << std::endl
    <<"# Materials: " << volume.size() << std::endl;
  for (size_t i=0;i<volume.size();++i)
    f << "#\t Material" << i+1 << ": " << volume[i]*100 << "%" << std::endl;
  f << "#" << std::endl
    << "######################################################################" << std::endl
    << C << std::endl;
}

//! \brief Main solution loop. Allows templating over the AMG type
  template<class GridType, class AMG>
int run(Params& p)
{
  try {
    static const int dim = 3;

    Opm::time::StopWatch watch;
    watch.start();

    GridType grid;
    if (p.file == "uniform") {
      std::array<int,3> cells;
      cells[0] = p.cellsx;
      cells[1] = p.cellsy;
      cells[2] = p.cellsz;
      std::array<double,3> cellsize;
      cellsize[0] = cellsize[1] = cellsize[2] = 1.f; 
      grid.createCartesian(cells,cellsize);
    } else
      grid.readEclipseFormat(p.file,p.ctol,false);

    typedef typename GridType::ctype ctype;
    Opm::Elasticity::ElasticityUpscale<GridType, AMG> upscale(grid, p.ctol, 
                                                              p.Emin, p.file,
                                                              p.rocklist,
                                                              p.verbose);
    if (p.max[0] < 0 || p.min[0] < 0) {
      std::cout << "determine side coordinates..." << std::endl;
      upscale.findBoundaries(p.min,p.max);
      std::cout << "  min " << p.min[0] << " " << p.min[1] << " " << p.min[2] << std::endl;
      std::cout << "  max " << p.max[0] << " " << p.max[1] << " " << p.max[2] << std::endl;
    }
    if (p.n1 == -1 || p.n2 == -1) {
      p.n1 = grid.logicalCartesianSize()[0];
      p.n2 = grid.logicalCartesianSize()[1];
    }
    if (p.linsolver.zcells == -1) {
      double lz = p.max[2]-p.min[2];
      int nz = grid.logicalCartesianSize()[2];
      double hz = lz/nz;
      double lp = sqrt((double)(p.max[0]-p.min[0])*(p.max[1]-p.min[1]));
      int np = std::max(grid.logicalCartesianSize()[0],
                        grid.logicalCartesianSize()[1]);
      double hp = lp/np;
      p.linsolver.zcells = (int)(2*hp/hz+0.5);
    }
    std::cout << "logical dimension: " << grid.logicalCartesianSize()[0]
              << "x"                   << grid.logicalCartesianSize()[1]
              << "x"                   << grid.logicalCartesianSize()[2]
              << std::endl;

    if (p.method == UPSCALE_MPC) {
      std::cout << "using MPC couplings in all directions..." << std::endl;
      upscale.periodicBCs(p.min, p.max);
      std::cout << "preprocessing grid..." << std::endl;
      upscale.A.initForAssembly();
    } else if (p.method == UPSCALE_MORTAR) {
      std::cout << "using Mortar couplings.." << std::endl;
      upscale.periodicBCsMortar(p.min, p.max, p.n1, p.n2,
                                p.lambda[0], p.lambda[1]);
    } else if (p.method == UPSCALE_NONE) {
      std::cout << "no periodicity approach applied.." << std::endl;
      upscale.fixCorners(p.min, p.max);
      upscale.A.initForAssembly();
    }
    Dune::FieldMatrix<double,6,6> C;
    Dune::VTKWriter<typename GridType::LeafGridView>* vtkwriter=0;
    if (!p.vtufile.empty())
      vtkwriter = new Dune::VTKWriter<typename GridType::LeafGridView>(grid.leafView());
    Opm::Elasticity::Vector field[6];
    std::cout << "assembling elasticity operator..." << "\n";
    upscale.assemble(-1,true);
    std::cout << "setting up linear solver..." << std::endl;
    upscale.setupSolvers(p.linsolver);

//#pragma omp parallel for schedule(static)
    for (int i=0;i<6;++i) {
      std::cout << "processing case " << i+1 << "..." << std::endl;
      std::cout << "\tassembling load vector..." << std::endl;
      upscale.assemble(i,false);
      std::cout << "\tsolving..." << std::endl;
      upscale.solve(i);
      upscale.A.expandSolution(field[i],upscale.u[i]);
#define CLAMP(x) (fabs(x)<1.e-4?0.0:x)
      for (size_t j=0;j<field[i].size();++j) {
        double val = field[i][j];
        field[i][j] = CLAMP(val);
      }
      Dune::FieldVector<double,6> v;
      upscale.averageStress(v,upscale.u[i],i);
      for (int j=0;j<6;++j)
        C[i][j] = CLAMP(v[j]);
    }
    for (int i=0;i<6;++i) {
      std::stringstream str;
      str << "sol " << i+1;
      if (vtkwriter)
        vtkwriter->addVertexData(field[i], str.str().c_str(), dim);
    }
    if (vtkwriter) {
      vtkwriter->write(p.vtufile);
      delete vtkwriter;
    }
    // voigt notation
    for (int j=0;j<6;++j)
      std::swap(C[3][j],C[5][j]);
    for (int j=0;j<6;++j)
      std::swap(C[j][3],C[j][5]);
    std::cout << "---------" << std::endl;
    std::cout << C << std::endl;
    if (!p.output.empty()) {
      writeOutput(p,watch,grid.size(0),upscale.volumeFractions,C);
    }

    return 0;
  }
  catch (Dune::Exception &e) {
    throw e;
  }
  catch (...) {
    throw;
  }
  return 1;
}

//! \brief Main driver
int main(int argc, char** argv)
try
{
  try {
    if (argc < 2 || strcmp(argv[1],"-h") == 0 
                 || strcmp(argv[1],"--help") == 0
                 || strcmp(argv[1],"-?") == 0) {
      syntax(argv);
      exit(1);
    }

    Dune::MPIHelper& mpi=Dune::MPIHelper::instance(argc, argv);
    const int size = mpi.size();
    if (size != 1) {
      std::cerr << "This program does not support MPI parallelization" << std::endl;
      return 2;
    }

    Params p;
    parseCommandLine(argc,argv,p);

    if (p.linsolver.fastamg)
      return run<Dune::CpGrid, Opm::Elasticity::FastAMG>(p);
    else
      return run<Dune::CpGrid, Opm::Elasticity::AMG>(p);
  }
  catch (Dune::Exception &e) {
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...) {
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
  return 1;
}
catch (const std::exception &e) {
  std::cerr << "Program threw an exception: " << e.what() << "\n";
  throw;
}
