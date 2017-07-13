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

#include <opm/common/utility/platform_dependent/disable_warnings.h>

#include <dune/common/exceptions.hh> // We use exceptions
#include <dune/common/version.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/matrixmarket.hh>
#include <dune/common/fmatrix.hh>

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include <dune/common/version.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Parser/ParseContext.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>

#include <opm/common/utility/platform_dependent/reenable_warnings.h>

#include <opm/core/utility/StopWatch.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>

#include <opm/elasticity/elasticity_upscale.hpp>
#include <opm/elasticity/matrixops.hpp>

#include <cstring>
#include <iostream>

#include <unistd.h>

using namespace Opm::Elasticity;

#ifndef M_PI
#define M_PI std::acos(-1.0)
#endif

//! \brief Display the available command line parameters
void syntax(char** argv)
{
  std::cerr << "Usage: " << argv[0] << " gridfilename=filename.grdecl [method=]" << std::endl
            << "\t[xmax=] [ymax=] [zmax=] [xmin=] [ymin=] [zmin=] [linsolver_type=]" << std::endl
            <<" \t[Emin=] [ctol=] [ltol=] [rock_list=] [vtufilename=] [output=] [verbose=]" << std::endl << std::endl
            << "\t gridfilename             - the grid file. can be 'uniform'" << std::endl
            << "\t vtufilename              - save results to vtu file" << std::endl
            << "\t rock_list                - use material information from given rocklist" << std::endl
            << "\t output                   - output results in text report format" << std::endl
            << "\t resultfilename           - result template filename" << std::endl
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
            << "\t linsolver_type=direct    - use the SuperLU or UMFPACK sparse direct solvers" << std::endl
            << "\t verbose                  - set to true to get verbose output" << std::endl
            << "\t linsolver_pre            - preconditioner for elasticity block. amg, fastamg or schwarz" << std::endl
            << "\t linsolver_restart        - number of iterations before gmres is restarted" << std::endl
            << "\t linsolver_presteps       - number of pre-smooth steps in the AMG" << std::endl
            << "\t linsolver_poststeps      - number of post-smooth steps in the AMG" << std::endl
            << "\t linsolver_smoother       - smoother used in the AMG" << std::endl
            << "\t linsolver_report         - print report at end of solution phase" << std::endl
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
  //! \brief Result filenames
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
  //! \brief Dip angle for wavespeeds
  double dip;
  //! \brief Azimuth angle for wavespeeds
  double azimuth;

  //! \brief Linear solver parameters
  LinSolParams linsolver;

  // Debugging options

  //! \brief cell in x
  int cellsx;
  //! \brief cell in y
  int cellsy;
  //! \brief cell in z
  int cellsz;
  //! \brief verbose output
  bool verbose;
  //! \brief Run a inspection only, currently 'mesh, results, load'
  std::string inspect;
  //! \brief Result template filename (input/output)
  std::string resultfilename;
};

//! \brief Parse the command line arguments
void parseCommandLine(int argc, char** argv, Params& p)
{
  Opm::ParameterGroup param(argc, argv);
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
  p.dip      = param.getDefault<double>("dip_angle", M_PI/2);
  p.azimuth  = param.getDefault<double>("azimuth_angle", 0.0);
  p.ctol     = param.getDefault<double>("ctol",1.e-6);
#ifndef HAVE_OLD_CPGRID_API
  if (p.ctol != 1e-6)
    std::cerr << "WARNING: Use PINCH in input file to control nodal collapse!" << std::endl;
#endif
  p.file     = param.get<std::string>("gridfilename");
  p.rocklist = param.getDefault<std::string>("rock_list","");
  p.vtufile  = param.getDefault<std::string>("vtufilename","");
  p.resultfilename  = param.getDefault<std::string>("resultfilename","");
  p.output   = param.getDefault<std::string>("output","");
  p.verbose  = param.getDefault<bool>("verbose",false);
  p.inspect  = param.getDefault<std::string>("inspect","");
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
                 const std::vector<double>& volume, bool bySat,
                 const Dune::FieldMatrix<double,6,6>& C,
                 double upscaledRho, const Dune::FieldVector<double,3>& speeds)
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
    << "#\t linsolver_type: " << (p.linsolver.type==DIRECT?"direct":"iterative")
                              << std::endl;
  if (p.linsolver.type == ITERATIVE)
    f << "#\t           ltol: " << p.linsolver.tol << std::endl;
  if (p.file == "uniform") {
    f << "#\t          cellsx: " << p.cellsx << std::endl
      << "#\t          cellsy: " << p.cellsy << std::endl
      << "#\t          cellsz: " << p.cellsz << std::endl;
  }
  f << "#" << std::endl;
  if (bySat) {
    f << "# SATNUM: " << volume.size() << std::endl;
    for (size_t i=0;i<volume.size();++i)
      f << "#\t SATNUM " << i+1 << ": " << volume[i]*100 << "%" << std::endl;
  } else {
    f <<"# Materials: " << volume.size() << std::endl;
    for (size_t i=0;i<volume.size();++i)
      f << "#\t Material" << i+1 << ": " << volume[i]*100 << "%" << std::endl;
  }
  if (upscaledRho > 0) {
    f << "#" << std::endl << "# Upscaled density: " << upscaledRho << std::endl;
    f << "#" << std::endl << "# Wave speeds: " << speeds << std::endl;
  }

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
      cellsize[0] = p.max[0] > -1?p.max[0]/cells[0]:1.0;
      cellsize[1] = p.max[1] > -1?p.max[1]/cells[1]:1.0;
      cellsize[2] = p.max[2] > -1?p.max[2]/cells[2]:1.0;
      grid.createCartesian(cells,cellsize);
    } else {
        Opm::ParseContext parseContext;
        Opm::Parser parser;
        auto deck = parser.parseFile(p.file , parseContext);
        Opm::EclipseGrid inputGrid(deck);
        grid.processEclipseFormat(inputGrid, false);
    }
    ElasticityUpscale<GridType, AMG> upscale(grid, p.ctol, p.Emin, p.file,
                                             p.rocklist, p.verbose);

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

    if (p.inspect == "mesh")
      return 0;

    if (p.linsolver.pre == UNDETERMINED) {
      double hx = (p.max[0]-p.min[0])/grid.logicalCartesianSize()[0];
      double hy = (p.max[1]-p.min[1])/grid.logicalCartesianSize()[1];
      double hz = (p.max[2]-p.min[2])/grid.logicalCartesianSize()[2];
      double aspect = sqrt(hx*hy)/hz;
      std::cout << "Estimated cell aspect ratio: " << aspect;
      if (aspect > 80) {
        p.linsolver.pre = TWOLEVEL;
        std::cout << " => Using two level preconditioner" << std::endl;
      } else {
        p.linsolver.pre = Opm::Elasticity::AMG;
        std::cout << " => Using AMG preconditioner" << std::endl;
      }
    }

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
    Opm::Elasticity::Vector field[6];
    std::cout << "assembling elasticity operator..." << "\n";
    upscale.assemble(-1,true);
    std::cout << "setting up linear solver..." << std::endl;
    upscale.setupSolvers(p.linsolver);

    if (p.inspect == "load")
      Dune::storeMatrixMarket(upscale.A.getOperator(), "A.mtx");

    // the uzawa solver cannot run multithreaded
#ifdef HAVE_OPENMP
    if (p.linsolver.uzawa) {
      std::cout << "WARNING: disabling multi-threaded solves due to uzawa" << std::endl;
      omp_set_num_threads(1);
    }
#endif

#pragma omp parallel for schedule(static)
    for (int i=0;i<6;++i) {
      std::cout << "processing case " << i+1 << "..." << std::endl;
      if (p.inspect == "results") {
        char temp[1024];
        sprintf(temp, p.resultfilename.c_str(), "x", i+1);
        Dune::loadMatrixMarket(upscale.u[i], temp);
      } else {
        std::cout << "\tassembling load vector..." << std::endl;
        upscale.assemble(i,false);
        if (p.inspect == "load") {
          char temp[1024];
          sprintf(temp, p.resultfilename.c_str(), "b", i+1);
          Dune::storeMatrixMarket(upscale.b[i], temp);
        }
        std::cout << "\tsolving..." << std::endl;
        upscale.solve(i);
        if (p.inspect == "load") {
          char temp[1024];
          sprintf(temp, p.resultfilename.c_str(), "x", i+1);
          Dune::storeMatrixMarket(upscale.u[i], temp);
        }
      }
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

    if (!p.vtufile.empty()) {
      Dune::VTKWriter<typename GridType::LeafGridView> vtkwriter(grid.leafGridView());

      for (int i=0;i<6;++i) {
        std::stringstream str;
        str << "sol " << i+1;
        vtkwriter.addVertexData(field[i], str.str().c_str(), dim);
      }
      vtkwriter.write(p.vtufile);
    }

    Dune::FieldVector<double,3> speeds;
    if (upscale.upscaledRho > -1) {
      speeds = Opm::Elasticity::waveSpeeds(C, p.dip, p.azimuth, 1.0);
      std::cout << "Wave speeds: " << speeds << std::endl;
    }

    // voigt notation
    for (int j=0;j<6;++j)
      std::swap(C[3][j],C[5][j]);
    for (int j=0;j<6;++j)
      std::swap(C[j][3],C[j][5]);
    std::cout << "---------" << std::endl;
    std::cout << C << std::endl;
    if (!p.output.empty())
      writeOutput(p, watch, grid.size(0), upscale.volumeFractions,
                  upscale.bySat, C, upscale.upscaledRho, speeds);

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

  template<template<class Smoother> class AMG>
int runAMG(Params& p)
{
  if (p.linsolver.smoother == SMOOTH_SCHWARZ)
    return run< Dune::CpGrid, AMG<SchwarzSmoother> >(p);
  else if (p.linsolver.smoother == SMOOTH_JACOBI)
    return run< Dune::CpGrid, AMG<JACSmoother> >(p);
  else if (p.linsolver.smoother == SMOOTH_ILU)
    return run< Dune::CpGrid, AMG<ILUSmoother> >(p);
  else
    return run<Dune::CpGrid, AMG<SSORSmoother> >(p);
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

    if (p.linsolver.pre == FASTAMG)
      return run<Dune::CpGrid, FastAMG>(p);
    else if (p.linsolver.pre == SCHWARZ)
      return run<Dune::CpGrid, Schwarz>(p);
    else if (p.linsolver.pre == TWOLEVEL)
      return runAMG<AMG2Level>(p);
    else
      return runAMG<AMG1>(p);
  } catch (const std::exception &e) {
    throw e;
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
