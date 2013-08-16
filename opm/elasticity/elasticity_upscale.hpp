//==============================================================================
//!
//! \file elasticity_upscale.hpp
//!
//! \date Nov 9 2011
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Elasticity upscale class
//!
//==============================================================================
#ifndef ELASTICITY_UPSCALE_HPP_
#define ELASTICITY_UPSCALE_HPP_

#include <dune/common/fmatrix.hh>
#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/istl/ilu.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/grid/CpGrid.hpp>
#include <opm/elasticity/shapefunctions.hpp>

#include <opm/elasticity/asmhandler.hpp>
#include <opm/elasticity/boundarygrid.hh>
#include <opm/elasticity/elasticity.hpp>
#include <opm/elasticity/logutils.hpp>
#include <opm/elasticity/materials.hh>
#include <opm/elasticity/mpc.hh>
#include <opm/elasticity/mortar_schur.hpp>
#include <opm/elasticity/mortar_utils.hpp>
#include <opm/elasticity/mortar_evaluator.hpp>
#include <opm/elasticity/mortar_schur_precond.hpp>
#include <opm/elasticity/uzawa_solver.hpp>

#include <dune/istl/superlu.hh>
#include <dune/istl/paamg/amg.hh>
#include <dune/istl/paamg/fastamg.hh>
#include <dune/istl/overlappingschwarz.hh>

#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>

namespace Opm {
namespace Elasticity {


//! \brief An enumeration of available linear solver classes
enum Solver {
  DIRECT,
  ITERATIVE 
};

enum Preconditioner {
  AMG,
  FASTAMG,
  SCHWARZ
};

//! \brief An enumeration of the available preconditioners for multiplier block
enum MultiplierPreconditioner {
  SCHUR,
  SIMPLE,
};

//! \brief Smoother used in the AMG
enum Smoother {
  SMOOTH_SSOR    = 0,
  SMOOTH_SCHWARZ = 1,
  SMOOTH_JACOBI  = 2,
  SMOOTH_ILU     = 4
};

struct LinSolParams {
  //! \brief The linear solver to employ
  Opm::Elasticity::Solver type;

  //! \brief Number of iterations in GMRES before restart
  int restart;

  //! \brief Max number of iterations
  int maxit;

  //! \brief The tolerance for the iterative linear solver
  double tol;

  //! \brief Use MINRES instead of GMRES (and thus symmetric preconditioning)
  bool symmetric;

  //! \brief Use a Uzawa approach
  bool uzawa;

  //! \brief The number of pre/post steps in the AMG
  int steps[2];

  //! \brief Coarsening target in the AMG
  int coarsen_target;

  //! \brief Number of cells in z to collapse in each cell
  int zcells;

  //! \brief Give a report at end of solution phase
  bool report;

  //! \brief Smoother type used in the AMG
  Smoother smoother;

  //! \brief Preconditioner for elasticity block
  Opm::Elasticity::Preconditioner pre;

  //! \brief Preconditioner for mortar block
  Opm::Elasticity::MultiplierPreconditioner mortarpre;

  //! \brief Parse command line parameters
  //! \param[in] param The parameter group to parse
  void parse(Opm::parameter::ParameterGroup& param)
  {
    std::string solver = param.getDefault<std::string>("linsolver_type","iterative");
    if (solver == "iterative")
      type = Opm::Elasticity::ITERATIVE;
    else
      type = Opm::Elasticity::DIRECT;
    restart = param.getDefault<int>("linsolver_restart", 1000);
    tol    = param.getDefault<double>("ltol",1.e-8);
    maxit   = param.getDefault<int>("linsolver_maxit", 10000);
    steps[0] = param.getDefault<int>("linsolver_presteps", 2);
    steps[1] = param.getDefault<int>("linsolver_poststeps", 2);
    coarsen_target = param.getDefault<int>("linsolver_coarsen", 5000);
    symmetric = param.getDefault<bool>("linsolver_symmetric", true);
    report = param.getDefault<bool>("linsolver_report", false);
    solver = param.getDefault<std::string>("linsolver_pre","amg");

    if (solver == "schwarz")
      pre = Opm::Elasticity::SCHWARZ;
    else if (solver == "fastamg")
      pre = Opm::Elasticity::FASTAMG;
    else
      pre = Opm::Elasticity::AMG;

    solver = param.getDefault<std::string>("linsolver_mortarpre","schur");
    if (solver == "schur")
      mortarpre = Opm::Elasticity::SCHUR;
    else
      mortarpre = Opm::Elasticity::SIMPLE;

    uzawa = param.getDefault<bool>("linsolver_uzawa", false);
    zcells = param.getDefault<int>("linsolver_zcells", 2);

    solver = param.getDefault<std::string>("linsolver_smoother","ssor");
    if (solver == "schwarz")
      smoother = SMOOTH_SCHWARZ;
    else if (solver == "ilu")
      smoother = SMOOTH_ILU;
    else if (solver == "jacobi")
      smoother = SMOOTH_JACOBI;
    else {
      if (solver != "ssor")
        std::cerr << "WARNING: Invalid smoother specified, falling back to SSOR" << std::endl;
      smoother = SMOOTH_SSOR;
    }

    if (symmetric)
      steps[1] = steps[0];
  }
};

//! \brief SSOR AMG smoother
typedef Dune::SeqSSOR<Matrix, Vector, Vector> SSORSmoother;
//! \brief GJ AMG smoother
typedef Dune::SeqJac<Matrix, Vector, Vector> JACSmoother;
//! \brief ILU0 AMG smoother
typedef Dune::SeqILU0<Matrix, Vector, Vector> ILUSmoother;
//! \brief Schwarz + ILU0 AMG smoother
typedef Dune::SeqOverlappingSchwarz<Matrix,Vector,
                              Dune::AdditiveSchwarzMode> SchwarzSmoother;

//! \brief The coupling metric used in the AMG
typedef Dune::Amg::FirstDiagonal CouplingMetric;

//! \brief The coupling criterion used in the AMG
typedef Dune::Amg::SymmetricCriterion<Matrix, CouplingMetric> CritBase;

//! \brief The coarsening criterion used in the AMG
typedef Dune::Amg::CoarsenCriterion<CritBase> Criterion;

//! \brief A linear operator
typedef Dune::MatrixAdapter<Matrix,Vector,Vector> Operator;

//! \brief A FastAMG for an elasticity operator
typedef Dune::Amg::FastAMG<Operator, Vector> FastAMG;

//! \brief Overlapping Schwarz preconditioner
typedef Dune::SeqOverlappingSchwarz<Matrix, Vector,
                                    Dune::AdditiveSchwarzMode,
                     Dune::SuperLU<Matrix> > Schwarz;

//! \brief Setup preconditioner
//! \param[in] pre The number of pre-smoothing steps
//! \param[in] post The number of post-smoothing steps
//! \param[in] target The coarsening target
//! \param[in] zcells The wanted number of cells to collapse in z per level
//! \param[in] op The linear operator
//! \param[in] gv The cornerpoint grid
//! \param[out] thread Whether or not to clone for threads
  template<class EAMG>
EAMG* setupPC(int pre, int post, int target, int zcells,
              std::shared_ptr<Operator>& op, const Dune::CpGrid& gv,
              ASMHandler<Dune::CpGrid>& A, bool& copy);

//! \brief The main driver class
  template<class GridType, class EAMG>
class ElasticityUpscale
{
  public:
    //! \brief Dimension of our grid
    static const int dim = GridType::dimension;

    //! \brief A basic number
    typedef typename GridType::LeafGridView::ctype ctype;

    //! \brief A vectorial node value
    typedef Dune::FieldVector<double,dim> NodeValue;

    //! \brief A global coordinate
    typedef typename GridType::LeafGridView::template Codim<1>::Geometry::GlobalCoordinate GlobalCoordinate;

    //! \brief A set of indices
    typedef typename GridType::LeafGridView::IndexSet LeafIndexSet;

    //! \brief An iterator over grid cells
    typedef typename GridType::LeafGridView::template Codim<0>::Iterator LeafIterator;

    //! \brief The linear operator
    ASMHandler<GridType> A;

    //! \brief The solution vectors
    Vector u[6];
    //! \brief The load vectors
    Vector b[6];

    //! \brief Vector holding the volume fractions for materials
    std::vector<double> volumeFractions;

    //! \brief Main constructor
    //! \param[in] gv_ The grid to operate on
    //! \param[in] tol_ The tolerance to use when deciding whether or not a coordinate falls on a plane/line/point. \sa tol
    //! \param[in] Escale_ A scale value for E-moduluses to avoid numerical issues
    //! \param[in] file The eclipse grid file
    //! \param[in] rocklist If not blank, file is a rocklist
    //! \param[in] verbose If true, give verbose output
    ElasticityUpscale(const GridType& gv_, ctype tol_, ctype Escale_, 
                      const std::string& file, const std::string& rocklist,
                      bool verbose_)
      :  A(gv_), gv(gv_), tol(tol_), Escale(Escale_), E(gv_), verbose(verbose_)
    {
      if (rocklist.empty())
        loadMaterialsFromGrid(file);
      else
        loadMaterialsFromRocklist(file,rocklist);
    }

    //! \brief Find boundary coordinates
    //! \param[out] min The miminum coordinates of the grid
    //! \param[out] max The maximum coordinates of the grid
    void findBoundaries(double* min, double* max);

    //! \brief Add a MPC equation
    //! \param[in] dir The direction of the MPC
    //! \param[in] slave The slave node index
    //! \param[in] m The vertices on the master grid
    void addMPC(Direction dir, int slave, 
                const BoundaryGrid::Vertex& m);

    //! \brief Establish periodic boundaries using the MPC approach
    //! \param[in] min The minimum coordinates of the grid
    //! \param[in] max The maximum coordinates of the grid
    void periodicBCs(const double* min, const double* max);

    //! \brief Establish periodic boundaries using the mortar approach
    //! \param[in] min The minimum coordinates of the grid
    //! \param[in] max The maximum coordinates of the grid
    //! \param[in] n1 The number of elements on the lambda grid in the X direction
    //! \param[in] n2 The number of elements on the lambda grid in the Y direction
    //! \param[in] p1 The order of multipliers in the X direction
    //! \param[in] p2 The order of multipliers in the Y direction
    void periodicBCsMortar(const double* min,
                           const double* max, int n1, int n2,
                           int p1, int p2);

    //! \brief Fix corner nodes
    //! \param[in] min The minimum coordinates on the grid
    //! \param[in] max The maximum coordinates on the grid
    void fixCorners(const double* min, const double* max);

    //! \brief Assemble (optionally) stiffness matrix A and load vector
    //! \param[in] loadcase The strain load case. Set to -1 to skip
    //! \param[in] matrix Whether or not to assemble the matrix
    void assemble(int loadcase, bool matrix);

    //! \brief Calculate the average stress vector for the given loadcase
    //! \param[out] sigma The stress vector
    //! \param[in] u The displacement vector
    //! \param[in] loadcase The strain load case considered
      template <int comp>
    void averageStress(Dune::FieldVector<ctype,comp>& sigma,
                       const Vector& u, int loadcase);

    //! \brief Solve Au = b for u
    //! \param[in] loadcase The load case to solve
    void solve(int loadcase);

    //! \param[in] params The linear solver parameters
    void setupSolvers(const LinSolParams& params);

  private:
    //! \brief An iterator over grid vertices
    typedef typename GridType::LeafGridView::template Codim<dim>::Iterator LeafVertexIterator;

    //! \brief A reference to our grid
    const GridType& gv;

    //! \brief Tolerance used to decide whether or not a coordinate falls on a plane/line/point.
    ctype tol;

    //! \brief Minimum E-modulus (scaling factor)
    ctype Escale;

    //! \brief Minimum real E for materials
    ctype Emin;

    //! \brief Check if the given coordinate falls on a given plane
    //! \param[in] plane The plane of interest
    //! \param[in] coord The coordinates to check
    //! \param[in] value The constant coordinate describing the plane
    bool isOnPlane(Direction plane, GlobalCoordinate coord, ctype value);

    //! \brief Check if the given coordinate falls on a given line 
    //! \param[in] dir The line direction 
    //! \param[in] coord The coordinates to check
    //! \param[in] x The first coordinate of the line
    //! \param[in] y The second coordinate of the line
    bool isOnLine(Direction dir, GlobalCoordinate coord, ctype x, ctype y);

    //! \brief Check if the given coordinate is the given point
    //! \param[in] coord The coordinates to check
    //! \param[in] point The point
    bool isOnPoint(GlobalCoordinate coord, GlobalCoordinate point);

    //! \brief Vector holding material parameters for each active grid cell
    std::vector< std::shared_ptr<Material> > materials;

    //! \brief Extract the vertices on a given face
    //! \param[in] dir The direction of the face normal
    //! \param[in] coord The coordinate of the face plane
    //! \returns A vector holding the matching vertices
    std::vector<BoundaryGrid::Vertex> extractFace(Direction dir, ctype coord);

    //! \brief An enumeration used to indicate which side to extract from a cube
    enum SIDE {
      LEFT,
      RIGHT
    };

    //! \brief Extract a quad grid over a given face
    //! \param[in] dir The direction of the face normal
    //! \param[in] coord the coordinate of the face plance
    //! \param[in] side Extract left or right side
    //! \param[in] dc If true, order vertices in dune convention
    //! \returns A quad grid spanning the face
    BoundaryGrid extractMasterFace(Direction dir, ctype coord,
                                   SIDE side=LEFT, bool dc=false);

    //! \brief Find and establish master/slave grids (MPC)
    //! \param[in] min The minimum coordinates of the grid
    //! \param[in] max The maximum coordinates of the grid
    void determineSideFaces(const double* min, const double* max);

    //! \brief Establish periodic boundary conditions on a given plane using MPC couplings
    //! \param[in] plane The direction of the plane normal
    //! \param[in] dir The coordinate direction to enforce periodicity in
    //! \param[in] slave The slave point grid
    //! \param[in] master The master quad grid
    void periodicPlane(Direction plane, Direction dir,
                       const std::vector<BoundaryGrid::Vertex>& slave,
                       const BoundaryGrid& master);

    //! \brief Fix the DOFs in a given point on the grid
    //! \param[in] dir The coordinate direction to fix in
    //! \param[in] coord The coordinates of the node to fix
    //! \param[in] value The values to fix the given DOFs to
    void fixPoint(Direction dir, GlobalCoordinate coord,
                  const NodeValue& value = NodeValue(0));

    //! \brief Fix the DOFs in a given line on the grid
    //! \param[in] dir The coordinate direction to fix in
    //! \param[in] x The first coordinate of the line
    //! \param[in] y The second coordinate of the line
    //! \param[in] value The values to fix the given DOFs to
    void fixLine(Direction dir, ctype x, ctype y,
                 const NodeValue& value = NodeValue(0));

    //! \brief A point grid
    typedef std::map<int,BoundaryGrid::Vertex> SlaveGrid;

    //! \brief Add a block to the B matrix associated with mortar couplings
    //! \param[in] b1 The primal boundary to match
    //! \param[in] interface The interface/multiplier grid
    //! \param[in] dir The normal direction on the boundary (0/1)
    //! \param[in] n1 The multipler order in the first direction
    //! \param[in] n2 The multipler order in the second direction
    //! \param[in] colofs The column offset (multiplier number)
    //! \returns Number of multipliers DOFs added
    int addBBlockMortar(const BoundaryGrid& b1,
                        const BoundaryGrid& interface, int dir,
                        int n1, int n2, int colofs);

    //! \brief Assemble part of the B block associated with mortar couplings
    //! \param[in] b1 The primal boundary to match
    //! \param[in] interface The interface/multiplier grid
    //! \param[in] dir The normal direction on the boundary (0/1)
    //! \param[in] n1 The multipler order in the first direction
    //! \param[in] n2 The multipler order in the second direction
    //! \param[in] colofs The column offset (first multiplier unknown)
    //! \param[in] alpha Scaling for matrix elements
    void assembleBBlockMortar(const BoundaryGrid& b1,
                              const BoundaryGrid& interface, int dir,
                              int n1, int n2, int colofs, double alpha=1.0);

    //! \brief This function loads and maps materials to active grid cells
    //! \param[in] file The eclipse grid to read materials from
    void loadMaterialsFromGrid(const std::string& file);

    //! \brief This function loads and maps materials to active grid cells
    //! \param[in] file The grid file to read SATNUM from
    //! \param[in] rocklist The rock list to read materials from
    void loadMaterialsFromRocklist(const std::string& file,
                                   const std::string& rocklist);

    //! \brief Master grids
    std::vector<BoundaryGrid> master;

    //! \brief Slave point grids
    std::vector< std::vector<BoundaryGrid::Vertex> > slave;

    //! \brief Lagrangian multiplier block
    AdjacencyPattern Bpatt;
    Matrix B;

    Matrix P; //!< Preconditioner for multiplier block

    //! \brief Linear solver
    typedef std::shared_ptr<Dune::InverseOperator<Vector, Vector> > SolverPtr;
    SolverPtr solver;

    //! \brief Matrix adaptor for the elasticity block
    std::shared_ptr<Operator> op;

    //! \brief Preconditioner for multiplier block
    typedef MortarBlockEvaluator<Dune::Preconditioner<Vector, Vector> > SchurPreconditioner;

    //! \brief Evaluator for multiplier block
    typedef MortarBlockEvaluator<Dune::InverseOperator<Vector, Vector> > SchurEvaluator;

    //! \brief Outer evaluator, used with uzawa
    std::shared_ptr<SchurEvaluator> op2;

    //! \brief The preconditioner for the elasticity operator
    std::shared_ptr<EAMG> upre;

    //! \brief An LU solve as a preconditioner
    typedef Dune::InverseOperator2Preconditioner<Dune::SuperLU<Matrix>,
                                        Dune::SolverCategory::sequential> SeqLU;
    //! \brief The preconditioner for the multiplier block (used with uzawa)
    std::shared_ptr<SeqLU> lpre;
    std::shared_ptr< Dune::SuperLU<Matrix> > lprep;

    //! \brief Preconditioner for the Mortar system
    typedef std::shared_ptr< MortarSchurPre<EAMG> > MortarAmgPtr;
    MortarAmgPtr tmpre;

    //! \brief Evaluator for the Mortar system
    std::shared_ptr<MortarEvaluator> meval;

    //! \brief Elasticity helper class
    Elasticity<GridType> E;

    //! \brief Verbose output
    bool verbose;
};

}} // namespace Opm, Elasticity

#include "elasticity_upscale_impl.hpp"

#endif
