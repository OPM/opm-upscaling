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
#pragma once

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <dune/common/fmatrix.hh>
#include <opm/core/eclipse/EclipseGridParser.hpp>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/istl/ilu.hh>
#include <dune/istl/solvers.hh>
#if HAVE_SUPERLU
#include <dune/istl/superlu.hh>
#endif
#include <dune/istl/preconditioners.hh>
#include <dune/istl/overlappingschwarz.hh>
#include <dune/grid/CpGrid.hpp>
#include <dune/elasticity/shapefunctions.hpp>

#include <dune/elasticity/asmhandler.hpp>
#include <dune/elasticity/boundarygrid.hh>
#include <dune/elasticity/elasticity.hpp>
#include <dune/elasticity/materials.hh>
#include <dune/elasticity/mpc.hh>

namespace Opm {
namespace Elasticity {

//! \brief An enumeration of available linear solvers
enum Solver {
  SLU,
   CG 
};

typedef Dune::SeqOverlappingSchwarz<Matrix,Vector,Dune::SymmetricMultiplicativeSchwarzMode> OverlappingSchwarz;

//! \brief The main driver class
  template<class GridType>
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
    //! \param[in] rocklist If true, file is a rocklist
    ElasticityUpscale(const GridType& gv_, ctype tol_, ctype Escale_, 
                      const std::string& file, const std::string& rocklist)
      : gv(gv_), tol(tol_), A(gv_), E(gv_), Escale(Escale_)
    {
      if (rocklist.empty())
        loadMaterialsFromGrid(file);
      else
        loadMaterialsFromRocklist(file,rocklist);
#if HAVE_SUPERLU
      slu = 0;
#endif
      cgsolver = 0;
      op = 0;
      ilu = 0;
      ovl = 0;
    }

    //! \brief The destructor
    ~ElasticityUpscale()
    {
      // sort the pointers so unique can do its job
      std::sort(materials.begin(),materials.end());
      // this reorders the vector so we only get one entry per pointer
      std::vector<Material*>::iterator itend = std::unique(materials.begin(),materials.end());
      // now delete the pointers
      for (std::vector<Material*>::iterator it  = materials.begin();
                                            it != itend; ++it)
        delete *it;

#if HAVE_SUPERLU
      delete slu;
#endif
      delete cgsolver;
      delete ilu;
      delete op;
      delete ovl;
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

    //! \brief Establish periodic boundaries using the LLM approach
    //! \param[in] min The minimum coordinates of the grid
    //! \param[in] max The maximum coordinates of the grid
    //! \param[in] n1 The number of elements on the interface grid in the X direction
    //! \param[in] n2 The number of elements on the interface grid in the Y direction
    void periodicBCsLLM(const double* min,
                        const double* max,
                        int n1, int n2);

    //! \brief Establish periodic boundaries using the mortar approach
    //! \param[in] min The minimum coordinates of the grid
    //! \param[in] max The maximum coordinates of the grid
    //! \param[in] n1 The number of elements on the lambda grid in the X direction
    //! \param[in] n2 The number of elements on the lambda grid in the Y direction
    void periodicBCsMortar(const double* min,
                           const double* max, int n1, int n2);

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
    //! \param[in] solver The linear equation solver to employ
    //! \param[in] tol The tolerance for iterative solvers
    void solve(Solver solver, double tol, int loadcase);

    void setupSolvers(Solver solver);
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
    std::vector<Material*> materials;

    //! \brief Fix corner nodes
    //! \param[in] min The minimum coordinates on the grid
    //! \param[in] max The maximum coordinates on the grid
    void fixCorners(const double* min, const double* max);

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

    //! \brief Find the B matrix associated with a particular slave grid (LLM)
    //! \param[in] slave The slave grid to extract
    //! \returns The B matrix associated with the given slave grid
    Matrix findBMatrixLLM(const SlaveGrid& slave);

    //! \brief Find the L matrix associated with a particular slave grid (LLM)
    //! \param[in] slave The slave grid we want to couple
    //! \param[in] master The interface grid we want to couple to
    //! \returns The L matrix describing the coupling between slave and master 
    Matrix findLMatrixLLM(const SlaveGrid& slave, const BoundaryGrid& master);

    //! \brief Find the L matrix associated with mortar couplings
    //! \param[in] b1 The primal boundary to match
    //! \param[in] interface The interface/multiplier grid
    //! \param[in] dir The normal direction on the boundary (0/1)
    Matrix findLMatrixMortar(const BoundaryGrid& b1,
                             const BoundaryGrid& interface, int dir);

    //! \brief This function loads and maps materials to active grid cells
    //! \param[in] file The eclipse grid to read materials from
    void loadMaterialsFromGrid(const std::string& file);

    //! \brief This function loads and maps materials to active grid cells
    //! \param[in] file The grid file to read SATNUM from
    //! \param[in] rocklist The rock list to read materials from
    void loadMaterialsFromRocklist(const std::string& file,
                                   const std::string& rocklist);

    //! \brief Setup an overlapping schwarz preconditioner with one 
    //!        subdomain per pillar.
    void setupPreconditioner();

    //! \brief Master grids
    std::vector<BoundaryGrid> master;

    //! \brief Slave point grids
    std::vector< std::vector<BoundaryGrid::Vertex> > slave;

    //! \brief Vector of matrices holding boolean coupling matrices (LLM)
    std::vector<Matrix> B;
    
    //! \brief Vector of matrices holding Lagrangian multipliers
    std::vector<Matrix> L;

#if HAVE_SUPERLU
    //! \brief SuperLU solver
    Dune::SuperLU<Matrix>* slu;
#endif
    //! \brief CG solver
    Dune::CGSolver<Vector>* cgsolver;
    //! \brief The linear operator
    Dune::MatrixAdapter<Matrix,Vector,Vector>* op;

    //! \brief ILU preconditioner
    Dune::SeqILU0<Matrix,Vector,Vector>* ilu;
    //! \brief Overlapping Schwarz preconditioner
    OverlappingSchwarz* ovl;

    //! \brief Elasticity helper class
    Elasticity<GridType> E;
};

#include "elasticity_upscale_impl.hpp"

}
}
