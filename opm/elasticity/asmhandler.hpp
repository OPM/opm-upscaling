//==============================================================================
//!
//! \file asmhandler.hpp
//!
//! \date Nov 9 2011
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Class handling finite element assembly
//!
//==============================================================================
#ifndef ASMHANDLER_HPP_
#define ASMHANDLER_HPP_

#include <opm/common/utility/platform_dependent/disable_warnings.h>

#include <dune/geometry/referenceelements.hh>
#include <dune/common/fmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/common/fvector.hh>

#include <opm/common/utility/platform_dependent/reenable_warnings.h>


#include <opm/elasticity/logutils.hpp>
#include <opm/elasticity/mpc.hh>
#include <opm/elasticity/matrixops.hpp>

namespace Opm {
namespace Elasticity {


//!\brief Class handling finite element assembly
  template<class GridType>
class ASMHandler {
  public:
    //! \brief The dimension of the grid
    static const int dim = GridType::dimension;

    //! \brief A set of indices
    typedef typename GridType::LeafGridView::IndexSet LeafIndexSet;

    //! \brief An iterator over grid cells
    typedef typename GridType::LeafGridView::template Codim<0>::Iterator LeafIterator;

    //! \brief The default constructor
    //! \param[in] gv_ The grid the operator is assembled over
    ASMHandler(const GridType& gv_) : gv(gv_), maxeqn(0)
    {
    }

    //! \brief Destructor
   ~ASMHandler()
   {
     for (MPCMap::iterator it=mpcs.begin(); it != mpcs.end();++it)
       delete it->second;
   }

    //! \brief A vectorial node value
    typedef Dune::FieldVector<double,dim> NodeValue;

    //! \brief Get the number of equations in the system
    //! \returns The number of equations in the system
    size_t getEqns() const
    {
      return maxeqn;
    }

    //! \brief Get the equation number for a given dof on a given node
    //! \param[in] node The node number
    //! \param[in] dof The DOF
    //! \returns The equation number if active, -1 otherwise
    int getEquationForDof(int node, int dof)
    {
      return meqn[node*dim+dof];
    }

    //! \brief Obtain a reference to the linear operator
    //! \returns Reference to linear operator
    Matrix& getOperator()
    {
      return A;
    }

    //! \brief Obtain a reference to the load vector
    //! \returns Reference to load vector
    Vector& getLoadVector()
    {
      return b;
    }

    //! \brief This function needs to be called before starting
    //!        the element assembly. 
    void initForAssembly();

    //! \brief Add an element matrix/vector to the system
    //! \param[in] K Pointer to the element matrix. Can be NULL
    //! \param[in] S Pointer to the element load vector. Can be NULL
    //! \param[in] cell An iterator pointing to the cell we're assembling for
    //! \param[in] b Vector to add contributions to. If not given, 
    //!              we use the internal vector
      template<int esize>
    void addElement(const Dune::FieldMatrix<double,esize,esize>* K,
                    const Dune::FieldVector<double,esize>* S,
                    const LeafIterator& cell,
                    Vector* b=NULL);
    void addMatElement(int i,int j,double val){ A[i][j] += val;}  
    //! \brief Extract values corresponding to cell
    //! \param[in] u The global load vector
    //! \param[in] it An iterator to the cell we want to extract values for
    //! \param[out] v Vector holding the values requested
      template<int comp>
    void extractValues(Dune::FieldVector<double,comp>& v, 
                       const Vector& u, const LeafIterator& it);

    //! \brief Expand a system vector to a solution vector
    void expandSolution(Vector& result, const Vector& u);

    //! \brief Add a MPC
    //! \param[in] mpc Pointer to the MPC to add.
    //! \note This class handles destruction
    void addMPC(MPC* mpc);

    //! \brief Look for and return a MPC for a specified node+dof
    //! \param[in] node The requested node
    //! \param[in] dof The requested DOF at given node
    //! \returns The MPC for the node/dof if found, else NULL
    MPC* getMPC(int node, int dof);

    //! \brief Update/add a fixed node
    //! \param[in] node The node number
    //! \param[in] entry The fixed values
    void updateFixedNode(int node, 
                         const std::pair<Direction,NodeValue>& entry);

    //! \brief Check if a node is marked as fixed (in any direction)
    //! \param[in] node The node to query for
    bool isFixed(int node)
    {
      return (fixedNodes.find(node) != fixedNodes.end());
    }

    //! \brief Print the current operator
    void printOperator() const;

    //! \brief Print the current load vector
    void printLoadVector() const;

    //! \brief Access current adjacency pattern
    //! \details Can be used to add extra entries, such as other blocks
    AdjacencyPattern& getAdjacencyPattern()
    {
      return adjacencyPattern;
    }
  protected:
    //! \brief Resolve chained MPCs
    void resolveMPCChains()
    {
      for (MPCMap::iterator it  = mpcs.begin();
                            it != mpcs.end();++it)
        resolveMPCChain(it->second);
    }

    //! \brief Internal function. Handles a single MPC
    //! \param[in] mpc Pointer to the MPC to resolve
    void resolveMPCChain(MPC* mpc);

    //! \brief Internal function. Generate meqn for registered MPC/fixed nodes
    void preprocess();

    //! \brief Internal function. Generate adjacency pattern for a given node
    //! \param[in] it Iterator pointing to the cell in of the node
    //! \param[in] vertexsize Number of vertices in the cell
    //! \param[in] row The equation number/row in matrix
    void nodeAdjacency(const LeafIterator& it, int vertexsize, int row);

    //! \brief Internal function. Calculate adjacency pattern 
    void determineAdjacencyPattern();

    //! \brief Internal function. Assemble entries for a single DOF
    //! \param[in] row The row in the global matrix
    //! \param[in] erow The row in the element matrix
    //! \param[in] K Pointer to the element matrix. Can be NULL
    //! \param[in] S Pointer to the element load vector. Can be NULL
    //! \param[in] set The index set
    //! \param[in] cell An iterator pointing to the cell we're assembling for
    //! \param[in] b Vector to add contributions to
    //! \param[in] scale Scale for elements. Used with MPC couplings
      template<int esize>
    void addDOF(int row, int erow,
                const Dune::FieldMatrix<double,esize,esize>* K,
                const Dune::FieldVector<double,esize>* S,
                const LeafIndexSet& set,
                const LeafIterator& cell,
                Vector* b,
                double scale=1.f);

    //! \brief The set of MPC
    MPCMap mpcs;

    //! \brief Vector of (interleaved) dof<->eqn mapping
    std::vector<int> meqn;

    //! \brief Fixed nodes
    typedef std::pair<Direction,NodeValue> fixEntry;

    //! \brief A mapping from dof to fix value info
    typedef std::map<int,fixEntry> fixMap;

    //! \brief Iterator over a fixmap
    typedef typename fixMap::iterator fixIt;

    //! \brief The map holding information about our fixed nodes
    fixMap fixedNodes;

    //! \brief Holds the adjacency pattern of the sparse matrix
    AdjacencyPattern adjacencyPattern;

    //! \brief The linear operator
    Matrix A;

    //! \brief The load vector
    Vector b;

    //! \brief A reference to the grid in use
    const GridType& gv;

    //! \brief The number of equations in the system
    size_t maxeqn;
  private:
    //! \brief No copying of this class
    ASMHandler(const ASMHandler&) {}
    //! \brief No copying of this class
    ASMHandler& operator=(const ASMHandler&) {}
};

}
}

#include "asmhandler_impl.hpp"

#endif
