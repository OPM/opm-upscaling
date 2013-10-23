//==============================================================================
//!
//! \file asmhandler_impl.hpp
//!
//! \date Nov 9 2011
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Class handling finite element assembly - template implementations
//!
//==============================================================================
#ifndef ASMHANDLER_IMPL_HPP_
#define ASMHANDLER_IMPL_HPP_

#include <dune/common/version.hh>
#include <iostream>

namespace Opm {
namespace Elasticity {

  template<class GridType>
void ASMHandler<GridType>::initForAssembly()
{
  resolveMPCChains();
  preprocess();
  determineAdjacencyPattern();

#if !DUNE_VERSION_NEWER(DUNE_ISTL, 2, 3)
  // workaround a bug in bcrs matrix
  A.setBuildMode(Matrix::random);
  A.endrowsizes();
#endif

  MatrixOps::fromAdjacency(A,adjacencyPattern,
                           adjacencyPattern.size(),adjacencyPattern.size());
  b.resize(adjacencyPattern.size());
  b = 0;
  adjacencyPattern.clear();

  // print some information
  std::cout << "\tNumber of nodes: " << gv.size(dim) << std::endl;
  std::cout << "\tNumber of elements: " << gv.size(0) << std::endl;
  std::cout << "\tNumber of constraints: " << mpcs.size() << std::endl;
  int fixedDofs=0;
  for (fixIt it = fixedNodes.begin(); it != fixedNodes.end(); ++it) {
    if (it->second.first & X)
      fixedDofs++;
    if (it->second.first & Y)
      fixedDofs++;
    if (it->second.first & Z)
      fixedDofs++;
  }
  std::cout << "\tNumber of fixed dofs: " << fixedDofs << std::endl;
}

  template<class GridType>
    template<int esize>
void ASMHandler<GridType>::addDOF(int row, int erow,
                              const Dune::FieldMatrix<double,esize,esize>* K,
                              const Dune::FieldVector<double,esize>* S,
                              const LeafIndexSet& set,
                              const LeafIterator& cell, 
                              Vector* b,
                              double scale)
{
  if (row == -1)
    return;
  if (K) {
    for (int j=0;j<esize/dim;++j) {
      int index2 = set.subIndex(*cell,j,dim);
      for (int l=0;l<dim;++l) {
        MPC* mpc = getMPC(index2,l);
        if (mpc) {
          for (size_t n=0;n<mpc->getNoMaster();++n) {
            int idx = meqn[mpc->getMaster(n).node*dim+mpc->getMaster(n).dof-1];
            if (idx != -1)
              A[row][idx] += scale*mpc->getMaster(n).coeff*(*K)[erow][j*dim+l];
          }
        } else if (meqn[index2*dim+l] != -1) {
          A[row][meqn[index2*dim+l]] += scale*(*K)[erow][j*dim+l];
        }
      }
    }
  }
  if (S && b)
    (*b)[row] += scale*(*S)[erow];
}

  template<class GridType>
    template<int esize>
void ASMHandler<GridType>::addElement(
                    const Dune::FieldMatrix<double,esize,esize>* K,
                    const Dune::FieldVector<double,esize>* S,
                    const LeafIterator& cell,
                    Vector* b2)
{
  if (!b2)
    b2 = &b;
  const LeafIndexSet& set = gv.leafView().indexSet();
  for (int i=0;i<esize/dim;++i) {
    int index1 = set.subIndex(*cell,i,dim);
    fixIt it = fixedNodes.find(index1);
    if (it != fixedNodes.end() && it->second.first == XYZ)
      continue;
    for (int k=0;k<dim;++k) {
      MPC* mpc = getMPC(index1,k);
      if (mpc) {
        for (size_t n=0;n<mpc->getNoMaster();++n) {
          int idx = meqn[mpc->getMaster(n).node*dim+mpc->getMaster(n).dof-1];
          addDOF(idx,i*dim+k,K,S,set,cell,b2,mpc->getMaster(n).coeff);
        }
      } else
        addDOF(meqn[index1*dim+k],i*dim+k,K,S,set,cell,b2);
    }
  }
}

  template<class GridType>
    template<int comp>
void ASMHandler<GridType>::extractValues(Dune::FieldVector<double,comp>& v,
                                         const Vector& u,
                                         const LeafIterator& it)
{
  v = 0;
  const LeafIndexSet& set = gv.leafView().indexSet();
  Dune::GeometryType gt = it->type();

#if DUNE_VERSION_NEWER(DUNE_GRID, 2, 3)
  const Dune::template ReferenceElement<double,dim> &ref =
                      Dune::ReferenceElements<double,dim>::general(gt);
#else
  const Dune::template GenericReferenceElement<double,dim> &ref =
                      Dune::GenericReferenceElements<double,dim>::general(gt);
#endif
  int vertexsize = ref.size(dim);
  int l=0;
  for (int i=0;i<vertexsize;++i) {
    int indexi = set.subIndex(*it,i,dim);
    fixIt it2 = fixedNodes.find(indexi);
    for (int n=0;n<dim;++n) {
      MPC* mpc = getMPC(indexi,n);
      if (it2 != fixedNodes.end() && it2->second.first & (1 << n))
        v[l++] = it2->second.second[n];
      else if (mpc) {
        for (size_t m=0;m<mpc->getNoMaster();++m) {
          int idx = meqn[mpc->getMaster(m).node*dim+mpc->getMaster(m).dof-1];
          if (idx != -1)
            v[l] += u[idx]*mpc->getMaster(m).coeff;
        }
        l++;
      } else
        v[l++] = u[meqn[indexi*dim+n]];
    }
  }
}

  template<class GridType>
void ASMHandler<GridType>::expandSolution(Vector& result, const Vector& u)
{
  int nodes = gv.size(dim);
  result.resize(nodes*dim);
  result = 0;
  int l=0;
  for (int i=0;i<nodes;++i) {
    fixIt it = fixedNodes.find(i);
    Direction dir;
    if (it == fixedNodes.end())
      dir = NONE;
    else
      dir = it->second.first;

    int flag=1;
    for (int j=0;j<dim;++j) {
      if (dir & flag)
        result[l] = it->second.second[j];
      else if (meqn[l] != -1)
        result[l] = u[meqn[l]];
      l++;
      flag *= 2;
    }
  }
  // second loop - handle MPC couplings
  l = 0;
  for (int i=0;i<nodes;++i) {
    for (int j=0;j<dim;++j) {
      MPC* mpc = getMPC(i,j);
      if (mpc) {
        for (size_t n=0;n<mpc->getNoMaster();++n) {
          int idx = mpc->getMaster(n).node*dim+mpc->getMaster(n).dof-1;
          if (meqn[idx] != -1)
            result[l] += u[meqn[idx]]*mpc->getMaster(n).coeff;
        }
      }
      ++l;
    }
  }
}

  template<class GridType>
void ASMHandler<GridType>::addMPC(MPC* mpc)
{
  if (!mpc)
    return;

  int slaveNode = mpc->getSlave().node*dim+mpc->getSlave().dof-1;
  fixIt it = fixedNodes.find(mpc->getSlave().node);
  int flag = 1 << (mpc->getSlave().dof-1);
  if (it == fixedNodes.end() || 
      !(it->second.first & flag)) {
    mpcs.insert(std::make_pair(slaveNode,mpc));
    return;
  }

  delete mpc;
}

  template<class GridType>
MPC* ASMHandler<GridType>::getMPC(int node, int dof)
{
  if (mpcs.find(node*dim+dof) != mpcs.end())
    return mpcs[node*dim+dof];

  return NULL;
}

  template<class GridType>
void ASMHandler<GridType>::updateFixedNode(int node,
                                  const std::pair<Direction,NodeValue>& entry)
{
  fixIt it = fixedNodes.find(node);
  // same type or new - update values/add and return
  if (it == fixedNodes.end() || it->second.first == entry.first) {
    fixedNodes[node] = entry;
    return;
  }
  int temp = it->second.first;
  temp |= entry.first;
  it->second.first = (Direction)temp;
  int flag = 1;
  for (int i=0;i<dim;++i) {
    if (entry.first & flag)
      it->second.second[i] = entry.second[i];
    flag *= 2;
  }
}

  template<class GridType>
void ASMHandler<GridType>::printOperator() const
{
  MatrixOps::print(A);
}

  template<class GridType>
void ASMHandler<GridType>::printLoadVector() const
{
  for (int i=0;i<b.size();++i) {
    double val = b[i];
    std::cout << (fabs(val)>1.e-12?val:0.f) << std::endl;
  }
}

  template<class GridType>
void ASMHandler<GridType>::resolveMPCChain(MPC* mpc)
{
  size_t nMaster = mpc->getNoMaster();
  if (nMaster == 0) return; // no masters, prescribed displacement only

  for (size_t i = 0; i < nMaster; i++)
  {
    // Check if the master node has a local coordinate system attached. If yes,
    // the slave DOF might be coupled to all (local) DOFs of the master node.
    const MPC::DOF& master = mpc->getMaster(i);
    Dune::FieldVector<double,dim> coeff;
    coeff = 0;
    coeff[master.dof-1] = master.coeff;

    int removeOld = 0;
    for (int d = 1; d <= dim; d++)
      if (fabs(coeff[d-1]) > 1.0e-8)
      {
        MPC* mpc2 = getMPC(mpc->getMaster(i).node,d-1);
        if (mpc2)
        {
          // We have a master DOF which is a slave in another MPC.
          // Invoke resolveMPCchain recursively to ensure that all master DOFs
          // of that equation are not slaves themselves.
          resolveMPCChain(mpc2);

          // Remove current master specification, unless it has been updated
          if (!removeOld) removeOld = 1;

          // Add constant offset from the other equation
          mpc->addOffset(mpc2->getSlave().coeff*coeff[d-1]);

          // Add masters from the other equations
          for (size_t j = 0; j < mpc2->getNoMaster(); j++)
            mpc->addMaster(mpc2->getMaster(j).node,
                           mpc2->getMaster(j).dof,
                           mpc2->getMaster(j).coeff*coeff[d-1]);
        }
        else
          // The master node is free, but has a local coordinate system
          if (d != mpc->getMaster(i).dof)
            // Add coupling to the other local DOFs of this master node.
            mpc->addMaster(mpc->getMaster(i).node,d,coeff[d-1]);
          else if (coeff[d-1] != mpc->getMaster(i).coeff)
          {
            // Update the coupling coefficient of this master DOF
            // due to the local-to-global transformation
            mpc->updateMaster(i,coeff[d-1]);
            removeOld = -1;
          }
      }
      else if (d == mpc->getMaster(i).dof && !removeOld)
        // The coefficient of the current master DOF is zero,
        removeOld = 1; // so remove it from the contraint equation

    if (removeOld == 1)
    {
      // Remove the old master DOF specification when it has been replaced
      mpc->removeMaster(i--);
      nMaster--; // we don't need to check the added masters
    }
  }
}

  template<class GridType>
void ASMHandler<GridType>::preprocess()
{
  int nodes = gv.size(dim);
  meqn.resize(nodes*dim);

  // iterate over nodes
  for (int indexi=0;indexi<nodes;++indexi) {
    fixIt it2 = fixedNodes.find(indexi);
    if (it2 == fixedNodes.end()) {
      for (int i=0;i<dim;++i) {
        MPC* mpc = getMPC(indexi,i);
        if (!mpc)
          meqn[indexi*dim+i] = maxeqn++;
        else
          meqn[indexi*dim+i] = -1;
      }
    } else {
      int flag=1;
      for (int i=0;i<dim;++i) {
        if (it2->second.first & flag)
          meqn[indexi*dim+i] = -1;
        else {
          MPC* mpc = getMPC(indexi,i);
          if (!mpc)
            meqn[indexi*dim+i] = maxeqn++;
          else
            meqn[indexi*dim+i] = -1;
        }
        flag *= 2;
      }
    }
  }
  std::cout << "\tnumber of equations: " << maxeqn << std::endl;
}

  template<class GridType>
void ASMHandler<GridType>::nodeAdjacency(const LeafIterator& it,
                                         int vertexsize, int row)
{
  if (row == -1)
    return;
  const LeafIndexSet& set = gv.leafView().indexSet();
  for (int j=0;j<vertexsize;++j) {
    int indexj = set.subIndex(*it,j,dim);
    for (int l=0;l<dim;++l) {
      MPC* mpc = getMPC(indexj,l);
      if (mpc) {
        for (size_t i=0;i<mpc->getNoMaster();++i) {
          int idx = meqn[mpc->getMaster(i).node*dim+
            mpc->getMaster(i).dof-1];
          if (idx != -1)
            adjacencyPattern[row].insert(idx);
        }
      } else if (meqn[indexj*dim+l] != -1)
        adjacencyPattern[row].insert(meqn[indexj*dim+l]);
    }
  }
}

  template<class GridType>
void ASMHandler<GridType>::determineAdjacencyPattern()
{
  adjacencyPattern.resize(maxeqn);
  std::cout << "\tsetting up sparsity pattern..." << std::endl;
  LoggerHelper help(gv.size(0), 5, 50000);

  const LeafIndexSet& set = gv.leafView().indexSet();
  LeafIterator itend = gv.leafView().template end<0>();

  // iterate over cells
  int cell=0;
  for (LeafIterator it = gv.leafView().template begin<0>(); it != itend; ++it, ++cell) {
    Dune::GeometryType gt = it->type();

#if DUNE_VERSION_NEWER(DUNE_GRID, 2, 3)
    const Dune::template ReferenceElement<double,dim>& ref =
      Dune::ReferenceElements<double,dim>::general(gt);
#else
    const Dune::template GenericReferenceElement<double,dim>& ref =
      Dune::GenericReferenceElements<double,dim>::general(gt);
#endif

    int vertexsize = ref.size(dim);
    for (int i=0; i < vertexsize; i++) {
      int indexi = set.subIndex(*it,i,dim);
      for (int k=0;k<dim;++k) {
        MPC* mpc = getMPC(indexi,k);
        if (mpc) {
          for (size_t l=0;l<mpc->getNoMaster();++l) {
            nodeAdjacency(it,vertexsize,
                          meqn[mpc->getMaster(l).node*dim+
                               mpc->getMaster(l).dof-1]);
          }
        } else
          nodeAdjacency(it,vertexsize,meqn[indexi*dim+k]);
      }
    }
    if (cell % 10000 == 0)
      help.log(cell, "\t\t... still processing ... cell ");
  }
}

}} // namespace Opm, Elasticity

#endif
