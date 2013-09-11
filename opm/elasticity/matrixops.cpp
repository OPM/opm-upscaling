//==============================================================================
//!
//! \file matrixops.cpp
//!
//! \date Nov 9 2011
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Helper class with some matrix operations
//!
//==============================================================================

#if HAVE_CONFIG_H
#include "config.h"
#endif

#include "matrixops.hpp"

#include <iostream>
#include <fstream>

namespace Opm {
namespace Elasticity {

void MatrixOps::fromAdjacency(Matrix& A, const std::vector< std::set<int> >& adj,
                              int rows, int cols)
{
  size_t sum=0;
  for (size_t i=0;i<adj.size();++i)
    sum += adj[i].size();
  A.setSize(rows, cols, sum);
  A.setBuildMode(Matrix::random);

  for (size_t i = 0; i < adj.size(); i++)
    A.setrowsize(i,adj[i].size());
  A.endrowsizes();

  for (size_t i = 0; i < adj.size(); i++) {
    std::set<int>::iterator setend = adj[i].end();
    for (std::set<int>::iterator setit = adj[i].begin();
        setit != setend; ++setit) {
      A.addindex(i,*setit);
    }
  }
  A.endindices();
  A = 0;
}

Matrix MatrixOps::fromDense(const Dune::DynamicMatrix<double>& T)
{
  AdjacencyPattern a;
  a.resize(T.N());
  for (size_t i=0; i < T.N(); ++i) {
    for (size_t j=0; j < T.M(); ++j) {
      if (fabs(T[i][j]) > 1.e-14)
        a[i].insert(j);
    }
  }

  Matrix result;
  fromAdjacency(result, a, T.N(), T.M());

  for (Matrix::ConstIterator it  = result.begin();
                             it != result.end(); ++it) {
    for (Matrix::ConstColIterator it2  = it->begin();
                                  it2 != it->end();++it2)
      result[it.index()][it2.index()] = T[it.index()][it2.index()];
  }

  return result;
}

void MatrixOps::print(const Matrix& A)
{
  for (Matrix::ConstIterator it  = A.begin();
                             it != A.end(); ++it) {
    for (Matrix::ConstColIterator it2  = it->begin();
                                  it2 != it->end();++it2) {
      double val = *it2; 
      if (fabs(val) < 1.e-14)
        continue;
      std::cout << it.index() << " " << it2.index() << " : " << val << std::endl;
    }
  }
}

Matrix MatrixOps::Axpy(const Matrix& A, const Matrix& B, double alpha)
{
  assert(A.M() == B.M() && A.N() == B.N());

  // establish union adjacency pattern
  AdjacencyPattern adj;
  adj.resize(A.N());
  for (Matrix::ConstIterator it  = A.begin();
                             it != A.end(); ++it) {
    for (Matrix::ConstColIterator it2  = it->begin();
                                  it2 != it->end();++it2)
      adj[it.index()].insert(it2.index());
  }
  for (Matrix::ConstIterator it  = B.begin();
                             it != B.end(); ++it) {
    for (Matrix::ConstColIterator it2  = it->begin();
                                  it2 != it->end();++it2)
      adj[it.index()].insert(it2.index());
  }
  Matrix result;
  fromAdjacency(result,adj,A.N(),A.M());
  // now insert elements from A
  for (Matrix::ConstIterator it  = A.begin();
                             it != A.end(); ++it) {
    for (Matrix::ConstColIterator it2  = it->begin();
                                  it2 != it->end();++it2)
      result[it.index()][it2.index()] = *it2;
  }
  // and subtract elements from B
  for (Matrix::ConstIterator it  = B.begin();
                             it != B.end(); ++it) {
    for (Matrix::ConstColIterator it2  = it->begin();
                                  it2 != it->end();++it2)
      result[it.index()][it2.index()] += alpha*(*it2);
  }

  return result;
}

Matrix MatrixOps::augment(const Matrix& A, const Matrix& B,
                     size_t r0, size_t c0, bool symmetric)
{
  std::cout << "Augmenting matrix of dimension " << A.N() << "x" << A.M()
            << " with matrix of dimension " << B.N() << "x" << B.M() << std::endl;
  size_t nrow = A.N();
  size_t ncol = A.M();
  if (r0+B.N() > nrow) nrow = r0+B.N();
  if (symmetric && r0+B.N() > ncol) ncol = r0+B.N();
  if (c0+B.M() > ncol) ncol = c0+B.M();
  if (symmetric && c0+B.M() > nrow) nrow = c0+B.M();
  std::cout << "Resulting size: " << nrow << "x" << ncol << std::endl;

  AdjacencyPattern adj;
  adj.resize(nrow);
  for (Matrix::ConstIterator it  = A.begin();
                             it != A.end();++it) {
    for (Matrix::ConstColIterator it2  = it->begin(); 
                                  it2 != it->end();++it2) {
      adj[it.index()].insert(it2.index());
    }
  }
  for (Matrix::ConstIterator it  = B.begin();
                             it != B.end();++it) {
    for (Matrix::ConstColIterator it2  = it->begin(); 
                                  it2 != it->end();++it2) {
      adj[it.index()+r0].insert(it2.index()+c0);
      if (symmetric)
        adj[it2.index()+c0].insert(it.index()+r0);
    }
  }
  if (symmetric) {
    // always establish diagonal elements or superLU crashes
    for (size_t i=0;i<nrow;++i)
      adj[i].insert(i);
  }
  Matrix result;
  fromAdjacency(result,adj,nrow,ncol);
  for (Matrix::ConstIterator it  = A.begin();
                             it != A.end();++it) {
    for (Matrix::ConstColIterator it2  = it->begin(); 
                                  it2 != it->end();++it2) {
      result[it.index()][it2.index()] = *it2;
    }
  }
  for (Matrix::ConstIterator it  = B.begin();
                             it != B.end();++it) {
    for (Matrix::ConstColIterator it2  = it->begin(); 
                                  it2 != it->end();++it2) {
      result[it.index()+r0][it2.index()+c0] = *it2;
      if (symmetric)
        result[it2.index()+c0][it.index()+r0] = *it2;
    }
  }
  
  return result;
}

Matrix MatrixOps::extractDiagonal(const Matrix& A)
{
  AdjacencyPattern adj;
  adj.resize(A.M());
  for (size_t i=0;i<A.M();++i)
    adj[i].insert(i);
  Matrix result;
  fromAdjacency(result,adj,A.M(),A.M());
  for (size_t i=0;i<A.M();++i)
    result[i][i] = A[i][i];

  return result;
}

Matrix MatrixOps::diagonal(size_t N)
{
  AdjacencyPattern adj;
  adj.resize(N);
  for (size_t i=0;i<N;++i)
    adj[i].insert(i);
  Matrix result;
  fromAdjacency(result,adj,N,N);

  return result;
}

void MatrixOps::saveAsc(const Matrix& A, const std::string& file)
{
  std::ofstream f;
  f.open(file.c_str());
  f << "% " << A.N() << " " << A.M() << std::endl;
  int prevrow=-1;
  for (Matrix::ConstIterator it  = A.begin();
                             it != A.end(); ++it) {
    for (int i=0;i<int(it.index())-prevrow-1;++i) {
      for (size_t j=0;j<A.M();++j)
        f << "0 ";
      f << std::endl;
    }
    int prevcol=-1;
    for (Matrix::ConstColIterator it2  = it->begin();
                                  it2 != it->end();++it2) {
      for (int j=0;j<int(it2.index())-prevcol-1;++j)
        f << "0 ";
      double val = *it2;
      f << val << " ";
      prevcol = it2.index();
    }
    for (int j=0;j<int(A.M())-prevcol-1;++j)
      f << "0 ";
    prevrow = it.index();
    f << std::endl;
  }
  f.close();
}

Matrix MatrixOps::extractBlock(const Matrix& A, size_t r0, size_t N,
                               size_t c0, size_t M)
{
  // establish adjacency pattern
  AdjacencyPattern adj;
  adj.resize(N);
  for (Matrix::ConstIterator it  = A.begin();
                             it != A.end(); ++it) {
    if (it.index() < r0 || it.index() >= N+r0)
      continue;
    for (Matrix::ConstColIterator it2  = it->begin();
                                  it2 != it->end(); ++it2) {
      if (it2.index() < c0 || it2.index() >= M+c0)
        continue;
      adj[it.index()-r0].insert(it2.index()-c0);
    }
  }

  Matrix result;
  fromAdjacency(result, adj, N, M);

  // now insert elements from A
  for (Matrix::ConstIterator it  = A.begin();
                             it != A.end(); ++it) {
    if (it.index() < r0 || it.index() >= N+r0)
      continue;
    for (Matrix::ConstColIterator it2  = it->begin();
                                  it2 != it->end(); ++it2) {
      if (it2.index() < c0 || it2.index() >= M+c0)
        continue;
      result[it.index()-r0][it2.index()-c0] = *it2;
    }
  }

  return result;
}

}}
