//==============================================================================
//!
//! \file matrixops_impl.hpp
//!
//! \date Nov 9 2011
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Helper class with some matrix operations - template implementations
//!
//==============================================================================

void MatrixOps::fromAdjacency(Matrix& A, const std::vector< std::set<int> >& adj,
                              int rows, int cols)
{  
  size_t sum=0;
  for (size_t i=0;i<adj.size();++i)
    sum += adj[i].size();
  A.setSize(rows, cols, sum);
  A.setBuildMode(Matrix::random);

  for (int i = 0; i < adj.size(); i++)
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
  std::vector<std::set<int> > adj;
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
    for (int i=0;i<nrow;++i)
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
