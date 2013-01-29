/*
  Copyright 2011 IRIS - International Research Institute of Stavanger.

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


#ifndef IRISDuneISTLInterfaceH
#define IRISDuneISTLInterfaceH

#include<dune/common/fmatrix.hh>
#include<dune/istl/bcrsmatrix.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>

class IRISDuneISTLInterface
{
  public:
    typedef Dune::FieldMatrix<double, 1, 1> MatrixBlock;
    typedef Dune::BCRSMatrix<MatrixBlock> Matrix;
    typedef Dune::FieldVector<double, 1> VectorBlock;
    typedef Dune::BlockVector<VectorBlock> Vector; 
    
    IRISDuneISTLInterface(int nRow);
    
    template<class IRIDune::SGrid>
    void setupMatrix(IRIDune::SGrid& grid, const int& codim);

    void reset();
    
    Matrix::block_type& matrixElement(Matrix::size_type i, Matrix::size_type j) {return mtx_[i][j];}
    
    Vector::block_type& rhsElement(Vector::size_type i) {return rhs_[i];}
    
    Vector::block_type& solutionElement(Vector::size_type i) {return xx_[i];}
    
    void solveLinearSystem();
  
  private:
    int nRow_;
    Matrix mtx_;
    Vector rhs_;
    Vector xx_;
    	  
};

IRISDuneISTLInterface::IRISDuneISTLInterface(int nRow)
:nRow_(nRow)
,mtx_(nRow,nRow,IRISDuneISTLInterface::Matrix::random)
,rhs_(nRow)
,xx_(nRow)
{
}

template<class IRIDune::SGrid>
void IRISDuneISTLInterface::setupMatrix(IRIDune::SGrid& grid, const int& codim)
{
  if (codim == 0)
  {
    // Generate look-up table of element to element neighbour relations ...
    std::vector<std::set<int> > elm2elmIndx(nRow_);
    for (typename IRIDune::SGrid::ElementIterator it = grid.template setEntityPointerToFirst<typename IRISGrid::ElementIterator>(0); 
	                                it != grid.template entityPointerIsAtEnd<typename IRIDune::SGrid::ElementIterator>(0); ++it)
	 {
	   typename IRIDune::SGrid::ElementPointer ep = it;
	   int eCurrentIndx = grid.getGlobalIndexFromEntityPointer(ep);
	   int nVrtx = Dune::GenericReferenceElements<typename IRIDune::SGrid::ctype,IRISGrid::dim>::general(ep->type()).size(IRISGrid::VertexPointer::codim);
	   for (int v=0; v<nVrtx; ++v)
	   {
	     typename IRIDune::SGrid::VertexPointer vp = ep->template subEntity<IRISGrid::VertexPointer::codim>(v);
	     std::vector<typename IRIDune::SGrid::ElementPointer> elemNgbs;
	     grid.getNeighbours(vp,elemNgbs);
	     for (unsigned int e=0; e<elemNgbs.size(); ++e)
	     {
	       int eNeighbourIndx = grid.getGlobalIndexFromEntityPointer(elemNgbs[e]);
	       if (elm2elmIndx[eNeighbourIndx].count(eCurrentIndx))
	         continue;
	       elm2elmIndx[eNeighbourIndx].insert(eCurrentIndx);
	     }
	   } 
	 } 
	 
	 // Matrix generation - First pass - Allocate row sizes:
	 for (unsigned int i=0; i<elm2elmIndx.size(); ++i)
	 {
	   mtx_.setrowsize(i,elm2elmIndx[i].size());
	 }
	 mtx_.endrowsizes(); // finalize row setup phase
	 
	 // Matrix generation - Second pass - Column entries: 
	 for (unsigned int i=0; i<elm2elmIndx.size(); ++i)
	 {
	   for (std::set<int>::iterator itSet=elm2elmIndx[i].begin(); itSet!=elm2elmIndx[i].end(); ++itSet)
	   {
	     mtx_.addindex(i,*itSet);
	   }
	 }
	 mtx_.endindices(); // finalize column setup phase
	 

	 /*	 
	 // Print test:
	 for (unsigned int i=0; i<mtx_.N(); ++i)
	 {
	   for (unsigned int j=0; j<mtx_.M(); ++j)
	   {
	     if (mtx_.exists(i,j))
	       std::cout << "*";
	     else
	       std::cout << "-";
	   }
	   std::cout << std::endl;
	 }
	 */
  }
  else
  {
    assert(false); // Currently only implemented for codim=0
  }
}



void IRISDuneISTLInterface::reset()
{
  //Simply makes sure that all (relevant) elements in the matrix, right hand side
  //and solution vector are put to zero.

  mtx_ *= 0.0;
  rhs_ *= 0.0;
  xx_ *= 0.0;
}


void
IRISDuneISTLInterface::solveLinearSystem()
{
  Dune::SeqILU0<Matrix, Vector, Vector> pre(mtx_, double(1.0));
  Dune::MatrixAdapter<Matrix,Vector,Vector> op(mtx_);
  Dune::RestartedGMResSolver<Vector> sol(op,pre,1.0e-6,10,100,2);      
  Dune::InverseOperatorResult res;      
  sol.apply(xx_,rhs_,res);
	      
  std::cout << "iterations: " << res.iterations
           << "\nreduction: " << res.reduction
           << "\nconv_rate: " << res.conv_rate
           << "\nelapsed: " << res.elapsed << std::endl;
  std::cout << "\n xx: " << std::endl;         
  for (Vector::size_type i=0; i<xx_.size(); ++i)
  {
    std::cout << xx_[i] << std::endl;
  }      
}

#endif
