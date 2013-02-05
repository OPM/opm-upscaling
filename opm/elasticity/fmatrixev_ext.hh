#ifndef DUNE_FMATRIXEIGENVALUES_EXT_HH
#define DUNE_FMATRIXEIGENVALUES_EXT_HH

namespace Opm {

namespace FMatrixHelp {

// defined in fmatrixev_ext.cpp
extern void eigenValuesNonsymLapackCall(
    const char* jobvl, const char* jobvr, const long
    int* n, double* a, const long int* lda, double* wr, double* wi, double* vl,
    const long int* ldvl, double* vr, const long int* ldvr, double* work,
    const long int* lwork, const long int* info);

template <int dim, typename K, class C>
static void eigenValuesNonSym(const FieldMatrix<K, dim, dim>& matrix,
                              Dune::FieldVector<C, dim>& eigenValues)
{
  {
    const long int N = dim ;
    const char jobvl = 'n';
    const char jobvr = 'n';

    // length of matrix vector 
    const long int w = N * N ;

    // matrix to put into dgeev 
    double matrixVector[dim * dim]; 

    // copy matrix  
    int row = 0;
    for(int i=0; i<dim; ++i) 
    {
      for(int j=0; j<dim; ++j, ++row) 
      {
        matrixVector[ row ] = matrix[ i ][ j ];
      }
    }

    // working memory 
    double eigenR[dim]; 
    double eigenI[dim];
    double work[3*dim];

    // return value information 
    long int info = 0;
    long int lwork = 3*dim;

    // call LAPACK routine (see fmatrixev_ext.cc)
    eigenValuesNonsymLapackCall(&jobvl, &jobvr, &N, &matrixVector[0], &N, 
                                &eigenR[0], &eigenI[0], 0, &N, 0, &N, &work[0],
                                &lwork, &info);

    if( info != 0 ) 
    {
      std::cerr << "For matrix " << matrix << " eigenvalue calculation failed! " << std::endl;
      DUNE_THROW(Dune::InvalidStateException,"eigenValues: Eigenvalue calculation failed!");
    }
    for (int i=0;i<N;++i) {
      eigenValues[i].real = eigenR[i];
      eigenValues[i].imag = eigenI[i];
    }
  }
}

}

}

#endif
