// $Id: fvector.hh 5262 2008-09-07 09:03:38Z christi $
#ifndef OPM_MATVECTRAITS_HH
#define OPM_MATVECTRAITS_HH

namespace Opm {

/**
   @addtogroup DenseMatVec
   \brief Type Traits to retrieve types associated with an implementation of Opm::DenseVector or Opm::DenseMatrix

   you have to specialize this class for every implementation of DenseVector or DenseMatrix.

   \code
   //! export the type of the derived class (e.g. FieldVector<K,SIZE>)
   typedef ... derived_type;
   //! export the type of the stored values
   typedef ... value_type;
   //! export the type representing the size information
   typedef ... size_type;
   \endcode

*/
template<class T>
struct DenseMatVecTraits {};

} // end namespace Opm

#endif // OPM_FTRAITS_HH
