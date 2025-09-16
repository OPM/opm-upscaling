//==============================================================================
//!
//! \file mortar_utils.hpp
//!
//! \date Nov 9 2012
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Mortar helper class
//!
//==============================================================================
#ifndef MORTAR_UTILS_HPP_
#define MORTAR_UTILS_HPP_

#include <opm/elasticity/matrixops.hpp>

namespace Opm {
  namespace Elasticity {

class MortarUtils {
  public:
    //! \brief Extract a range of indices from a vector
    //! \param[out] x The vector with the extracted data
    //! \param[in] y The vector
    //! \param[in] len The number of indices to extract
    //! \param[in] start The first index in the range
    static void extractBlock(Vector& x, const Vector& y, int len, int start = 0)
    {
      x.resize(len);
      std::copy(y.begin() + start, y.begin() + len + start, x.begin());
    }

    //! \brief Inject a range of indices into a vector
    //! \param[in,out] x The vector to inject into
    //! \param[in] y The vector with the data to inject
    //! \param[in] len The number of indices to inject
    //! \param[in] start The first index in the range
    static void injectBlock(Vector& x, const Vector& y, int len, int start = 0)
    {
      std::copy(y.begin(), y.begin() + len, x.begin() + start);
    }
};

}
}

#endif
