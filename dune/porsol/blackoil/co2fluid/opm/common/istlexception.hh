#ifndef OPM_ISTLEXC_HH
#define OPM_ISTLEXC_HH

#include "exceptions_.hh"

namespace Opm {
   
    /** 
		@addtogroup ISTL
		@{
     */

  //! derive error class from the base class in common
  class ISTLError : public Opm::MathError {};

  /** @} end documentation */

} // end namespace

#endif
