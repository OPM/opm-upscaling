//==============================================================================
//!
//! \file logutils.hpp
//!
//! \date Nov 9 2011
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Logging helper utilities
//!
//==============================================================================
#ifndef LOGUTILS_HPP_
#define LOGUTILS_HPP_

#include <iostream>

class LoggerHelper {
  public:
    LoggerHelper(int max_, int chunks, int minsize)
      : per_chunk(max_/chunks), max(max_)
    {
      if (max < minsize)
        per_chunk=-1;
    }

    void log(int it, const std::string& prefix)
    {
      if (per_chunk > 0 && it > 0 && (it % per_chunk) == 0 && max-it > per_chunk)
        std::cout << prefix << it  << '/' << max << std::endl;
    }
  protected:
    int per_chunk; //!< Will log for each per_chunk processed
    int max; //!< Total number of its
};

#endif
