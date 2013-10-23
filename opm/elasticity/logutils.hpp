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
#ifdef HAVE_OPENMP
#include <omp.h>
#endif

class LoggerHelper {
  public:
    LoggerHelper(int max_, int chunks, int minsize)
      : per_chunk(max_/chunks), max(max_)
    {
      std::vector<int> sizes;
      for (int i=0;i<chunks;++i)
        sizes.push_back(max/chunks);
      for(int i=chunks-max%chunks;i<chunks;++i)
        sizes[i]++;
      groups.resize(chunks);
      groups[0] = 0;
      for (int i=1;i<chunks;++i)
        groups[i] = groups[i-1]+sizes[i-1];
      groups.push_back(max);
      if (max < minsize)
        per_chunk=-1;
    }

    std::pair<int, int> group(int group)
    {
      return std::make_pair(groups[group],groups[group+1]);
    }

    void log(int it, const std::string& prefix)
    {
      if(per_chunk == -1)
        return;
#ifdef HAVE_OPENMP
      if (omp_get_num_threads() == 1)
#endif
        std::cout << prefix << it  << '/' << max << std::endl;
    }
  protected:
    std::vector<int> groups; //!< Group start/end offsets
    int per_chunk; //!< Will log for each per_chunk processed
    int max; //!< Total number of its
};

#endif
