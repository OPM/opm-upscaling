//!
//! \file meshcolorizer.cpp
//!
//! \date Jul 26 2013
//!
//! \author Arne Morten Kvarving / SINTEF, Knut Morten Okstad / SINTEF
//!
//! \brief Mesh colorizer class - template specializations
//!
//==============================================================================
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "meshcolorizer.hpp"

#include <dune/grid/CpGrid.hpp>

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

  template<>
void MeshColorizer<Dune::CpGrid>::calcGroups()
{
  int nel1 = grid.logicalCartesianSize()[0];
  int nel2 = grid.logicalCartesianSize()[1];
  int nel3 = grid.logicalCartesianSize()[2];
  int threads=1;
  int stripsize=0;
  int remainder=0;
  int dir=0, mul=1;

#ifdef HAVE_OPENMP
  threads = omp_get_max_threads();
  int parts = threads > 1 ? 2*threads : 1;
  dir = getStripDirection(nel1,nel2,nel3,parts);
  mul = dir == 0 ? 1 : nel1*(dir == 1 ? 1 : nel2);
  int els = dir == 0 ? nel1 : (dir == 1 ? nel2 : nel3);

  while (threads > 1 && (stripsize = els/parts) < 2) {
    threads --;
    parts -= 2;
  }

  if (threads > 1)
    remainder = els - stripsize*parts;
  else
    stripsize = els;

  std::cout << "Grid partitioning summary:" << std::endl;
  std::cout << "\twe have " << threads << " threads available"
    << "\n\tstripsize " << stripsize
    << "\n\t# of strips " << els/stripsize
    << "\n\tremainder " << remainder << std::endl;
#endif

  if (threads == 1) {
    tg[0].resize(1);
    tg[0].reserve(grid.size(0));
    for (int i = 0; i < grid.size(0); ++i)
      tg[0][0].push_back(i);
    tg[1].clear();
  } else {
    // invert compressed cell array
    std::vector<int> globalActive;
    globalActive.resize(nel1*nel2*nel3, -1);
    for (size_t i=0;i<grid.globalCell().size();++i)
      globalActive[grid.globalCell()[i]] = i;
    IntVec stripsizes[2];
    stripsizes[0].resize(threads,stripsize);
    stripsizes[1].resize(threads,stripsize);
    for (int i = 1; i <= remainder; ++i)
      stripsizes[i%2][threads-(i+1)/2]++;

    IntVec startelms[2];
    int offs;
    for (int i = offs = 0; i < threads; ++i) {
      startelms[0].push_back(offs*mul);
      offs += stripsizes[0][i];
      startelms[1].push_back(offs*mul);
      offs += stripsizes[1][i];
    }

    for (int i = 0; i < 2; ++i) { // loop over groups
      tg[i].resize(threads);
      for (int t = 0; t < threads; ++t) { // loop over threads
        int maxx = dir == 0 ? stripsizes[i][t] : nel1;
        int maxy = dir == 1 ? stripsizes[i][t] : nel2;
        int maxz = dir == 2 ? stripsizes[i][t] : nel3;
        for (int i3 = 0; i3 < maxz; ++i3) {
          for (int i2 = 0; i2 < maxy; ++i2) {
            for (int i1 = 0; i1 < maxx; ++i1) {
              int elem = startelms[i][t]+i1+nel1*(i2+nel2*i3);
              if (globalActive[elem] > -1)
                tg[i][t].push_back(globalActive[elem]);
            }
          }
        }
      }

#ifdef VERBOSE
      std::cout << "group " << i << std::endl;
      for (size_t j = 0; j < tg[i].size(); ++j) {
        std::cout << "\t thread " << j << " (" << tg[i][j].size() << "): ";
        for (size_t k = 0; k < tg[i][j].size(); ++k)
          std::cout << tg[i][j][k] << " ";
        std::cout << std::endl;
      }
#endif
    }
  }
}
