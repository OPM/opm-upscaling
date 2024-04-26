//==============================================================================
//!
//! \file meshcolorizer.hpp
//!
//! \date Jul 26 2013
//!
//! \author Arne Morten Kvarving / SINTEF, Knut Morten Okstad / SINTEF
//!
//! \brief Mesh colorizer class
//!
//==============================================================================

#ifndef MESHCOLORIZER_HPP_
#define MESHCOLORIZER_HPP_

#include <vector>

/*! \brief Generate a coloring of a mesh suitable for threaded assembly.
 *         The mesh is assumed structured, and is sliced in strips of
 *         alternating colors (two colors). The two color groups can then
 *         be assembled in parallel since they will never add to the same
 *         DOFs in the linear system. Currently it can only be instanced
 *         for a CpGrid.
 */
  template<class GridType>
class MeshColorizer {
  typedef std::vector<int> IntVec; //!< A vector of integers
  typedef std::vector<IntVec> IntMat; //!< A matrix of integers

  public:
    //! \brief Default constructor
    //! \param[in] grid_ The grid to colorize
    MeshColorizer(const GridType& grid_) :
      grid(grid_)
    {
      calcGroups();
    }

    //! \brief Return a color group
    //! \param[in] i The color to return, either 0 or 1
    const IntMat& operator[](unsigned int i)
    {
      return tg[i];
    }

    //! \brief Calculate the coloring
    void calcGroups(){
        tg[0].resize(1);
        tg[0].reserve(grid.size(0));
        for (int i = 0; i < grid.size(0); ++i)
            tg[0][0].push_back(i);
        tg[1].clear();
      
    }
  private:
    IntMat tg[2]; //!< The color groups
    const GridType& grid; //!< Reference to grid being colored

    //! \brief Find the optimal stripping direction (x, y or z)
    //! \param[in] nel1 Number of elements in x direction
    //! \param[in] nel2 Number of elements in y direction
    //! \param[in] nel3 Number of elements in z direction
    //! \param[in] parts Number of parts (threads) to use
    int getStripDirection (int nel1, int nel2, int nel3, int parts)
    {
      int s1 = nel1 / parts;
      int s2 = nel2 / parts;
      int s3 = nel3 / parts;
      int r1 = nel1 - s1*parts;
      int r2 = nel2 - s2*parts;
      int r3 = nel3 - s3*parts;

      if (r1*nel2*nel3 < nel1*r2*nel3 && r1*nel2*nel3 < nel1*nel2*r3)
        return 0; // strips along x axis
      else if (nel1*r2*nel3 < r1*nel2*nel3 && nel1*r2*nel3 < nel1*nel2*r3)
        return 1; // strips along y axis
      else if (nel1*nel2*r3 < r1*nel2*nel3 && nel1*nel2*r3 < nel1*r2*nel3)
        return 2; // strips along z axis

      // The number of left-over elements is not smallest in one direction only
      if (r1*nel2*nel3 > nel1*r2*nel3)
        return nel2 > nel3 ? 1 : 2;
      else if (nel1*r2*nel3 > nel1*nel2*r3)
        return nel1 > nel3 ? 0 : 2;
      else if (nel1*nel2*r3 > r1*nel2*nel3)
        return nel1 > nel2 ? 0 : 1;

      // The number of left-over elements is the same in all three directions
      if (nel1 >= nel2 && nel1 >= nel3)
        return 0;
      else if (nel2 >= nel1 && nel2 >= nel3)
        return 1;
      else
        return 2;
    }
};

#endif
