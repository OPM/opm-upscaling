//==============================================================================
//!
//! \file boundarygrid.hh
//!
//! \date Nov 9 2011
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Class describing 2D quadrilateral grids
//!
//==============================================================================
#ifndef BOUNDARYGRID_HH_
#define BOUNDARYGRID_HH_

#include <dune/common/version.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/genericgeometry/matrixhelper.hh>
#include <dune/grid/common/mcmgmapper.hh>

#include <vector>

namespace Opm {
namespace Elasticity {

//! \brief A class describing a quad grid
class BoundaryGrid {
  public:
    //! \brief A coordinate on the underlying 3D grid
    typedef Dune::FieldVector<double,3> GlobalCoordinate;
    //! \brief A coordinate on the quad grid
    typedef Dune::FieldVector<double,2> FaceCoord;

    //! \brief Establish an uniform quad grid
    //! \param[in] min Lower left corner
    //! \param[in] max Upper right corner
    //! \param[in] k1 Number of elements in the first direction
    //! \param[in] k2 Number of elements in the second direction
    //! \param[in] dc If true, order quads according to dune conventions
    //! \returns A quad grid spanning the given area. Nodes are numbered
    //              in natural order along the first direction
    static BoundaryGrid uniform(const FaceCoord& min, const FaceCoord& max,
                                int k1, int k2, bool dc=false);

    //! \brief Default (empty) constructor
    BoundaryGrid() {}

    //! \brief Default (empty) destructor
    virtual ~BoundaryGrid() {}

    //! \brief Holds the indices and relevant coordinates of the vertices 
    //         on a boundary
    class Quad;

    //! \brief A class describing a 2D vertex
    class Vertex {
      public:
        //! \brief Default constructor
        Vertex() : i(-1), c(0), q(0), fixed(false) {}

        //! \brief Index of the vertex
        int i;

        //! \brief Coordinates of the vertex (2D)
        FaceCoord c;

        //! \brief The quad containing the vertex (if search has been done)
        Quad* q;

        //! \brief Whether or not this node is fixed
        bool fixed;

        //! \brief Comparison operator
        bool operator==(const Vertex& v2)
        {
          return hypot(v2.c[0]-c[0],v2.c[1]-c[1]) < 1.e-8;
        }
    };
    //! \brief A class describing a linear, quadrilateral element
    class Quad {
      public:
        //! \brief Default constructor
        Quad()
        {
          v[0].i = v[1].i = v[2].i = v[3].i = 0;
          v[0].c = v[1].c = v[2].c = v[3].c = 0.f;
        }
        //! \brief Return the physical coordinates corresponding to the
        //!        given local coordinates
        //! \param[in] xi The local coordinate in the first direction
        //! \param[in] eta The local coordinate in the second direction
        FaceCoord pos(double xi, double eta) const;

        //! \brief Evaluate the basis functions
        //! \param[in] xi The local coordinate in the first direction
        //! \param[in] eta The local coordinate in the second direction
        std::vector<double> evalBasis(double xi, double eta) const;

        //! \brief Vertices
        Vertex v[4];
        //! \brief Bounding box
        FaceCoord bb[2];
      protected:
        //! \brief Print to a stream
        friend std::ostream& operator <<(std::ostream& os, const Quad& q)
        {
          os << "Nodes: " << q.v[0].i << "," << q.v[1].i << "," << q.v[2].i << "," << q.v[3].i << std::endl;
          os << "(" << q.v[0].c << ")(" << q.v[1].c << ")(" << q.v[2].c << ")(" << q.v[3].c << ")";
          return os;
        }
    };

    //! \brief Add a quad to the grid
    //! \param[in] quad The quad to add
    void add(const Quad& quad);

    void addToColumn(size_t col, const Quad& quad)
    {
      if (col >= colGrids.size())
        colGrids.resize(col+1);
      colGrids[col].push_back(quad);
    }

    //! \brief Obtain a reference to a quad
    //! \param[in] index The index of the requested quad
    //! \returns A reference to the requested quad
    Quad& operator[](int index)
    {
      return grid[index];
    }

    //! \brief Obtain a const reference to a quad
    //! \param[in] index The index of the requested quad
    //! \returns A const reference to the requested quad
    const Quad& operator[](int index) const
    {
      return grid[index];
    }

    const Quad& getQuad(int col, int index) const
    {
      return colGrids[col][index];
    }

    //! \brief Obtain the number of quads in the grid
    size_t size() const
    {
      return grid.size();
    }

    size_t colSize(int i) const
    {
      return colGrids[i].size();
    }

    //! \brief Return the total number of nodes on the grid when known
    //! \sa uniform
    size_t totalNodes() const
    {
      return nodes;
    }

    //! \brief Check if a given node is marked as fixed
    //! \param[in] node The requested node
    //! \returns Whether or not the node is marked as fixed
    bool isFixed(int node) const
    {
      return fixNodes[node];
    }

    //! \brief Locate the position of a vertex on the grid
    //! \param[in] coord The coordinate of the vertex
    //! \param[out] res The resulting coordinates
    bool find(Vertex& res, const Vertex& coord) const;

    //! \brief Helper function for extracting given 2D coordinates from a 3D coordinate
    //! \param[in] coord The 3D coordinates of the vertex
    //! \param[in] dir The direction to ignore
    //! \param[out] res The resulting coordinates
    static void extract(FaceCoord& res,
                        const GlobalCoordinate& coord, int dir);

    //! \brief Helper function for extracting given 2D coordinates from a 3D vector
    //! \param[in] coord The 3D coordinates of the vertex
    //! \param[in] dir The direction to ignore
    //! \param[out] res The resulting coordinates
    static void extract(Vertex& res,
                        const Dune::FieldVector<double,3>& coord, int dir);

    //! \brief Predicate for sorting vertices
    struct VertexLess {
        //! \brief Default constructor.
        //! \param[in] comp Direction to use for comparison. -1 to use index
        VertexLess(int comp) : dir(comp) {}

        //! \brief The comparison operator
        bool operator()(const Vertex& q1, const Vertex& q2)
        {
          if (dir >= 0)
            return q1.c[dir] < q2.c[dir];
          return q1.i < q2.i;
        }

        //! \brief Compare using this direction, if -1 compare using index values
        int dir;
    };

    //! \brief Predicate for checking if a vertex falls within a quads bounding box
    struct BoundedPredicate {
      //! \brief Default constructor
      //! \param[in] coord_ The coordinates to check
      BoundedPredicate(const FaceCoord& coord_) : coord(coord_) {}

      //! \brief The comparison operator
      bool operator()(const Quad& q)
      {
        double eps = 1.e-8;
        return (coord[0] >= q.bb[0][0]-eps && 
                coord[0] <= q.bb[1][0]+eps &&
                coord[1] >= q.bb[0][1]-eps &&
                coord[1] <= q.bb[1][1]+eps);
      }

      //! \brief The coordinates to check
      FaceCoord coord;
    };
  protected:
    //! \brief Our quadrilateral elements
    std::vector<Quad> grid;
    std::vector<std::vector<Quad> > colGrids;

    //! \brief Whether or not a given node is marked as fixed
    std::vector<bool> fixNodes;

    //! \brief Total number of nodes on grid
    size_t nodes;

    //! \brief Print to a stream
    friend std::ostream& operator <<(std::ostream& os, const BoundaryGrid& g)
    {
      for (size_t i=0;i<g.size();++i)
        os << g[i] << std::endl;
      return os;
    }

    //! \brief Solves a bi-linear set of equations in x and y.
    //!        A1 * x*y  +  A2 * x  +  A3 * y  =  A4
    //!        B1 * x*y  +  B2 * x  +  B3 * y  =  B4
    //! \param[in] epsilon The tolerance for equality checks with zero
    //! \param[in] order The expected order of the solution (used for unique checks)
    //! \param[in] A The coefficients of the first equation
    //! \param[in] B The coefficients of the second equation
    //! \param[out] X The first component of the solutions
    //! \param[out] Y The second component of the solutions
    bool bilinearSolve(double epsilon, double order,
                       const double* A, const double* B,
                       std::vector<double>& X,
                       std::vector<double>& Y) const;

    //! \brief Solves the cubic equation A*x^3+B*x^2+C*x+D=0
    //! \param[in] eps The tolerance for equality checks with zero
    //! \param[in] A Equation coefficient
    //! \param[in] B Equation coefficient
    //! \param[in] C Equation coefficient
    //! \param[in] D Equation coefficient
    //! \param[out] X The obtained solutions
    bool cubicSolve(double eps, double A, double B, double C,
                    double D, std::vector<double>& X) const;

    //! \brief Find the local coordinates of a given point within a given quad
    //! \param[in] q The quad to search within
    //! \param[in] coord The coordinates to search for
    //! \param[in] epsZero The tolerance for equality checks with zero
    //! \param[in] epsOut The tolerance check for outside checks
    //! \param[out] res The obtained result
    int Q4inv(FaceCoord& res, const Quad& q, const FaceCoord& coord,
              double epsZero, double epsOut) const;
};

//! \brief Geometry class for general hexagons
  template<int mydim, int cdim, class GridImp>
class HexGeometry
{
};

//! \brief Specialization for 2D quadrilaterals
  template<int cdim, class GridImp>
class HexGeometry<2, cdim, GridImp>
{
  public:
    //! \brief The dimension of the grid.
    enum {   dimension = 2};

    //! \brief Dimension of the domain space
    enum { mydimension = 2};

    //! \brief Dimension of the range space
    enum { coorddimension = cdim };

    //! \brief World dimension of underlying grid
    enum {dimensionworld = 2 };

    //! \brief Coordinate element type
    typedef double ctype;

    //! \brief Domain type
    typedef Dune::FieldVector<ctype,mydimension> LocalCoordinate;

    //! \brief Range type
    typedef Dune::FieldVector<ctype,coorddimension> GlobalCoordinate;

    //! \brief Type of Jacobian matrix
    typedef Dune::FieldMatrix<ctype,coorddimension,mydimension> Jacobian;

    //! \brief Type of transposed Jacobian matrix
    typedef Dune::FieldMatrix<ctype,mydimension,coorddimension> JacobianTransposed;

    //! \brief Construct integration element extracted from a 3D grid
    //! \param[in] q Quad describing element
    //! \param[in] gv Underlying 3D grid quads are extracted from
    //! \param[in] dir The direction of the normal vector on the face
    HexGeometry(const BoundaryGrid::Quad& q, const GridImp& gv, int dir)
    {
      Dune::LeafMultipleCodimMultipleGeomTypeMapper<GridImp,
                                           Dune::MCMGVertexLayout> mapper(gv);
      typename GridImp::LeafGridView::template Codim<3>::Iterator start=gv.leafGridView().template begin<3>();
      const typename GridImp::LeafGridView::template Codim<3>::Iterator itend = gv.leafGridView().template end<3>();
      for (int i=0;i<4;++i) {
        typename GridImp::LeafGridView::template Codim<3>::Iterator it=start;
        for (; it != itend; ++it) {
          if (mapper.map(*it) == q.v[i].i)
            break;
        }
        BoundaryGrid::extract(c[i],it->geometry().corner(0),dir);
      }

      m_volume = (c[1][0]-c[0][0])*(c[2][1]-c[0][1]);
    }

    //! \brief Construct integration element
    //! \param[in] q Quad describing element
    HexGeometry(const BoundaryGrid::Quad& q)
    {
      for (int i=0;i<4;++i)
        c[i] = q.v[i].c;
      m_volume = (c[1][0]-c[0][0])*(c[2][1]-c[0][1]);
    }

    //! \brief Returns entity type (a 2D cube)
    Dune::GeometryType type() const
    {
      Dune::GeometryType t;
      t.makeCube(mydimension);
      return t;
    }

    //! \brief Returns number of corners
    int corners() const
    {
      return 4;
    }

    //! \brief Returns volume (area) of quadrilateral
    ctype volume() const
    {
      return m_volume;
    }

    //! \brief Returns center of quadrilateral
    GlobalCoordinate center() const
    {
      LocalCoordinate local;
      local = .5f;
      return Global(local);
    }

    //! \brief Returns coordinates to requested corner
    //! \param[in] cor The requested corner (0..3)
    GlobalCoordinate corner(int cor) const
    {
      return c[cor];
    }

    //! \brief Map from local coordinates to global coordinates
    //! \param[in] local The local coordinates
    GlobalCoordinate global(const LocalCoordinate& local) const
    {
      // uvw = { (1-u, 1-v, 1-w), (u, v, w) }
      LocalCoordinate uvw[2] = { LocalCoordinate(1.0), local };
      uvw[0] -= local;
      // Access pattern for uvw matching ordering of corners.
      const int pat[4][2] = {{ 0, 0 },
                             { 1, 0 },
                             { 0, 1 },
                             { 1, 1 }};
      GlobalCoordinate xyz(0.0);
      for (int i = 0; i < 4; ++i) {
        GlobalCoordinate corner_contrib = corner(i);
        double factor = 1.0;
        for (int j = 0; j < 2; ++j) {
          factor *= uvw[pat[i][j]][j];
        }
        corner_contrib *= factor;
        xyz += corner_contrib;
      }
      return xyz;
    }

    //! \brief Map from global coordinates to local coordinates
    //! \param[in] y The global coordinates
    LocalCoordinate local(const GlobalCoordinate& y) const
    {
      const ctype epsilon = 1e-10;
      const Dune::ReferenceElement< ctype , 2 > & refElement =
        Dune::ReferenceElements< ctype, 2 >::general(type());
      LocalCoordinate x = refElement.position(0,0);
      LocalCoordinate dx;
      do {
        using namespace Dune::GenericGeometry;
        // DF^n dx^n = F^n, x^{n+1} -= dx^n
        JacobianTransposed JT = jacobianTransposed(x);
        GlobalCoordinate z = global(x);
        z -= y;
        MatrixHelper<DuneCoordTraits<double> >::template xTRightInvA<2, 2>(JT, z, dx );
        x -= dx;
      } while (dx.two_norm2() > epsilon*epsilon);
      return x;
    }

    //! \brief Return the transposed jacobian
    //! \param[in] local The local coordinates
    const Dune::FieldMatrix<ctype, mydimension, coorddimension>
      jacobianTransposed(const LocalCoordinate& local) const
    {
      // uvw = { (1-u, 1-v, (u, v) }
      LocalCoordinate uvw[2] = { LocalCoordinate(1.0), local };
      uvw[0] -= local;
      // Access pattern for uvw matching ordering of corners.
      const int pat[4][3] = {{ 0, 0 },
                             { 1, 0 },
                             { 0, 1 },
                             { 1, 1 }};
      Dune::FieldMatrix<ctype, mydimension, coorddimension> Jt(0.0);
      for (int i = 0; i < 4; ++i) {
        for (int deriv = 0; deriv < 2; ++deriv) {
          // This part contributing to dg/du_{deriv}
          double factor = 1.0;
          for (int j = 0; j < 2; ++j) {
            factor *= (j != deriv) ? uvw[pat[i][j]][j]
              : (pat[i][j] == 0 ? -1.0 : 1.0);
          }
          GlobalCoordinate corner_contrib = corner(i);
          corner_contrib *= factor;
          Jt[deriv] += corner_contrib; // using FieldMatrix row access.
        }
      }
      return Jt;
    }

    //! \brief Returns the inverse, transposed Jacobian
    //! \param[in] local The local coordinates
    const Dune::FieldMatrix<ctype, coorddimension, mydimension> 
      jacobianInverseTransposed(const LocalCoordinate& local) const
    {
      Dune::FieldMatrix<ctype, coorddimension, mydimension> Jti = jacobianTransposed(local);
      Jti.invert();
      return Jti;
    }

    //! \brief Returns the integration element (|J'*J|)^(1/2)
    //! \param[in] local The local coordinates
    ctype integrationElement(const LocalCoordinate& local) const
    {
      Dune::FieldMatrix<ctype, coorddimension, mydimension> Jt = jacobianTransposed(local);
      using namespace Dune::GenericGeometry;
      return MatrixHelper<DuneCoordTraits<double> >::template sqrtDetAAT<2, 2>(Jt);
    }
  private:
    //! \brief The coordinates of the corners
    GlobalCoordinate c[4];

    //! \brief The volume (area) of the quadrilateral
    ctype m_volume;
};

//! \brief Find the vertex in the vector with minimum X and minimum Y
//! \returns The requested vertex
BoundaryGrid::Vertex minXminY(std::vector<BoundaryGrid::Vertex>& in);

//! \brief Find the vertex in the vector with maximum X and minimum Y
//! \returns The requested vertex
BoundaryGrid::Vertex maxXminY(std::vector<BoundaryGrid::Vertex>& in);

//! \brief Find the vertex in the vector with maximum X and maximum Y
//! \returns The requested vertex
BoundaryGrid::Vertex maxXmaxY(std::vector<BoundaryGrid::Vertex>& in);

//! \brief Find the vertex in the vector with minimum X and maximum Y
//! \returns The requested vertex
BoundaryGrid::Vertex minXmaxY(std::vector<BoundaryGrid::Vertex>& in);

}
}

#endif
