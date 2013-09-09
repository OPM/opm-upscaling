//==============================================================================
//!
//! \file boundarygrid.cpp
//!
//! \date Nov 9 2011
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Class describing 2D quadrilateral grids
//!
//==============================================================================
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "boundarygrid.hh"

#include <iostream>


namespace Opm {
namespace Elasticity {

BoundaryGrid BoundaryGrid::uniform(const FaceCoord& min, const FaceCoord& max,
                                   int k1, int k2, bool dc)
{
  double dx = (max[0]-min[0])/k1;
  double dy = (max[1]-min[1])/k2;
  BoundaryGrid result;
  result.fixNodes.resize((k1+1)*(k2+1));
  for (int i=0;i<k1;++i) {
    for (int j=0;j<k2;++j) {
      Quad q;
      int k = j*(k1+1)+i;
      result.fixNodes[k] = false;
      q.v[0].i = k;
      q.v[0].c[0] = min[0]+i*dx;
      q.v[0].c[1] = min[1]+j*dy;
      q.v[1].i = k+1;
      result.fixNodes[k+1] = false;
      q.v[1].c[0] = q.v[0].c[0]+dx;
      int y=dc?3:2;
      q.v[y].i = k+k1+2;
      result.fixNodes[k+k1+2] = false;
      q.v[y].c[0] = q.v[0].c[0]+dx;
      q.v[y].c[1] = q.v[0].c[1]+dy;
      y=dc?2:3;
      q.v[y].i = k+k1+1;
      result.fixNodes[k+k1+1] = false;
      q.v[y].c[0] = q.v[0].c[0];
      q.v[y].c[1] = q.v[0].c[1]+dy;
      if (!dc) {
        if (i == 0 && j == 0)
          result.fixNodes[q.v[0].i] = true;
        if (i == k1-1 && j == 0)
          result.fixNodes[q.v[1].i] = true;
        if (i == 0 && j == k2-1)
          result.fixNodes[q.v[3].i] = true;
        if (i == k1-1 && j == k2-1)
          result.fixNodes[q.v[2].i] = true;
      }
      result.add(q);
    }
  }

  result.nodes = (k1+1)*(k2+1);
  return result;
}

void BoundaryGrid::extract(FaceCoord& res,
                           const GlobalCoordinate& coord, int dir)
{
  int j=0;
  for (int i=0;i<3;++i) {
    if (i != dir)
      res[j++] = coord[i];
  }
}

void BoundaryGrid::extract(Vertex& res, const Dune::FieldVector<double,3>& coord, int dir)
{
  extract(res.c,coord,dir);
}

void BoundaryGrid::add(const Quad& quad)
{
  grid.push_back(quad);
  Quad& q = grid.back();

  // establish bounding box
  q.bb[0] = q.v[0].c;
  q.bb[1] = q.v[0].c;
  for (int k=1;k<4;++k) {
    if (q.v[k].c[0] < q.bb[0][0])
      q.bb[0][0] = q.v[k].c[0];
    if (q.v[k].c[0] > q.bb[1][0])
      q.bb[1][0] = q.v[k].c[0];
    if (q.v[k].c[1] < q.bb[0][1])
      q.bb[0][1] = q.v[k].c[1];
    if (q.v[k].c[1] > q.bb[1][1])
      q.bb[1][1] = q.v[k].c[1];
  }
  if ((q.bb[1][0]-q.bb[0][0])*(q.bb[1][1]-q.bb[0][1]) < 1.e-7)
    grid.pop_back();
}

//! \brief Check that two points are sufficiently equal
//! \param[in] x First point
//! \param[in] y Second point
//! \param[in] tol Tolerance of comparison
inline bool EQUAL2(const BoundaryGrid::FaceCoord& x,
                   const BoundaryGrid::FaceCoord& y, double tol)
{
  return hypot(x[0]-y[0],x[1]-y[1]) < tol;
}

bool BoundaryGrid::find(Vertex& res, const Vertex& coord) const
{
  // find first quad with coord within bounding box 
  std::vector<Quad>::const_iterator it = std::find_if(grid.begin(),grid.end(),
                                                      BoundedPredicate(coord.c));

  res.i = -1;
  while (it != grid.end()) {
    res.q = const_cast<Quad*>(&(*it));
    // check if we have exactly a node
    bool ok=false;
    for (int i=0;i<4;++i) {
      if (EQUAL2(coord.c,it->v[i].c,1.e-8)) {
        res.i = it->v[i].i;
        ok = true;
        break;
      }
    } 
    if (!ok && Q4inv(res.c,*it,coord.c,1.e-8,1.e-8) > 0) {
      ok = true;
    }

    if (ok)
      break;
    it = std::find_if(it+1,grid.end(),BoundedPredicate(coord.c));
  }
  if (it == grid.end()) {
    std::cout << " failed to locate " << coord.c << std::endl;
    assert(0);
  }
  return it != grid.end();
}

int BoundaryGrid::Q4inv(FaceCoord& res, const Quad& q,
                        const FaceCoord& coord,
                        double epsZero, double epsOut) const
{
  double A[4];
  double B[4];

  /* Find coefficients of the bi-linear equation */
  A[0] =        ( q.v[0].c[0] - q.v[1].c[0] +
                  q.v[2].c[0] - q.v[3].c[0]);
  A[1] =        ( -q.v[0].c[0]+q.v[1].c[0]);
  A[2] =        ( -q.v[0].c[0]+q.v[3].c[0]);
  A[3] = coord[0]-q.v[0].c[0];

  B[0] =        ( q.v[0].c[1] - q.v[1].c[1] +
                  q.v[2].c[1] - q.v[3].c[1]);
  B[1] =        ( -q.v[0].c[1]+q.v[1].c[1]);
  B[2] =        ( -q.v[0].c[1]+q.v[3].c[1]);
  B[3] = coord[1]-q.v[0].c[1];

  // We have to solve the following set of equations:

  //  A1*XI*ETA + A2*XI + A3*ETA = A4
  //  B1*XI*ETA + B2*XI + B3*ETA = B4

  //  The way that we may solve this nonlinear (in XI and ETA)
  //  set of equations depends on the coefficients Ai,Bi.
  //  The solution is unique for proper input.
  std::vector<double> xi;
  std::vector<double> eta;
  bilinearSolve(epsZero,1,A,B,xi,eta);

  // check that obtained solutions are inside element
  double tol = 1+epsOut;
  int nInside=0;
  for (size_t i=0;i<xi.size();++i) {
    if (xi[i] < tol && eta[i] < tol) {
      if (++nInside > 1) {
        std::cout << "multiple solutions" << std::endl;
        FaceCoord old = q.pos(res[0],res[1]);
        FaceCoord ny  = q.pos(xi[i],eta[i]);
        double d1 = hypot(coord[0]-ny[0],coord[1]-ny[1]);
        double d2 = hypot(coord[0]-old[0],coord[1]-old[1]);
        if (d2 < d1)
          continue;
      }
    } else if (nInside == 0) {
      if (i > 0) {
        FaceCoord old = q.pos(res[0],res[1]);
        FaceCoord ny  = q.pos(xi[i],eta[i]);
        double d1 = hypot(coord[0]-ny[0],coord[1]-ny[1]);
        double d2 = hypot(coord[0]-old[0],coord[1]-old[1]);
        if (d2 < d1)
          continue;
      }
    }
    res[0] = xi[i];
    res[1] = eta[i];
  }

  return nInside;
}

bool BoundaryGrid::bilinearSolve(double epsilon, double order,
                                 const double* A, const double* B,
                                 std::vector<double>& X,
                                 std::vector<double>& Y) const
{
  double tol = 0;
  // geometric tolerance ?
  for (int i=0;i<4;++i) {
    double det = fabs(A[i]);
    if (det > tol) tol = det;
    det = fabs(B[i]);
    if (det > tol) tol = det;
  }
  tol *= epsilon;

  double det = A[1]*B[2]-B[1]*A[2];
  if (fabs(A[0]) < tol && fabs(B[0]) < tol) {
    // linear eqs
    if (fabs(det) < tol*tol) return false;
    X.push_back((B[2]*A[3]-A[2]*B[3])/det);
    Y.push_back((-B[1]*A[3]+A[1]*B[3])/det);
    return true;
  }
  // second order
  double Q0 = B[2]*A[3]-A[2]*B[3];
  double Q1 = B[0]*A[3]-A[0]*B[3]-det;
  double Q2 = A[0]*B[1]-B[0]*A[1];
  std::vector<double> Z;
  cubicSolve(epsilon,0,Q2,Q1,Q0,Z);
  for (size_t i=0;i<Z.size();++i) {
    Q0 = A[0]*Z[i]+A[2];
    if (fabs(Q0) > tol) {
      X.push_back(Z[i]);
      Y.push_back((A[3]-A[1]*Z[i])/Q0);
    }
  }
  Q0 = B[1]*A[3]-A[1]*B[3];
  Q1 = B[0]*A[3]-A[0]*B[3]+det;
  Q2 = A[0]*B[2]-B[0]*A[2];
  Z.clear();
  cubicSolve(epsilon,0,Q2,Q1,Q0,Z);
  for (size_t i=0;i<Z.size();++i) {
    Q0 = A[0]*Z[i]+A[1];
    if (fabs(Q0) > tol) {
      size_t j=0;
      for (j=0;j<Y.size();++j)
        if (fabs(Y[j]-Z[i]) <= epsilon*order) break;
      if (j == Y.size()) {
        X.push_back((A[3]-A[2]*Z[i])/Q0);
        Y.push_back(Z[i]);
      }
    }
  }

  return X.size() > 0;
}

bool BoundaryGrid::cubicSolve(double eps, double A, double B, double C,
                              double D, std::vector<double>& X) const
{
  if (fabs(A) > eps) { // cubic
    double epsmall = pow(eps,6.f);
    double P = (C-B*B/(3*A))/(3*A);
    double Q = ((2*B*B/(27*A)-C/3)*B/A+D)/(2*A);
    double W = Q*Q+P*P*P;
    if (W <= -epsmall && P < 0) {
      double FI = acos(-Q/sqrt(-P*P*P));
      X.push_back( 2*sqrt(-P)*cos(FI/3));
      X.push_back(-2*sqrt(-P)*cos((FI+M_PI)/3));
      X.push_back(-2*sqrt(-P)*cos((FI-M_PI)/3));
    } else if (fabs(W) < epsmall && Q < 0) {
      X.push_back(2*pow(-Q,1.f/3));
      X.push_back(-.5f*X[0]);
      X.push_back(X[1]);
    } else if (W > -epsmall && Q+sqrt(W) < 0 && Q-sqrt(W) < 0) {
      X.push_back(pow(-Q+sqrt(W),1.f/3)+pow(-Q-sqrt(W),1.f/3));
      X.push_back(-.5f*X[0]);
      X.push_back(X[1]);
    } else if (W >= epsmall && fabs(Q) > epsmall && P > 0) {
      double FI = atan(sqrt(P*P*P)/Q);
      double KI = atan(pow(tan(.5f*FI),1.f/3));
      X.push_back(-2*sqrt(P)/tan(KI+KI));
      X.push_back(-.5f*X[0]);
      X.push_back(X[1]);
    } else if (W > -epsmall && fabs(Q) > epsmall && P < 0)
      return false;
    else
      return false;

    W = B/(3*A);
    X[0] -= W;
    X[1] -= W;
    X[2] -= W;

    return true;
  } else if (fabs(B) > eps) {
    double epsmall = pow(eps,4.f);
    double P = C*C-4*B*D;
    if (P > 0) {
      double Q = sqrt(P);
      X.push_back((-C+Q)/(B+B));
      X.push_back((-C-Q)/(B+B));
    } else if (P > -epsmall) {
      X.push_back(-C/(B+B));
      X.push_back(X[0]);
    }
    else
      return false;
  } else if (fabs(C) > eps) {
    X.push_back(-D/C);
  } else
    return false;

  return true;
}

BoundaryGrid::FaceCoord BoundaryGrid::Quad::pos(double xi, double eta) const
{
  BoundaryGrid::FaceCoord res;
  std::vector<double> N = evalBasis(xi,eta);

  res[0] = N[0]*v[0].c[0]+N[1]*v[1].c[0]+N[2]*v[2].c[0]+N[3]*v[3].c[0];
  res[1] = N[0]*v[0].c[1]+N[1]*v[1].c[1]+N[2]*v[2].c[1]+N[3]*v[3].c[1];

  return res;
}

std::vector<double> BoundaryGrid::Quad::evalBasis(double xi, double eta) const
{
  std::vector<double> res;
  res.push_back((1-xi)*(1-eta));
  res.push_back(xi*(1-eta));
  res.push_back(xi*eta);
  res.push_back((1-xi)*eta);

  return res;
}

BoundaryGrid::Vertex minXminY(std::vector<BoundaryGrid::Vertex>& in)
{
  // find the nodes with minimal X
  // then find the minimum Y among these 
  std::vector<BoundaryGrid::Vertex> s(in);
  std::sort(s.begin(),s.end(),BoundaryGrid::VertexLess(0));
  std::sort(s.begin(),s.begin()+2,BoundaryGrid::VertexLess(1));
  return *s.begin();
}

BoundaryGrid::Vertex maxXminY(std::vector<BoundaryGrid::Vertex>& in)
{
  // find the nodes with maximum X
  // then find the minimum Y among these 
  std::vector<BoundaryGrid::Vertex> s(in);
  std::sort(s.begin(),s.end(),BoundaryGrid::VertexLess(0));
  std::sort(s.begin()+2,s.end(),BoundaryGrid::VertexLess(1));
  return *(s.end()-2);
}

BoundaryGrid::Vertex maxXmaxY(std::vector<BoundaryGrid::Vertex>& in)
{
  // find the nodes with maximum X
  // then find the maximum Y among these 
  std::vector<BoundaryGrid::Vertex> s(in);
  std::sort(s.begin(),s.end(),BoundaryGrid::VertexLess(0));
  std::sort(s.begin()+2,s.end(),BoundaryGrid::VertexLess(1));
  return *(s.end()-1);
}

BoundaryGrid::Vertex minXmaxY(std::vector<BoundaryGrid::Vertex>& in)
{
  // find the nodes with minimum X
  // then find the maximum Y among these 
  std::vector<BoundaryGrid::Vertex> s(in);
  std::sort(s.begin(),s.end(),BoundaryGrid::VertexLess(0));
  std::sort(s.begin(),s.begin()+2,BoundaryGrid::VertexLess(1));
  return *(s.begin()+1);
}

}
}
