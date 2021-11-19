//==============================================================================
//!
//! \file elasticity_upscale_impl.hpp
//!
//! \date Nov 9 2011
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Elasticity upscale class - template implementations
//!
//==============================================================================
#ifndef OPM_ELASTICITY_UPSCALE_IMPL_HPP
#define OPM_ELASTICITY_UPSCALE_IMPL_HPP

#include <iostream>

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include <opm/parser/eclipse/Deck/DeckKeyword.hpp>

namespace Opm {
namespace Elasticity {

#undef IMPL_FUNC
#define IMPL_FUNC(A,B) template<class GridType, class PC> \
                         A ElasticityUpscale<GridType, PC>::B

IMPL_FUNC(std::vector<BoundaryGrid::Vertex>, 
          extractFace(Direction dir, ctype coord))
{
  std::vector<BoundaryGrid::Vertex> result;
  const LeafVertexIterator itend = gv.leafGridView().template end<dim>();

  // make a mapper for codim dim entities in the leaf grid
  using LeafGridView = Dune::GridView<Dune::DefaultLeafGridViewTraits<GridType>>;
  Dune::MultipleCodimMultipleGeomTypeMapper<LeafGridView>  mapper(gv.leafGridView(), Dune::mcmgVertexLayout());
  // iterate over vertices and find slaves
  LeafVertexIterator start = gv.leafGridView().template begin<dim>();
  for (LeafVertexIterator it = start; it != itend; ++it) {
    if (isOnPlane(dir,it->geometry().corner(0),coord)) {
      BoundaryGrid::Vertex v;
      v.i = mapper.index(*it);
      BoundaryGrid::extract(v.c,it->geometry().corner(0),log2(float(dir)));
      result.push_back(v);
    }
  }

  return result;
}


IMPL_FUNC(BoundaryGrid, extractMasterFace(Direction dir,
                                          ctype coord,
                                          SIDE side,
                                          bool dc))
{
  static const int V1[3][4] = {{0,2,4,6},
                               {0,1,4,5},
                               {0,1,2,3}};
  static const int V2[3][4] = {{1,3,5,7},
                               {2,3,6,7},
                               {4,5,6,7}};
  const LeafIndexSet& set = gv.leafGridView().indexSet();

  int c = 0;
  int i = log2(float(dir));
  BoundaryGrid result;
  // we first group nodes into this map through the coordinate of lower left 
  // vertex. we then split this up into pillars for easy processing later
  std::map<double, std::vector<BoundaryGrid::Quad> > nodeMap;
  for (LeafIterator cell  = gv.leafGridView().template begin<0>(); 
                    cell != gv.leafGridView().template end<0>(); ++cell, ++c) {
    std::vector<BoundaryGrid::Vertex> verts;
    int idx=0; 
    if (side == LEFT)
     idx = set.subIndex(*cell,V1[i][0],dim);
    else if (side == RIGHT)
     idx = set.subIndex(*cell,V2[i][0],dim);
    Dune::FieldVector<double, 3> pos = gv.vertexPosition(idx);
    if (isOnPlane(dir,pos,coord)) {
      for (int j=0;j<4;++j) {
        if (side == LEFT)
          idx = set.subIndex(*cell,V1[i][j],dim);
        if (side == RIGHT)
          idx = set.subIndex(*cell,V2[i][j],dim);
        pos = gv.vertexPosition(idx);
        if (!isOnPlane(dir,pos,coord))
          continue;
        BoundaryGrid::Vertex v;
        BoundaryGrid::extract(v,pos,i);
        v.i = idx;
        verts.push_back(v);
      }
    }
    if (verts.size() == 4) {
      BoundaryGrid::Quad q;
      q.v[0] = minXminY(verts);
      q.v[1] = maxXminY(verts);
      if (dc) {
        q.v[2] = minXmaxY(verts);
        q.v[3] = maxXmaxY(verts);
      } else {
        q.v[2] = maxXmaxY(verts);
        q.v[3] = minXmaxY(verts);
      }
      std::map<double, std::vector<BoundaryGrid::Quad> >::iterator it;
      for (it  = nodeMap.begin(); it != nodeMap.end(); ++it) {
        if (fabs(it->first-q.v[0].c[0]) < 1.e-7) {
          it->second.push_back(q);
          break;
        }
      }
      if (it == nodeMap.end())
        nodeMap[q.v[0].c[0]].push_back(q);

      result.add(q);
    }
  }

  int p=0;
  std::map<double, std::vector<BoundaryGrid::Quad> >::const_iterator it;
  for (it = nodeMap.begin(); it != nodeMap.end(); ++it, ++p) {
    for (size_t ii=0;ii<it->second.size();++ii)
      result.addToColumn(p,it->second[ii]);
  }

  return result;
}

IMPL_FUNC(void, determineSideFaces(const double* min, const double* max))
{
  master.push_back(extractMasterFace(X,min[0]));
  master.push_back(extractMasterFace(Y,min[1]));
  master.push_back(extractMasterFace(Z,min[2]));

  slave.push_back(extractFace(X,max[0]));
  slave.push_back(extractFace(Y,max[1]));
  slave.push_back(extractFace(Z,max[2]));
}

IMPL_FUNC(void, findBoundaries(double* min, double* max))
{
  max[0] = max[1] = max[2] = -1e5;
  min[0] = min[1] = min[2] = 1e5;
  const LeafVertexIterator itend = gv.leafGridView().template end<dim>();

  // iterate over vertices and find slaves
  LeafVertexIterator start = gv.leafGridView().template begin<dim>();
  for (LeafVertexIterator it = start; it != itend; ++it) {
    for (int i=0;i<3;++i) {
      min[i] = std::min(min[i],it->geometry().corner(0)[i]);
      max[i] = std::max(max[i],it->geometry().corner(0)[i]);
    }
  }
}

IMPL_FUNC(void, addMPC(Direction dir, int slavenode,
                       const BoundaryGrid::Vertex& m))
{
  MPC* mpc = new MPC(slavenode,log2(float(dir))+1);
  if (m.i > -1) { // we matched a node exactly
    mpc->addMaster(m.i,log2(float(dir))+1,1.f);
  } else {
    std::vector<double> N = m.q->evalBasis(m.c[0],m.c[1]);
    for (int i=0;i<4;++i)
      mpc->addMaster(m.q->v[i].i,log2(float(dir))+1,N[i]);
  }
  A.addMPC(mpc);
}

IMPL_FUNC(void, periodicPlane(Direction /*plane */,
                              Direction /* dir */,
                              const std::vector<BoundaryGrid::Vertex>& slavepointgrid,
                              const BoundaryGrid& mastergrid))
{
  for (size_t i=0;i<slavepointgrid.size();++i) {
    BoundaryGrid::Vertex coord;
    if (mastergrid.find(coord,slavepointgrid[i])) {
      addMPC(X,slavepointgrid[i].i,coord);
      addMPC(Y,slavepointgrid[i].i,coord);
      addMPC(Z,slavepointgrid[i].i,coord);
    }
  }
}

//! \brief Static helper to renumber a Q4 mesh to a P_1/P_N lag mesh.
//! \param[in] b The boundary grid describing the Q4 mesh
//! \param[in] n1 Number of DOFS in the first direction for each element
//! \param[in] n2 Number of DOFS in the first direction for each element
//! \param[out] totalDOFs The total number of free DOFs
static std::vector< std::vector<int> > renumber(const BoundaryGrid& b,
                                                int n1, int n2,
                                                int& totalDOFs)
{
  std::vector<std::vector<int> > nodes;
  nodes.resize(b.size());
  // loop over elements
  totalDOFs = 0;

  // fix lower left multiplicator.
  // will be "transfered" to all corners through periodic conditions
  nodes[0].push_back(-1);
  for (size_t e=0; e < b.size(); ++e) {
    // first direction major ordered nodes within each element
    for (int i2=0; i2 < n2; ++i2) {
      if (e != 0)
        nodes[e].push_back(nodes[e-1][i2*n1+n1-1]);

      int start = (e==0 && i2 != 0)?0:1;

      // slave the buggers
      if (i2 == n2-1 && n2 > 2) {
        for (int i1=(e==0?0:1); i1 < n1; ++i1) {
          nodes[e].push_back(nodes[e][i1]);
        }
      } else {
        for (int i1=start; i1 < n1; ++i1) {
          if (e == b.size()-1)
            nodes[e].push_back(nodes[0][i2*n1]);
          else
            nodes[e].push_back(totalDOFs++);
        }
      }
    }
  }

  return nodes;
}

IMPL_FUNC(int, addBBlockMortar(const BoundaryGrid& b1,
                               const BoundaryGrid& interface,
                               int /* dir */, int n1, int n2,
                               int colofs))
{
  // renumber the linear grid to the real multiplier grid
  int totalEqns;
  std::vector<std::vector<int> > lnodes = renumber(interface, n1+1,
                                                   n2+1, totalEqns);
  if (Bpatt.empty())
    Bpatt.resize(A.getEqns());

  // process pillar by pillar
  for (size_t p=0;p<interface.size();++p) {
    for (size_t q=0;q<b1.colSize(p);++q) {
      for (size_t i=0;i<4;++i) {
        for (size_t d=0;d<3;++d) {
          const MPC* mpc = A.getMPC(b1.getQuad(p,q).v[i].i,d);
          if (mpc) {
            for (size_t n=0;n<mpc->getNoMaster();++n) {
              int dof = A.getEquationForDof(mpc->getMaster(n).node,d);
              if (dof > -1) {
                for (size_t j=0;j<lnodes[p].size();++j) {
                  int indexj = 3*lnodes[p][j]+d;
                  if (indexj > -1)
                    Bpatt[dof].insert(indexj+colofs);
                }
              }
            }
          } else {
            int dof = A.getEquationForDof(b1.getQuad(p,q).v[i].i,d);
            if (dof > -1) {
              for (size_t j=0;j<lnodes[p].size();++j) {
                int indexj = 3*lnodes[p][j]+d;
                if (indexj > -1)
                  Bpatt[dof].insert(indexj+colofs);
              }
            }
          }
        }
      }
    }
  }

  return 3*totalEqns;
}

IMPL_FUNC(void, assembleBBlockMortar(const BoundaryGrid& b1,
                                     const BoundaryGrid& interface,
                                     int dir, int n1,
                                     int n2, int colofs,
                                     double alpha))
{
  // get a set of P1 shape functions for the displacements
  P1ShapeFunctionSet<ctype,ctype,2> ubasis = 
                P1ShapeFunctionSet<ctype,ctype,2>::instance();

  // get a set of PN shape functions for the multipliers
  PNShapeFunctionSet<2> lbasis(n1+1, n2+1);

  // renumber the linear grid to the real multiplier grid
  int totalEqns;
  std::vector<std::vector<int> > lnodes = renumber(interface, n1+1,
                                                   n2+1, totalEqns);
  // get a reference element
  auto gt = Dune::GeometryTypes::cube(2);
  // get a quadrature rule
  int quadorder = std::max((1.0+n1+0.5)/2.0,(1.0+n2+0.5)/2.0);
  quadorder = std::max(quadorder, 2);
  const Dune::QuadratureRule<ctype,2>& rule = 
                  Dune::QuadratureRules<ctype,2>::rule(gt,quadorder);

  // do the assembly loop
  typename Dune::QuadratureRule<ctype,2>::const_iterator r;
  Dune::DynamicMatrix<ctype> lE(ubasis.n,(n1+1)*(n2+1),0.0);
  LoggerHelper help(interface.size(), 5, 1000);
  for (int g=0;g<5;++g) {
    for (int p=help.group(g).first;p<help.group(g).second;++p) {
      const BoundaryGrid::Quad& qi(interface[p]);
      HexGeometry<2,2,GridType> lg(qi);
      for (size_t q=0;q<b1.colSize(p);++q) {
        const BoundaryGrid::Quad& qu = b1.getQuad(p,q);
        HexGeometry<2,2,GridType> hex(qu,gv,dir);
        lE = 0;
        for (r = rule.begin(); r != rule.end();++r) {
          ctype detJ = hex.integrationElement(r->position());
          if (detJ < 0)
            assert(0);

          typename HexGeometry<2,2,GridType>::LocalCoordinate loc =
                                          lg.local(hex.global(r->position()));
          assert(loc[0] <= 1.0+1.e-4 && loc[0] >= 0.0 && loc[1] <= 1.0+1.e-4 && loc[1] >= 0.0);
          for (int i=0;i<ubasis.n;++i) {
            for (int j=0;j<lbasis.size();++j) {
              lE[i][j] += ubasis[i].evaluateFunction(r->position())*
                          lbasis[j].evaluateFunction(loc)*detJ*r->weight();
            }
          }
        }

        // and assemble element contributions
        for (int d=0;d<3;++d) {
          for (int i=0;i<4;++i) {
            const MPC* mpc = A.getMPC(qu.v[i].i,d);
            if (mpc) {
              for (size_t n=0;n<mpc->getNoMaster();++n) {
                int indexi = A.getEquationForDof(mpc->getMaster(n).node,d);
                if (indexi > -1) {
                  for (size_t j=0;j<lnodes[p].size();++j) {
                    int indexj = lnodes[p][j]*3+d;
                    if (indexj > -1)
                      B[indexi][indexj+colofs] += alpha*lE[i][j];
                  }
                }
              }
            } else {
              int indexi = A.getEquationForDof(qu.v[i].i,d);
              if (indexi > -1) {
                for (size_t j=0;j<lnodes[p].size();++j) {
                  int indexj = lnodes[p][j]*3+d;
                  if (indexj > -1)
                    B[indexi][indexj+colofs] += alpha*lE[i][j];
                }
              }
            }
          }
        }
      }
    }
    help.log(g, "\t\t\t... still processing ... pillar ");
  }
}

IMPL_FUNC(void, fixPoint(Direction dir,
                         GlobalCoordinate coord,
                         const NodeValue& value))
{
  typedef typename GridType::LeafGridView::template Codim<dim>::Iterator VertexLeafIterator;
  const VertexLeafIterator itend = gv.leafGridView().template end<dim>();

  // make a mapper for codim 0 entities in the leaf grid 
  using LeafGridView = Dune::GridView<Dune::DefaultLeafGridViewTraits<GridType>>;
  Dune::MultipleCodimMultipleGeomTypeMapper<LeafGridView> mapper(gv.leafGridView(), Dune::mcmgVertexLayout());

  // iterate over vertices
  for (VertexLeafIterator it = gv.leafGridView().template begin<dim>(); it != itend; ++it) {
    if (isOnPoint(it->geometry().corner(0),coord)) {
      int indexi = mapper.index(*it);
      A.updateFixedNode(indexi,std::make_pair(dir,value));
    }
  }
}

IMPL_FUNC(bool, isOnPlane(Direction plane,
                          GlobalCoordinate coord,
                          ctype value))
{
  if (plane < X || plane > Z)
    return false;
  int p = log2(float(plane));
  ctype delta = fabs(value-coord[p]);
  return delta < tol;
}

IMPL_FUNC(void, fixLine(Direction dir,
                        ctype x, ctype y,
                        const NodeValue& value))
{
  typedef typename GridType::LeafGridView::template Codim<dim>::Iterator VertexLeafIterator;
  const VertexLeafIterator itend = gv.leafGridView().template end<dim>();

  // make a mapper for codim 0 entities in the leaf grid
  using LeafGridView = Dune::GridView<Dune::DefaultLeafGridViewTraits<GridType>>;
  Dune::MultipleCodimMultipleGeomTypeMapper<LeafGridView> mapper(gv, Dune::mcmgVertexLayout());

  // iterate over vertices
  for (VertexLeafIterator it = gv.leafGridView().template begin<dim>(); it != itend; ++it) {
    if (isOnLine(dir,it->geometry().corner(0),x,y)) {
      int indexi = mapper.index(*it);
      A.updateFixedNode(indexi,std::make_pair(XYZ,value));
    }
  }
}

IMPL_FUNC(bool, isOnLine(Direction dir,
                         GlobalCoordinate coord,
                         ctype x, ctype y))
{
  if (dir < X || dir > Z)
    return false;
  int ix = int(log2(float(dir))+1) % 3;
  int iy = int(log2(float(dir))+2) % 3;
  ctype delta = x-coord[ix];
  if (delta > tol || delta < -tol)
    return false;
  delta = y-coord[iy];
  if (delta > tol || delta < -tol)
    return false;

  return true;
}

IMPL_FUNC(bool, isOnPoint(GlobalCoordinate coord,
                          GlobalCoordinate point))
{
  GlobalCoordinate delta = point-coord;
  return delta.one_norm() < tol;
}

IMPL_FUNC(void, assemble(int loadcase, bool matrix))
{
  const int comp = 3+(dim-2)*3;
  static const int bfunc = 4+(dim-2)*4;

  Dune::FieldVector<ctype,comp> eps0;
  eps0 = 0;
  if (loadcase > -1) {
    eps0[loadcase] = 1;
    b[loadcase] = 0;
  }
  if (matrix)
    A.getOperator() = 0;

  for (int i=0;i<2;++i) {
    if (color[1].size() && matrix)
      std::cout << "\tprocessing " << (i==0?"red ":"black ") << "elements" << std::endl;
#pragma omp parallel for schedule(static)
    for (size_t j=0;j<color[i].size();++j) {
      Dune::FieldMatrix<ctype,comp,comp> C;
      Dune::FieldMatrix<ctype,dim*bfunc,dim*bfunc> K;
      Dune::FieldMatrix<ctype,dim*bfunc,dim*bfunc>* KP=0;
      Dune::FieldVector<ctype,dim*bfunc> ES;
      Dune::FieldVector<ctype,dim*bfunc>* EP=0;
      if (matrix)
        KP = &K;
      if (loadcase > -1)
        EP = &ES;

      for (size_t k=0;k<color[i][j].size();++k) {
        LeafIterator it = gv.leafGridView().template begin<0>();
        for (int l=0;l<color[i][j][k];++l)
          ++it;
        materials[color[i][j][k]]->getConstitutiveMatrix(C);
        // determine geometry type of the current element and get the matching reference element
        Dune::GeometryType gt = it->type();

        Dune::FieldMatrix<ctype,dim*bfunc,dim*bfunc> Aq;
        K = 0;
        ES = 0;

        // get a quadrature rule of order two for the given geometry type
        const Dune::QuadratureRule<ctype,dim>& rule = Dune::QuadratureRules<ctype,dim>::rule(gt,2);
        for (typename Dune::QuadratureRule<ctype,dim>::const_iterator r = rule.begin();
            r != rule.end() ; ++r) {
          // compute the jacobian inverse transposed to transform the gradients
          Dune::FieldMatrix<ctype,dim,dim> jacInvTra =
            it->geometry().jacobianInverseTransposed(r->position());

          ctype detJ = it->geometry().integrationElement(r->position());
          if (detJ <= 1.e-5 && verbose) {
            std::cout << "cell " << color[i][j][k] << " is (close to) degenerated, detJ " << detJ << std::endl;
            double zdiff=0.0;
            for (int ii=0;ii<4;++ii)
              zdiff = std::max(zdiff, it->geometry().corner(ii+4)[2]-it->geometry().corner(ii)[2]);
            std::cout << " - Consider setting ctol larger than " << zdiff << std::endl;
          }

          Dune::FieldMatrix<ctype,comp,dim*bfunc> lB;
          E.getBmatrix(lB,r->position(),jacInvTra);

          if (matrix) {
            E.getStiffnessMatrix(Aq,lB,C,detJ*r->weight());
            K += Aq;
          }

          // load vector
          if (EP) {
            Dune::FieldVector<ctype,dim*bfunc> temp;
            temp = Dune::FMatrixHelp::multTransposed(lB,Dune::FMatrixHelp::mult(C,eps0));
            temp *= -detJ*r->weight();
            ES += temp;
          }
        }
        A.addElement(KP,EP,it,(loadcase > -1)?&b[loadcase]:NULL);
      }
    }
  }
}

IMPL_FUNC(template<int comp> void,
          averageStress(Dune::FieldVector<ctype,comp>& sigma,
                        const Vector& uarg, int loadcase))
{
  if (loadcase < 0 || loadcase > 5)
    return;

  static const int bfunc = 4+(dim-2)*4;

  const LeafIterator itend = gv.leafGridView().template end<0>();

  Dune::FieldMatrix<ctype,comp,comp> C;
  Dune::FieldVector<ctype,comp> eps0;
  eps0 = 0;
  eps0[loadcase] = 1;
  int m=0;
  sigma = 0;
  double volume=0;
  for (LeafIterator it = gv.leafGridView().template begin<0>(); it != itend; ++it) {
    materials[m++]->getConstitutiveMatrix(C);
    // determine geometry type of the current element and get the matching reference element
    Dune::GeometryType gt = it->type();

    Dune::FieldVector<ctype,bfunc*dim> v;
    A.extractValues(v,uarg,it);

    // get a quadrature rule of order two for the given geometry type
    const Dune::QuadratureRule<ctype,dim>& rule = Dune::QuadratureRules<ctype,dim>::rule(gt,2);
    for (typename Dune::QuadratureRule<ctype,dim>::const_iterator r = rule.begin();
        r != rule.end() ; ++r) {
      // compute the jacobian inverse transposed to transform the gradients
      Dune::FieldMatrix<ctype,dim,dim> jacInvTra =
        it->geometry().jacobianInverseTransposed(r->position());

      ctype detJ = it->geometry().integrationElement(r->position());

      volume += detJ*r->weight();

      Dune::FieldMatrix<ctype,comp,dim*bfunc> lB;
      E.getBmatrix(lB,r->position(),jacInvTra);

      Dune::FieldVector<ctype,comp> s;
      E.getStressVector(s,v,eps0,lB,C);
      s *= detJ*r->weight();
      sigma += s;
    }
  }
  sigma /= volume;
  if (Escale > 0)
    sigma /= Escale/Emin;
}

IMPL_FUNC(void, loadMaterialsFromGrid(const std::string& file))
{
  typedef std::map<std::pair<double,double>,
                   std::shared_ptr<Material> > MaterialMap;
  MaterialMap cache;
  std::vector<double> Emod;
  std::vector<double> Poiss;
  std::vector<int> satnum;
  std::vector<double> rho;
  upscaledRho = -1;

  if (file == "uniform") {
    int cells = gv.size(0);
    Emod.insert(Emod.begin(),cells,100.f);
    Poiss.insert(Poiss.begin(),cells,0.38f);
  } else {
    Opm::Parser parser;
    auto deck = parser.parseFile(file);
    if (deck.hasKeyword("YOUNGMOD") && deck.hasKeyword("POISSONMOD")) {
        Emod = deck["YOUNGMOD"].back().getRawDoubleData();
        Poiss = deck["POISSONMOD"].back().getRawDoubleData();
      std::vector<double>::const_iterator it = std::min_element(Poiss.begin(), Poiss.end());
      if (*it < 0) {
        std::cerr << "Auxetic material specified for cell " << it-Poiss.begin() << std::endl
                  << "Emod: "<< Emod[it-Poiss.begin()] << " Poisson's ratio: " << *it << std::endl << "bailing..." << std::endl;
        exit(1);
      }
    } else if (deck.hasKeyword("LAMEMOD") && deck.hasKeyword("SHEARMOD")) {
        std::vector<double> lame = deck["LAMEMOD"].back().getRawDoubleData();
        std::vector<double> shear = deck["SHEARMOD"].back().getRawDoubleData();
      Emod.resize(lame.size());
      Poiss.resize(lame.size());
      for (size_t i=0;i<lame.size();++i) {
        Emod[i]  = shear[i]*(3*lame[i]+2*shear[i])/(lame[i]+shear[i]);
        Poiss[i] = 0.5*lame[i]/(lame[i]+shear[i]);
      }
      std::vector<double>::const_iterator it = std::min_element(Poiss.begin(), Poiss.end());
      if (*it < 0) {
        std::cerr << "Auxetic material specified for cell " << it-Poiss.begin() << std::endl
                  << "Lamè modulus: " << lame[it-Poiss.begin()] << " Shearmodulus: " << shear[it-Poiss.begin()] << std::endl
                  << "Emod: "<< Emod[it-Poiss.begin()] << " Poisson's ratio: " << *it << std::endl << "bailing..." << std::endl;
        exit(1);
      }
    } else if (deck.hasKeyword("BULKMOD") && deck.hasKeyword("SHEARMOD")) {
        std::vector<double> bulk = deck["BULKMOD"].back().getRawDoubleData();
        std::vector<double> shear = deck["SHEARMOD"].back().getRawDoubleData();
      Emod.resize(bulk.size());
      Poiss.resize(bulk.size());
      for (size_t i=0;i<bulk.size();++i) {
        Emod[i]  = 9*bulk[i]*shear[i]/(3*bulk[i]+shear[i]);
        Poiss[i] = 0.5*(3*bulk[i]-2*shear[i])/(3*bulk[i]+shear[i]);
      }
      std::vector<double>::const_iterator it = std::min_element(Poiss.begin(), Poiss.end());
      if (*it < 0) {
        std::cerr << "Auxetic material specified for cell " << it-Poiss.begin() << std::endl
                  << "Bulkmodulus: " << bulk[it-Poiss.begin()] << " Shearmodulus: " << shear[it-Poiss.begin()] << std::endl
                 << "Emod: "<< Emod[it-Poiss.begin()] << " Poisson's ratio: " << *it << std::endl << "bailing..." << std::endl;
        exit(1);
      }
    } else if (deck.hasKeyword("PERMX") && deck.hasKeyword("PORO")) {
      std::cerr << "WARNING: Using PERMX and PORO for elastic material properties" << std::endl;
      Emod = deck["PERMX"].back().getRawDoubleData();
      Poiss = deck["PORO"].back().getRawDoubleData();
    } else {
      std::cerr << "No material data found in eclipse file, aborting" << std::endl;
      exit(1);
    }
    if (deck.hasKeyword("SATNUM"))
        satnum = deck["SATNUM"].back().getIntData();
    if (deck.hasKeyword("RHO"))
        rho = deck["RHO"].back().getRawDoubleData();
  }
  // scale E modulus of materials
  if (Escale > 0) {
    Emin = *std::min_element(Emod.begin(),Emod.end());
    for (size_t i=0;i<Emod.size();++i)
      Emod[i] *= Escale/Emin;
  }

  // make sure we only instance a minimal amount of materials.
  // also map the correct material to the correct cells.
  // their original ordering is as in the eclipse file itself
  // while globalCell holds the map of cells kept after preprocessing
  // the grid
  std::vector<int> cells = gv.globalCell();
  int j=0;
  std::map<Material*,double> volume;
  for (size_t i=0;i<cells.size();++i) {
    int k = cells[i];
    MaterialMap::iterator it;
    if ((it = cache.find(std::make_pair(Emod[k],Poiss[k]))) != cache.end()) {
      assert(gv.cellVolume(i) > 0);
      volume[it->second.get()] += gv.cellVolume(i);
      materials.push_back(it->second);
    } else {
      std::shared_ptr<Material> mat(new Isotropic(j++,Emod[k],Poiss[k]));
      cache.insert(std::make_pair(std::make_pair(Emod[k],Poiss[k]),mat));
      assert(gv.cellVolume(i) > 0);
      volume.insert(std::make_pair(mat.get(),gv.cellVolume(i)));
      materials.push_back(mat);
    }
    if (!satnum.empty()) {
      if (satnum[k] > (int)volumeFractions.size())
        volumeFractions.resize(satnum[k]);
      volumeFractions[satnum[k]-1] += gv.cellVolume(i);
    }
    if (!rho.empty())
      upscaledRho += gv.cellVolume(i)*rho[k];
  }
  std::cout << "Number of materials: " << cache.size() << std::endl;

  double totalvolume=0;
  for (std::map<Material*,double>::iterator it  = volume.begin(); 
                                            it != volume.end(); ++it) 
    totalvolume += it->second;

  // statistics
  if (satnum.empty()) {
    int i=0;
    for (MaterialMap::iterator it = cache.begin(); it != cache.end(); ++it, ++i) {
      std::cout << "  Material" << i+1 << ": " << 100.f*volume[it->second.get()]/totalvolume << '%' << std::endl;
      volumeFractions.push_back(volume[it->second.get()]/totalvolume);
    }
    bySat = false;
  } else {
    for (size_t jj=0; jj < volumeFractions.size(); ++jj) {
      volumeFractions[jj] /= totalvolume;
      std::cout << "SATNUM " << jj+1 << ": " << volumeFractions[jj]*100 << '%' << std::endl;
    }
    bySat = true;
  }
  if (upscaledRho > 0) {
    upscaledRho /= totalvolume;
    std::cout << "Upscaled density: " << upscaledRho << std::endl;
  }
}

IMPL_FUNC(void, loadMaterialsFromRocklist(const std::string& file,
                                          const std::string& rocklist))
{
  std::vector< std::shared_ptr<Material> > cache;
  // parse the rocklist
  std::ifstream f;
  f.open(rocklist.c_str());
  int mats;
  f >> mats;
  for (int i=0;i<mats;++i) {
    std::string matfile;
    f >> matfile;
    cache.push_back(std::shared_ptr<Material>(Material::create(i+1,matfile)));
  }

  // scale E modulus of materials
  if (Escale > 0) {
    Emin=1e10;
    for (size_t i=0;i<cache.size();++i)
      Emin = std::min(Emin,((Isotropic*)cache[i].get())->getE());
    for (size_t i=0;i<cache.size();++i) {
      double lE = ((Isotropic*)cache[i].get())->getE();
      ((Isotropic*)cache[i].get())->setE(lE*Escale/Emin);
    }
  }
  std::vector<double> volume;
  volume.resize(cache.size());
  if (file == "uniform") {
    if (cache.size() == 1) {
      for (int i=0;i<gv.size(0);++i)
        materials.push_back(cache[0]);
      volume[0] = 1;
    } else {
      for (int i=0;i<gv.size(0);++i) {
        materials.push_back(cache[i % cache.size()]);
        volume[i % cache.size()] += gv.cellVolume(i);
      }
    }
  } else {
    Opm::Parser parser;
    auto deck = parser.parseFile(file);
    std::vector<int> satnum = deck["SATNUM"].back().getIntData();
    std::vector<int> cells = gv.globalCell();
    for (size_t i=0;i<cells.size();++i) {
      int k = cells[i];
      if (satnum[k]-1 >= (int)cache.size()) {
        std::cerr << "Material " << satnum[k] << " referenced but not available. Check your rocklist." << std::endl;
        exit(1);
      }
      materials.push_back(cache[satnum[k]-1]);
      volume[satnum[k]-1] += gv.cellVolume(i);
    }
  }
  std::cout << "Number of materials: " << cache.size() << std::endl;
  // statistics
  double totalvolume = std::accumulate(volume.begin(),volume.end(),0.f);
  for (size_t i=0;i<cache.size();++i) {
    std::cout << " SATNUM " << i+1 << ": " << 100.f*volume[i]/totalvolume << '%' << std::endl;
    volumeFractions.push_back(volume[i]/totalvolume);
  }
  bySat = true;
}

IMPL_FUNC(void, fixCorners(const double* min, const double* max))
{
  ctype c[8][3] = {{min[0],min[1],min[2]},
                   {max[0],min[1],min[2]},
                   {min[0],max[1],min[2]},
                   {max[0],max[1],min[2]},
                   {min[0],min[1],max[2]},
                   {max[0],min[1],max[2]},
                   {min[0],max[1],max[2]},
                   {max[0],max[1],max[2]}};
  for (int i=0;i<8;++i) {
    GlobalCoordinate coord;
    coord[0] = c[i][0]; coord[1] = c[i][1]; coord[2] = c[i][2];
    fixPoint(XYZ,coord);
  }
}

IMPL_FUNC(void, periodicBCs(const double* min, const double* max))
{
  // this method
  // 1. fixes the primal corner dofs
  // 2. extracts establishes a quad grid for the left and front sides,
  //    while a point grid is created for the right and back sides.
  // 3. establishes strong couplings (MPC)

  // step 1
  fixCorners(min,max);

  // step 2
  determineSideFaces(min,max);
  std::cout << "Xslave " << slave[0].size() << " "
            << "Yslave " << slave[1].size() << " "
            << "Zslave " << slave[2].size() << std::endl;
  std::cout << "Xmaster " << master[0].size() << " "
            << "Ymaster " << master[1].size() << " "
            << "Zmaster " << master[2].size() << std::endl;

  // step 3
  periodicPlane(X,XYZ,slave[0],master[0]);
  periodicPlane(Y,XYZ,slave[1],master[1]);
  periodicPlane(Z,XYZ,slave[2],master[2]);
}

IMPL_FUNC(void, periodicBCsMortar(const double* min, 
                                  const double* max,
                                  int n1, int n2,
                                  int p1, int p2))
{
  // this method
  // 1. fixes the primal corner dofs
  // 2. establishes strong couplings (MPC) on top and bottom
  // 3. extracts and establishes a quad grid for the left/right/front/back sides
  // 4. establishes grids for the dual dofs
  // 5. calculates the coupling matrix between the left/right sides
  // 6. calculates the coupling matrix between the front/back sides
  //
  // The Mortar linear system is of the form
  // [A   B] [u]   [b]
  // [B'  0] [l] = [0]

  // step 1
  fixCorners(min,max);
  
  // step 2
  std::cout << "\textracting nodes on top face..." << std::endl;
  slave.push_back(extractFace(Z,max[2]));
  std::cout << "\treconstructing bottom face..." << std::endl;
  BoundaryGrid bottom = extractMasterFace(Z,min[2]);
  std::cout << "\testablishing couplings on top/bottom..." << std::endl;
  periodicPlane(Z,XYZ,slave[0],bottom);
  std::cout << "\tinitializing matrix..." << std::endl;
  A.initForAssembly();

  // step 3
  std::cout << "\treconstructing left face..." << std::endl;
  master.push_back(extractMasterFace(X, min[0], LEFT, true));
  std::cout << "\treconstructing right face..." << std::endl;
  master.push_back(extractMasterFace(X, max[0], RIGHT, true));
  std::cout << "\treconstructing front face..." << std::endl;
  master.push_back(extractMasterFace(Y, min[1], LEFT, true));
  std::cout << "\treconstructing back face..." << std::endl;
  master.push_back(extractMasterFace(Y, max[1], RIGHT, true));

  std::cout << "\testablished YZ multiplier grid with " << n2 << "x1" << " elements" << std::endl;

  // step 4
  BoundaryGrid::FaceCoord fmin,fmax;
  fmin[0] = min[1]; fmin[1] = min[2];
  fmax[0] = max[1]; fmax[1] = max[2];
  BoundaryGrid lambdax = BoundaryGrid::uniform(fmin, fmax, n2, 1, true);

  fmin[0] = min[0]; fmin[1] = min[2];
  fmax[0] = max[0]; fmax[1] = max[2];
  std::cout << "\testablished XZ multiplier grid with " << n1 << "x1" << " elements" << std::endl;
  BoundaryGrid lambday = BoundaryGrid::uniform(fmin, fmax, n1, 1, true);

  addBBlockMortar(master[0], lambdax, 0, 1, p2, 0);
  int eqns = addBBlockMortar(master[1], lambdax, 0, 1, p2, 0);
  addBBlockMortar(master[2], lambday, 1, 1, p1, eqns);
  int eqns2 = addBBlockMortar(master[3], lambday, 1, 1, p1, eqns);

  MatrixOps::fromAdjacency(B, Bpatt, A.getEqns(), eqns+eqns2);
  Bpatt.clear();

  // step 5
  std::cout << "\tassembling YZ mortar matrix..." << std::endl;
  assembleBBlockMortar(master[0], lambdax, 0, 1, p2, 0);
  assembleBBlockMortar(master[1], lambdax, 0, 1, p2, 0, -1.0);

  // step 6
  std::cout << "\tassembling XZ mortar matrix..." << std::endl;
  assembleBBlockMortar(master[2], lambday, 1, 1, p1, eqns);
  assembleBBlockMortar(master[3], lambday, 1, 1, p1, eqns, -1.0);

  master.clear();
  slave.clear();
}

 template<class M, class A>
static void applyMortarBlock(int i, const Matrix& B, M& T,
                             A& upre)
{
  Vector v, v2;
  v.resize(B.N());
  v2.resize(B.N());
  v = 0;
  v2 = 0;
  MortarBlockEvaluator<Dune::Preconditioner<Vector,Vector> > pre(upre, B);
  v[i] = 1;
  pre.apply(v, v2);
  for (size_t j=0; j < B.M(); ++j)
    T[j][i] = v2[j];
}

IMPL_FUNC(void, setupSolvers(const LinSolParams& params))
{
  int siz = A.getOperator().N(); // system size
  int numsolvers = 1;
#ifdef HAVE_OPENMP
   numsolvers = omp_get_max_threads();
#endif
          
  if (params.type == ITERATIVE) {
    op.reset(new Operator(A.getOperator()));
    bool copy;
    upre.push_back(PC::setup(params.steps[0], params.steps[1],
                             params.coarsen_target, params.zcells,
                             op, gv, A, copy));

    // Mortar in use
    if (B.N()) {
      siz += B.M();

      // schur system: B'*P^-1*B
      if (params.mortarpre >= SCHUR) {
        Dune::DynamicMatrix<double> T(B.M(), B.M());
        std::cout << "\tBuilding preconditioner for multipliers..." << std::endl;
        LoggerHelper help(B.M(), 5, 100);

        std::vector< std::shared_ptr<Dune::Preconditioner<Vector,Vector> > > pc;
        pc.resize(B.M());

        if (params.mortarpre == SCHUR ||
            (params.mortarpre == SCHURAMG &&
             params.pre == AMG) ||
            (params.mortarpre == SCHURSCHWARZ &&
             params.pre == SCHWARZ) ||
            (params.mortarpre == SCHURTWOLEVEL &&
             params.pre == TWOLEVEL)) {
          pc[0] = upre[0];
          if (copy) {
            for (size_t i=1;i<B.M();++i)
              pc[i].reset(new PCType(*upre[0]));
          }
        } else if (params.mortarpre == SCHURAMG) {
          std::shared_ptr<AMG1<SSORSmoother>::type> mpc;
          pc[0] = mpc = AMG1<SSORSmoother>::setup(params.steps[0],
                                                  params.steps[1],
                                                  params.coarsen_target,
                                                  params.zcells,
                                                  op, gv, A, copy);
          for (size_t i=1;i<B.M();++i)
            pc[i].reset(new AMG1<SSORSmoother>::type(*mpc));
        } else if (params.mortarpre == SCHURSCHWARZ) {
          std::shared_ptr<Schwarz::type> mpc;
          pc[0] = mpc = Schwarz::setup(params.steps[0],
                                       params.steps[1],
                                       params.coarsen_target,
                                       params.zcells,
                                       op, gv, A, copy);
        } else if (params.mortarpre == SCHURTWOLEVEL) {
          std::shared_ptr<AMG2Level<SSORSmoother>::type> mpc;
          pc[0] = mpc = AMG2Level<SSORSmoother>::setup(params.steps[0],
                                                       params.steps[1],
                                                       params.coarsen_target,
                                                       params.zcells,
                                                       op, gv, A, copy);
          for (size_t i=1;i<B.M();++i)
            pc[i].reset(new AMG2Level<SSORSmoother>::type(*mpc));
        }
        for (int t=0;t<5;++t) {
#pragma omp parallel for schedule(static)
          for (int i=help.group(t).first; i < help.group(t).second; ++i)
            applyMortarBlock(i, B, T, *pc[copy?i:0]);

          help.log(help.group(t).second,
                   "\t\t... still processing ... multiplier ");
        }
        P = MatrixOps::fromDense(T);
      } else if (params.mortarpre == SIMPLE) {
        Matrix D = MatrixOps::diagonal(A.getEqns());

        // scale by row sums
        size_t row=0;
        for (Matrix::ConstRowIterator it  = A.getOperator().begin();
                                      it != A.getOperator().end(); ++it, ++row) {
          double alpha=0;
          for (Matrix::ConstColIterator it2  = it->begin(); 
                                        it2 != it->end(); ++it2) {
            if (it2.index() != row)
              alpha += fabs(*it2);
          }
          D[row][row] = 1.0/(A.getOperator()[row][row]/alpha);
        }

        Matrix t1;
        // t1 = Ad*B
        Dune::matMultMat(t1, D, B);
        // P = B'*t1 = B'*Ad*B
        Dune::transposeMatMultMat(P,  B, t1);
      }

      if (params.uzawa) {
        std::shared_ptr<Dune::InverseOperator<Vector,Vector> > innersolver;

        innersolver.reset(new Dune::CGSolver<Vector>(*op, *upre[0], params.tol,
                                                     params.maxit,
                                                     verbose?2:(params.report?1:0)));
        op2.reset(new SchurEvaluator(*innersolver, B));
        lprep.reset(new LUSolver(P));
        lpre.reset(new SeqLU(*lprep));
        std::shared_ptr<Dune::InverseOperator<Vector,Vector> > outersolver;
        outersolver.reset(new Dune::CGSolver<Vector>(*op2, *lpre, params.tol*10,
                                                     params.maxit,
                                                     verbose?2:(params.report?1:0)));
        tsolver.push_back(SolverPtr(new UzawaSolver<Vector, Vector>(innersolver, outersolver, B)));
      } else {
        for (int i=0;i<numsolvers;++i) {
          if (copy && i != 0)
            upre.push_back(PCPtr(new PCType(*upre[0])));
          tmpre.push_back(MortarAmgPtr(new MortarSchurPre<PCType>(P, B,
                                                                *upre[copy?i:0],
                                                            params.symmetric)));
        }
        meval.reset(new MortarEvaluator(A.getOperator(), B));
        if (params.symmetric) {
          for (int i=0;i<numsolvers;++i)
           tsolver.push_back(SolverPtr(new Dune::MINRESSolver<Vector>(*meval, *tmpre[i],
                                                                      params.tol,
                                                                      params.maxit,
                                                                      verbose?2:(params.report?1:0))));
        } else {
          for (int i=0;i<numsolvers;++i)
            tsolver.push_back(SolverPtr(new Dune::RestartedGMResSolver<Vector>(*meval, *tmpre[i],
                                                                               params.tol,
                                                                               params.restart,
                                                                               params.maxit,
                                                                               verbose?2:(params.report?1:0))));
        }
      }
    } else {
      for (int i=0;i<numsolvers;++i) {
        if (copy && i != 0)
          upre.push_back(PCPtr(new PCType(*upre[0])));
        tsolver.push_back(SolverPtr(new Dune::CGSolver<Vector>(*op, *upre[copy?i:0],
                                                               params.tol,
                                                               params.maxit,
                                                               verbose?2:(params.report?1:0))));
      }
    }
  } else {
    if (B.N()) 
      A.getOperator() = MatrixOps::augment(A.getOperator(), B,
                                           0, A.getOperator().M(), true);
#if HAVE_UMFPACK || HAVE_SUPERLU
    tsolver.push_back(SolverPtr(new LUSolver(A.getOperator(),
                                             verbose?2:(params.report?1:0))));
    for (int i=1;i<numsolvers;++i)
      tsolver.push_back(tsolver.front());
#else
    std::cerr << "No direct solver available" << std::endl;
    exit(1);
#endif
    siz = A.getOperator().N();
  }

  for (int i=0;i<6;++i)
    b[i].resize(siz);
}

IMPL_FUNC(void, solve(int loadcase))
{
  try {
    Dune::InverseOperatorResult r;
#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 7)
    u[loadcase].resize(b[loadcase].size());
#else
    u[loadcase].resize(b[loadcase].size(), false);
#endif
    u[loadcase] = 0;
    int solver=0;
#ifdef HAVE_OPENMP
    solver = omp_get_thread_num();
#endif

    tsolver[solver]->apply(u[loadcase], b[loadcase], r);

    std::cout << "\tsolution norm: " << u[loadcase].two_norm() << std::endl;
  } catch (Dune::ISTLError& e) {
    std::cerr << "exception thrown " << e << std::endl;
  }
}

}} // namespace Opm, Elasticity

#endif
