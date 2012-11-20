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

  template<class GridType>
std::vector<BoundaryGrid::Vertex> ElasticityUpscale<GridType>::extractFace(Direction dir, ctype coord)
{
  std::vector<BoundaryGrid::Vertex> result;
  const LeafVertexIterator itend = gv.leafView().template end<dim>();

  // make a mapper for codim dim entities in the leaf grid 
  Dune::LeafMultipleCodimMultipleGeomTypeMapper<GridType,
                                            Dune::MCMGVertexLayout> mapper(gv);
  // iterate over vertices and find slaves
  LeafVertexIterator start = gv.leafView().template begin<dim>();
  for (LeafVertexIterator it = start; it != itend; ++it) {
    if (isOnPlane(dir,it->geometry().corner(0),coord)) {
      BoundaryGrid::Vertex v;
      v.i = mapper.map(*it);
      BoundaryGrid::extract(v.c,it->geometry().corner(0),log2(dir));
      result.push_back(v);
    }
  }

  return result;
}

  template<class GridType>
BoundaryGrid ElasticityUpscale<GridType>::extractMasterFace(Direction dir,
                                                            ctype coord,
                                                            SIDE side, bool dc)
{
  static const int V1[3][4] = {{0,2,4,6},
                               {0,1,4,5},
                               {0,1,2,3}};
  static const int V2[3][4] = {{1,3,5,7},
                               {2,3,6,7},
                               {4,5,6,7}};
  const LeafIndexSet& set = gv.leafView().indexSet();
  const LeafVertexIterator itend = gv.leafView().template end<dim>();

  // make a mapper for codim dim entities in the leaf grid 
  Dune::LeafMultipleCodimMultipleGeomTypeMapper<GridType,
                                            Dune::MCMGVertexLayout> mapper(gv);
  LeafVertexIterator start=gv.leafView().template begin<dim>();
  LeafIterator cellend = gv.leafView().template end<0>();
  int c = 0;
  int i = log2(dir);
  BoundaryGrid result;
  // we first group nodes into this map through the coordinate of lower left 
  // vertex. we then split this up into pillars for easy processing later
  std::map<double, std::vector<BoundaryGrid::Quad> > nodeMap;
  for (LeafIterator cell  = gv.leafView().template begin<0>(); 
                    cell != cellend; ++cell, ++c) {
    std::vector<BoundaryGrid::Vertex> verts;
    int idx; 
    if (side == LEFT)
     idx = set.subIndex(*cell,V1[i][0],dim);
    else if (side == RIGHT)
     idx = set.subIndex(*cell,V2[i][0],dim);
    LeafVertexIterator it=start;
    for (; it != itend; ++it) {
      if (mapper.map(*it) == idx)
        break;
    }
    if (isOnPlane(dir,it->geometry().corner(0),coord)) {
      for (int j=0;j<4;++j) {
        if (side == LEFT)
          idx = set.subIndex(*cell,V1[i][j],dim);
        if (side == RIGHT)
          idx = set.subIndex(*cell,V2[i][j],dim);
        LeafVertexIterator it=start;
        for (; it != itend; ++it) {
          if (mapper.map(*it) == idx)
            break;
        }
        if (!isOnPlane(dir,it->geometry().corner(0),coord))
          continue;
        BoundaryGrid::Vertex v;
        BoundaryGrid::extract(v,it->geometry().corner(0),i);
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
    for (size_t i=0;i<it->second.size();++i)
      result.addToColumn(p,it->second[i]);
  }

  return result;
}

  template<class GridType>
void ElasticityUpscale<GridType>::determineSideFaces(const double* min, 
                                                     const double* max)
{
  master.push_back(extractMasterFace(X,min[0]));
  master.push_back(extractMasterFace(Y,min[1]));
  master.push_back(extractMasterFace(Z,min[2]));

  slave.push_back(extractFace(X,max[0]));
  slave.push_back(extractFace(Y,max[1]));
  slave.push_back(extractFace(Z,max[2]));
}

  template<class GridType>
void ElasticityUpscale<GridType>::findBoundaries(double* min, 
                                                 double* max)
{
  max[0] = max[1] = max[2] = -1e5;
  min[0] = min[1] = min[2] = 1e5;
  const LeafVertexIterator itend = gv.leafView().template end<dim>();

  // iterate over vertices and find slaves
  LeafVertexIterator start = gv.leafView().template begin<dim>();
  for (LeafVertexIterator it = start; it != itend; ++it) {
    for (int i=0;i<3;++i) {
      min[i] = std::min(min[i],it->geometry().corner(0)[i]);
      max[i] = std::max(max[i],it->geometry().corner(0)[i]);
    }
  }
}


  template<class GridType>
void ElasticityUpscale<GridType>::addMPC(Direction dir, int slave,
                                         const BoundaryGrid::Vertex& m)
{
  MPC* mpc = new MPC(slave,log2(dir)+1);
  if (m.i > -1) { // we matched a node exactly
    mpc->addMaster(m.i,log2(dir)+1,1.f);
  } else {
    std::vector<double> N = m.q->evalBasis(m.c[0],m.c[1]);
    for (int i=0;i<4;++i)
      mpc->addMaster(m.q->v[i].i,log2(dir)+1,N[i]);
  }
  A.addMPC(mpc);
}

  template<class GridType>
void ElasticityUpscale<GridType>::periodicPlane(Direction plane, Direction dir, 
                                 const std::vector<BoundaryGrid::Vertex>& slave,
                                 const BoundaryGrid& master)
{
  for (size_t i=0;i<slave.size();++i) {
    BoundaryGrid::Vertex coord;
    if (master.find(coord,slave[i])) {
      addMPC(X,slave[i].i,coord);
      addMPC(Y,slave[i].i,coord);
      addMPC(Z,slave[i].i,coord);
    }
  }
}

  template<class GridType>
Matrix ElasticityUpscale<GridType>::findBMatrixLLM(const SlaveGrid& slave)
{
  std::vector< std::set<int> > adj;
  adj.resize(A.getEqns());
  int cols=0;
  std::vector<int> colmap;
  cols=0;
  for (std::map<int,BoundaryGrid::Vertex>::const_iterator it  = slave.begin();
                                                    it != slave.end();++it) {
    for (int k=0;k<dim;++k) {
      if (A.isFixed(it->first))
        continue;
      MPC* mpc = A.getMPC(it->first,k);
      if (mpc) {
        for (size_t n=0;n<mpc->getNoMaster();++n) {
          int idx = A.getEquationForDof(mpc->getMaster(n).node,
                                        mpc->getMaster(n).dof-1);
          if (idx != -1)
            adj[idx].insert(cols++);
        }
      } else {
        int idx = A.getEquationForDof(it->first,k);
        if (idx != -1)
          adj[idx].insert(cols++);
      }
    }
  }
  Matrix B;
  MatrixOps::fromAdjacency(B,adj,A.getEqns(),cols);
  int col=0;
  for (std::map<int,BoundaryGrid::Vertex>::const_iterator it  = slave.begin();
                                                    it != slave.end();++it) {
    for (int k=0;k<dim;++k) {
      if (A.isFixed(it->first))
        continue;
      MPC* mpc = A.getMPC(it->first,k);
      if (mpc) {
        for (size_t n=0;n<mpc->getNoMaster();++n) {
          int idx = A.getEquationForDof(mpc->getMaster(n).node,
                                        mpc->getMaster(n).dof-1);
          if (idx != -1)
            B[idx][col++] += 1;
        }
      } else {
        int idx = A.getEquationForDof(it->first,k);
        if (idx != -1)
          B[idx][col++] += 1;
      }
    }
  }

  return B;
}

  template<class GridType>
Matrix ElasticityUpscale<GridType>::findLMatrixLLM(const SlaveGrid& slave,
                                                   const BoundaryGrid& master)
{
  int nbeqn=0;
  std::vector<int> dofmap(master.totalNodes()*dim,-1);
  int col=0;
  for (size_t i=0;i<master.totalNodes();++i) {
    if (master.isFixed(i)) {
        dofmap[i*dim  ] = -1;
        dofmap[i*dim+1] = -1;
        dofmap[i*dim+2] = -1;
    }
    else {
      dofmap[i*dim  ] = col++;
      dofmap[i*dim+1] = col++;
      dofmap[i*dim+2] = col++;
    }
  }
  for (SlaveGrid::const_iterator it  = slave.begin();
                                 it != slave.end(); ++it) {
    if (!A.isFixed(it->first))
      nbeqn += 3;
  }
  std::vector< std::set<int> > adj;
  adj.resize(nbeqn);

  int row=0;
  for (SlaveGrid::const_iterator it  = slave.begin();
                                 it != slave.end();++it) {
    if (A.isFixed(it->first))
      continue;
    for (int k=0;k<dim;++k) {
      if (it->second.i == -1) {
        for (int i=0;i<4;++i) {
          int idx = dofmap[it->second.q->v[i].i*dim+k];
          if (idx > -1)
            adj[row].insert(idx);
        }
      }
      else {
        int idx = dofmap[it->second.i*dim+k];
        if (idx > -1)
          adj[row].insert(idx);
      }
      row++;
    }
  }
  Matrix L;
  MatrixOps::fromAdjacency(L,adj,nbeqn,col);
  row = 0;
  for (SlaveGrid::const_iterator it  = slave.begin();
                                 it != slave.end();++it) {
    if (A.isFixed(it->first))
      continue;
    for (int k=0;k<dim;++k) {
      if (it->second.i == -1) {
        std::vector<double> N = it->second.q->evalBasis(it->second.c[0],
                                                        it->second.c[1]);
        for (int i=0;i<4;++i) {
          int idx = dofmap[it->second.q->v[i].i*dim+k];
          if (idx > -1)
            L[row][idx] -= N[i];
        }
      }
      else {
        int idx = dofmap[it->second.i*dim+k];
        if (idx > -1)
          L[row][idx] -= 1;
      }
      row++;
    }
  }

  return L;
}

static std::vector< std::vector<int> > renumber(const BoundaryGrid& b,
                                                int n1, int n2)
{
  std::vector<std::vector<int> > nodes;
  nodes.resize(b.size());
  // loop over elements
  int ofs = 0;
  for (size_t e=0; e < b.size(); ++e) {
    // first direction major ordered nodes within each element
    for (int i2=0; i2 < n2; ++i2) {
      if (e != 0)
        nodes[e].push_back(nodes[e-1][i2*n1+n1-1]);
      for (int i1=(e==0?0:1); i1 < n1; ++i1)
        nodes[e].push_back(ofs++);
    }
  }

  return nodes;
}

  template<class GridType>
Matrix ElasticityUpscale<GridType>::findLMatrixMortar(const BoundaryGrid& b1,
                                                      const BoundaryGrid& interface,
                                                      int dir, int n1, int n2)
{
  std::vector< std::set<int> > adj;
  adj.resize(A.getEqns());

#define QINDEX(p,q,dir) (q)

  // get a set of P1 shape functions for the displacements
  P1ShapeFunctionSet<ctype,ctype,2> ubasis = 
                P1ShapeFunctionSet<ctype,ctype,2>::instance();

  // get a set of PN shape functions for the multipliers
  PNShapeFunctionSet<2> lbasis(n1+1, n2+1);

  std::vector<std::vector<int> > lnodes = renumber(interface,n1,n2);
  int totalNodes = lnodes.back().back()+1;

  // process pillar by pillar
  for (size_t p=0;p<interface.size();++p) {
    for (size_t q=0;q<b1.colSize(p);++q) {
      for (size_t i=0;i<ubasis.n;++i) {
        for (size_t d=0;d<3;++d) {
          MPC* mpc = A.getMPC(b1.getQuad(p,QINDEX(p,q,dir)).v[i].i,d);
          if (mpc) {
            for (size_t n=0;n<mpc->getNoMaster();++n) {
              int dof = A.getEquationForDof(mpc->getMaster(n).node,d);
              if (dof > -1) {
                for (size_t j=0;j<lnodes[p].size();++j)
                  adj[dof].insert(3*lnodes[p][j]+d);
              }
            }
          } else {
            int dof = A.getEquationForDof(b1.getQuad(p,QINDEX(p,q,dir)).v[i].i,d);
            if (dof > -1) {
              for (size_t j=0;j<lnodes[p].size();++j)
                adj[dof].insert(3*lnodes[p][j]+d);
            }
          }
        }
      }
    }
  }

  Matrix B;
  MatrixOps::fromAdjacency(B,adj,A.getEqns(),3*totalNodes);

  // get a reference element
  Dune::GeometryType gt;
  gt.makeCube(2);
  // get a quadrature rule
  int quadorder = std::max((1.0+n1+0.5)/2.0,(1.0+n2+0.5)/2.0);
  quadorder = std::max(quadorder, 2);
  const Dune::QuadratureRule<ctype,2>& rule = 
                  Dune::QuadratureRules<ctype,2>::rule(gt,quadorder);

  // do the assembly loop
  typename Dune::QuadratureRule<ctype,2>::const_iterator r;
  Dune::DynamicMatrix<ctype> E(ubasis.n,(n1+1)*(n2+1),0.0);
  for (size_t p=0;p<interface.size();++p) {
    const BoundaryGrid::Quad& qi(interface[p]);
    HexGeometry<2,2,GridType> lg(qi);
    for (size_t q=0;q<b1.colSize(p);++q) {
      const BoundaryGrid::Quad& qu = b1.getQuad(p,q);
      HexGeometry<2,2,GridType> hex(qu,gv,dir);
      E = 0;
      for (r = rule.begin(); r != rule.end();++r) {
        ctype detJ = hex.integrationElement(r->position());
        if (detJ < 1.e-4 && r == rule.begin())
          std::cerr << "warning: interface cell (close to) degenerated, |J|=" << detJ << std::endl;

        typename HexGeometry<2,2,GridType>::LocalCoordinate loc = 
                                        lg.local(hex.global(r->position()));
        assert(loc[0] <= 1.0+1.e-4 && loc[0] >= 0.0 && loc[1] <= 1.0+1.e-4 && loc[1] >= 0.0);
        for (int i=0;i<ubasis.n;++i) {
          for (int j=0;j<lbasis.size();++j) {
            E[i][j] += ubasis[i].evaluateFunction(r->position())*
                       lbasis[j].evaluateFunction(loc)*detJ*r->weight();
          }
        }
      }
      // and assemble element contributions
      for (int d=0;d<3;++d) {
        for (int i=0;i<4;++i) {
          MPC* mpc = A.getMPC(qu.v[i].i,d);
          if (mpc) {
            for (size_t n=0;n<mpc->getNoMaster();++n) {
              int indexi = A.getEquationForDof(mpc->getMaster(n).node,d);
              if (indexi > -1) {
                for (size_t j=0;j<lnodes[p].size();++j) {
                  int indexj = lnodes[p][j]*3+d;
                  B[indexi][indexj] += E[i][j];
                }
              }
            }
          } else {
            int indexi = A.getEquationForDof(qu.v[i].i,d);
            if (indexi > -1) {
              for (size_t j=0;j<lnodes[p].size();++j) {
                int indexj = lnodes[p][j]*3+d;
                B[indexi][indexj] += E[i][j];
              }
            }
          }
        }
      }
    }
  }

  return B;
}

  template<class GridType>
void ElasticityUpscale<GridType>::fixPoint(Direction dir,
                                           GlobalCoordinate coord,
                                           const NodeValue& value)
{
  typedef typename GridType::LeafGridView::template Codim<dim>::Iterator VertexLeafIterator;
  const VertexLeafIterator itend = gv.leafView().template end<dim>();

  // make a mapper for codim 0 entities in the leaf grid 
  Dune::LeafMultipleCodimMultipleGeomTypeMapper<GridType,
                                            Dune::MCMGVertexLayout> mapper(gv);

  // iterate over vertices
  for (VertexLeafIterator it = gv.leafView().template begin<dim>(); it != itend; ++it) {
    if (isOnPoint(it->geometry().corner(0),coord)) {
      int indexi = mapper.map(*it);
      A.updateFixedNode(indexi,std::make_pair(dir,value));
    }
  }
}

  template<class GridType>
bool ElasticityUpscale<GridType>::isOnPlane(Direction plane,
                                            GlobalCoordinate coord,
                                            ctype value)
{
  if (plane < X || plane > Z)
    return false;
  int p = log2(plane);
  ctype delta = fabs(value-coord[p]);
  return delta < tol;
}

  template<class GridType>
void ElasticityUpscale<GridType>::fixLine(Direction dir,
                                          ctype x, ctype y,
                                          const NodeValue& value)
{
  typedef typename GridType::LeafGridView::template Codim<dim>::Iterator VertexLeafIterator;
  const VertexLeafIterator itend = gv.leafView().template end<dim>();

  // make a mapper for codim 0 entities in the leaf grid 
  Dune::LeafMultipleCodimMultipleGeomTypeMapper<GridType,
                                            Dune::MCMGVertexLayout> mapper(gv);

  // iterate over vertices
  for (VertexLeafIterator it = gv.leafView().template begin<dim>(); it != itend; ++it) {
    if (isOnLine(dir,it->geometry().corner(0),x,y)) {
      int indexi = mapper.map(*it);
      A.updateFixedNode(indexi,std::make_pair(XYZ,value));
    }
  }
}

  template<class GridType>
bool ElasticityUpscale<GridType>::isOnLine(Direction dir,
                                           GlobalCoordinate coord,
                                           ctype x, ctype y)
{
  if (dir < X || dir > Z)
    return false;
  int ix = int(log2(dir)+1) % 3;
  int iy = int(log2(dir)+2) % 3;
  ctype delta = x-coord[ix];
  if (delta > tol || delta < -tol)
    return false;
  delta = y-coord[iy];
  if (delta > tol || delta < -tol)
    return false;

  return true;
}

  template<class GridType>
bool ElasticityUpscale<GridType>::isOnPoint(GlobalCoordinate coord,
                                            GlobalCoordinate point)
{
  GlobalCoordinate delta = point-coord;
  return delta.one_norm() < tol;
}

  template<class GridType>
void ElasticityUpscale<GridType>::assemble(int loadcase, bool matrix)
{
  const int comp = 3+(dim-2)*3;
  static const int bfunc = 4+(dim-2)*4;

  const LeafIterator itend = gv.leafView().template end<0>();

  Dune::FieldMatrix<ctype,comp,comp> C;
  Dune::FieldMatrix<ctype,dim*bfunc,dim*bfunc> K;
  Dune::FieldVector<ctype,dim*bfunc> ES;
  Dune::FieldVector<ctype,dim*bfunc>* EP=0;
  Dune::FieldVector<ctype,comp> eps0;
  eps0 = 0;
  if (loadcase > -1) {
    EP = &ES;
    eps0[loadcase] = 1;
    A.getLoadVector() = 0;
    b[loadcase] = 0;
  }
  int m=0;
  Dune::FieldMatrix<ctype,dim*bfunc,dim*bfunc>* KP=0;
  if (matrix) {
    A.getOperator() = 0;
    KP = &K;
  }
  for (LeafIterator it = gv.leafView().template begin<0>(); it != itend; ++it) {
    materials[m++]->getConstitutiveMatrix(C);
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
      if (detJ <= 1.e-4) {
        std::cout << "cell " << m << " is (close to) degenerated, detJ " << detJ << std::endl;
        double zdiff=0.0;
        for (int i=0;i<4;++i)
          zdiff = std::max(zdiff, it->geometry().corner(i+4)[2]-it->geometry().corner(i)[2]);
        std::cout << "  - Consider setting ctol larger than " << zdiff << std::endl;
      }

      Dune::FieldMatrix<ctype,comp,dim*bfunc> B;
      E.getBmatrix(B,r->position(),jacInvTra);

      if (matrix) {
        E.getStiffnessMatrix(Aq,B,C,detJ*r->weight());
        K += Aq;
      }

      // load vector
      if (EP) {
        Dune::FieldVector<ctype,dim*bfunc> temp;
        temp = Dune::FMatrixHelp::multTransposed(B,Dune::FMatrixHelp::mult(C,eps0));
        temp *= -detJ*r->weight();
        ES += temp;
      }
    }

    A.addElement(KP,EP,it,(loadcase > -1)?&b[loadcase]:NULL);
  }
}

  template<class GridType>
    template<int comp>
void ElasticityUpscale<GridType>::averageStress(Dune::FieldVector<ctype,comp>& sigma,
                                                const Vector& u, int loadcase)
{
  if (loadcase < 0 || loadcase > 5)
    return;

  static const int bfunc = 4+(dim-2)*4;

  const LeafIterator itend = gv.leafView().template end<0>();

  Dune::FieldMatrix<ctype,comp,comp> C;
  Dune::FieldVector<ctype,comp> eps0;
  eps0 = 0;
  eps0[loadcase] = 1;
  int m=0;
  sigma = 0;
  double volume=0;
  for (LeafIterator it = gv.leafView().template begin<0>(); it != itend; ++it) {
    materials[m++]->getConstitutiveMatrix(C);
    // determine geometry type of the current element and get the matching reference element
    Dune::GeometryType gt = it->type();

    Dune::FieldVector<ctype,bfunc*dim> v;
    A.extractValues(v,u,it);

    // get a quadrature rule of order two for the given geometry type
    const Dune::QuadratureRule<ctype,dim>& rule = Dune::QuadratureRules<ctype,dim>::rule(gt,2);
    for (typename Dune::QuadratureRule<ctype,dim>::const_iterator r = rule.begin();
        r != rule.end() ; ++r) {
      // compute the jacobian inverse transposed to transform the gradients
      Dune::FieldMatrix<ctype,dim,dim> jacInvTra =
        it->geometry().jacobianInverseTransposed(r->position());

      ctype detJ = it->geometry().integrationElement(r->position());
      if (detJ <= 0) // wtf ?
        continue;

      volume += detJ*r->weight();

      Dune::FieldMatrix<ctype,comp,dim*bfunc> B;
      E.getBmatrix(B,r->position(),jacInvTra);

      Dune::FieldVector<ctype,comp> s;
      E.getStressVector(s,v,eps0,B,C);
      s *= detJ*r->weight();
      sigma += s;
    }
  }
  sigma /= volume;
  if (Escale > 0)
    sigma /= Escale/Emin;
}

  template<class GridType>
void ElasticityUpscale<GridType>::loadMaterialsFromGrid(const std::string& file)
{
  typedef std::map<std::pair<double,double>, Material*> MaterialMap;
  MaterialMap cache;
  std::vector<double> Emod;
  std::vector<double> Poiss;
  if (file == "uniform") {
    int cells = gv.size(0);
    Emod.insert(Emod.begin(),cells,100.f);
    Poiss.insert(Poiss.begin(),cells,0.38f);
  } else {
    Opm::EclipseGridParser parser(file,false);
    if (parser.hasField("YOUNGMOD") && parser.hasField("POISSONMOD")) {
      Emod = parser.getFloatingPointValue("YOUNGMOD");
      Poiss = parser.getFloatingPointValue("POISSONMOD");
    } else if (parser.hasField("PERMX") && parser.hasField("PORO")) {
      Emod = parser.getFloatingPointValue("PERMX");
      Poiss = parser.getFloatingPointValue("PORO");
    } else if (parser.hasField("LAMEMOD") && parser.hasField("SHEARMOD")) {
      std::vector<double> lame = parser.getFloatingPointValue("LAMEMOD");
      std::vector<double> shear = parser.getFloatingPointValue("SHEARMOD");
      for (size_t i=0;i<lame.size();++i) {
        Emod[i]  = shear[i]*(3*lame[i]+2*shear[i])/(lame[i]+shear[i]);
        Poiss[i] = 0.5*lame[i]/(lame[i]+shear[i]);
      }
    } else if (parser.hasField("BULKMOD") && parser.hasField("SHEARMOD")) {
      std::vector<double> bulk = parser.getFloatingPointValue("BULKMOD");
      std::vector<double> shear = parser.getFloatingPointValue("SHEARMOD");
      for (size_t i=0;i<bulk.size();++i) {
        Emod[i]  = 9*bulk[i]*shear[i]/(3*bulk[i]+shear[i]);
        Poiss[i] = 0.5*(3*bulk[i]-2*shear[i])/(3*bulk[i]+shear[i]);
      }
    } else {
      std::cerr << "No material data found in eclipse file, aborting" << std::endl;
      exit(1);
    }
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
    if ((it = cache.find(std::make_pair(Emod[k],Poiss[k]))) != cache.end())
    {
      volume[it->second] += gv.cellVolume(i);
      materials.push_back(it->second);
    }
    else {
      Material* mat = new Isotropic(j++,Emod[k],Poiss[k]);
      cache.insert(std::make_pair(std::make_pair(Emod[k],Poiss[k]),mat));
      volume.insert(std::make_pair(mat,gv.cellVolume(i)));
      materials.push_back(mat);
    }
  }
  std::cout << "Number of materials: " << cache.size() << std::endl;
  // statistics
  double totalvolume=0;
  for (std::map<Material*,double>::iterator it  = volume.begin(); 
                                            it != volume.end(); ++it) 
    totalvolume += it->second;

  int i=0;
  for (MaterialMap::iterator it = cache.begin(); it != cache.end(); ++it, ++i) {
    std::cout << "  Material" << i+1 << ": " << 100.f*volume[it->second]/totalvolume << '%' << std::endl;
    volumeFractions.push_back(volume[it->second]/totalvolume);
  }
}

  template<class GridType>
void ElasticityUpscale<GridType>::loadMaterialsFromRocklist(const std::string& file,
                                                            const std::string& rocklist)
{
  std::vector<Material*> cache;
  // parse the rocklist
  std::ifstream f;
  f.open(rocklist.c_str());
  int mats;
  f >> mats;
  for (int i=0;i<mats;++i) {
    std::string file;
    f >> file;
    cache.push_back(Material::create(i+1,file));
  }

  // scale E modulus of materials
  if (Escale > 0) {
    Emin=1e10;
    for (size_t i=0;i<cache.size();++i)
      Emin = std::min(Emin,((Isotropic*)cache[i])->getE());
    for (size_t i=0;i<cache.size();++i) {
      double E = ((Isotropic*)cache[i])->getE();
      ((Isotropic*)cache[i])->setE(E*Escale/Emin);
    }
  }
  std::vector<double> volume;
  volume.resize(cache.size());
  if (file == "uniform") {
    for (int i=0;i<gv.size(0);++i)
      materials.push_back(cache[0]);
    volume[0] = 1;
  } else {
    Opm::EclipseGridParser parser(file,false);
    std::vector<int> satnum = parser.getIntegerValue("SATNUM");
    std::vector<int> cells = gv.globalCell();
    for (size_t i=0;i<cells.size();++i) {
      int k = cells[i];
      if (satnum[k]-1 >= cache.size()) {
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
    std::cout << "  Material" << i+1 << ": " << 100.f*volume[i]/totalvolume << '%' << std::endl;
    volumeFractions.push_back(volume[i]/totalvolume);
  }
}

  template<class GridType>
void ElasticityUpscale<GridType>::fixCorners(const double* min,
                                             const double* max)
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

  template<class GridType>
void ElasticityUpscale<GridType>::periodicBCs(const double* min, 
                                              const double* max)
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

  template<class GridType>
void ElasticityUpscale<GridType>::periodicBCsMortar(const double* min, 
                                                    const double* max,
                                                    int n1, int n2,
                                                    int p1, int p2)
{
  // this method
  // 1. fixes the primal corner dofs
  // 2. establishes strong couplings (MPC) on top and bottom
  // 3. extracts and establishes a quad grid for the left/right/front/back sides
  // 4. establishes grids for the dual dofs
  // 5. calculates the coupling matrix L1 between the left/right sides
  // 6. calculates the coupling matrix L2 between the front/back sides

  // step 1
  fixCorners(min,max);
  
  // step 2
  slave.push_back(extractFace(Z,max[2]));
  BoundaryGrid bottom = extractMasterFace(Z,min[2]);
  periodicPlane(Z,XYZ,slave[0],bottom);
  A.initForAssembly();

  // step 3
  master.push_back(extractMasterFace(X,min[0],LEFT,true));
  master.push_back(extractMasterFace(X,max[0],RIGHT,true));
  master.push_back(extractMasterFace(Y,min[1],LEFT,true));
  master.push_back(extractMasterFace(Y,max[1],RIGHT,true));

  std::cout << "Xsides " << master[0].size() << " " << master[1].size() << std::endl
            << "Ysides " << master[2].size() << " " << master[3].size() << std::endl
            << "Zmaster " << bottom.size() << std::endl;
  std::cout << "Establish YZ multiplier grid with " << n2 << "x1" << " elements" << std::endl;

  // step 4
  BoundaryGrid::FaceCoord fmin,fmax;
  fmin[0] = min[1]; fmin[1] = min[2];
  fmax[0] = max[1]; fmax[1] = max[2];
  BoundaryGrid lambdax = BoundaryGrid::uniform(fmin,fmax,n2,1,true);

  fmin[0] = min[0]; fmin[1] = min[2];
  fmax[0] = max[0]; fmax[1] = max[2];
  std::cout << "Establish XZ multiplier grid with " << n1 << "x1" << " elements" << std::endl;
  BoundaryGrid lambday = BoundaryGrid::uniform(fmin,fmax,n1,1,true);

  // step 5
  Matrix L1 = findLMatrixMortar(master[0],lambdax,0, p2, 1);
  Matrix L2 = findLMatrixMortar(master[1],lambdax,0, p2, 1);
  L.push_back(MatrixOps::Axpy(L1,L2,-1));

  // step 6
  Matrix L3 = findLMatrixMortar(master[2],lambday,1, p1, 1);
  Matrix L4 = findLMatrixMortar(master[3],lambday,1, p1, 1);
  L.push_back(MatrixOps::Axpy(L3,L4,-1));
}

  template<class GridType>
void ElasticityUpscale<GridType>::periodicBCsLLM(const double* min, 
                                                 const double* max,
                                                 int n1, int n2)
{
  // this method
  // 1. fixes the primal corner dofs
  // 2. fixes the primal dofs on the skeleton
  // 3. establishes strong couplings (MPC) on top and bottom
  // 4. extracts a point grid for the left/right/front/back sides,
  //    and establishes a uniform interface grid in each direction
  // 5. calculates the coupling matrices B1-4 and L1-4

  // step 1
  fixCorners(min,max);

  // step 2
  fixLine(X,min[1],min[2]);
  fixLine(X,max[1],min[2]);
  fixLine(X,min[1],max[2]);
  fixLine(X,max[1],max[2]);

  fixLine(Y,min[0],min[2]);
  fixLine(Y,max[0],min[2]);
  fixLine(Y,min[0],max[2]);
  fixLine(Y,max[0],max[2]);

  fixLine(Z,min[0],min[1]);
  fixLine(Z,max[0],min[1]);
  fixLine(Z,min[0],max[1]);
  fixLine(Z,max[0],max[1]);

  // step 3
  std::vector<BoundaryGrid::Vertex> Zs = extractFace(Z,max[2]);
  BoundaryGrid Zm = extractMasterFace(Z,min[2]);
  periodicPlane(Z,XYZ,Zs,Zm);
  A.initForAssembly();

  // step 4
  slave.push_back(extractFace(X,min[0]));
  slave.push_back(extractFace(X,max[0]));
  slave.push_back(extractFace(Y,min[1]));
  slave.push_back(extractFace(Y,max[1]));

  Dune::LeafMultipleCodimMultipleGeomTypeMapper<GridType,
                                            Dune::MCMGVertexLayout> mapper(gv);

  BoundaryGrid::FaceCoord fmin,fmax;
  // YZ plane
  fmin[0] = min[1]; fmin[1] = min[2];
  fmax[0] = max[1]; fmax[1] = max[2];
  std::cout << "Establish YZ interface grid with " << n2 << "x1" << " elements" << std::endl;
  master.push_back(BoundaryGrid::uniform(fmin,fmax,n2,1));
  // XZ plane
  fmin[0] = min[0]; fmin[1] = min[2];
  fmax[0] = max[0]; fmax[1] = max[2];
  std::cout << "Establish XZ interface grid with " << n1 << "x1" << " elements" << std::endl;
  master.push_back(BoundaryGrid::uniform(fmin,fmax,n1,1));

  // step 5
  std::map<int,BoundaryGrid::Vertex> m;
  // find matching coefficients
  for (size_t i=0;i<master.size();++i) {
    for (int j=0;j<2;++j) {
      m.clear();
      for (size_t k=0;k<slave[i*2+j].size();++k) {
        BoundaryGrid::Vertex c;
        if (master[i].find(c,slave[i*2+j][k])) {
          m.insert(std::make_pair(slave[i*2+j][k].i,c));
        } else
          assert(0);
      }
      B.push_back(findBMatrixLLM(m));
      L.push_back(findLMatrixLLM(m,master[i]));
    }
  }
}

 template<class GridType>
void ElasticityUpscale<GridType>::setupPreconditioner()
{
  // setup subdomain maps
  typename OverlappingSchwarz::subdomain_vector rows;
  const LeafIterator itend = gv.leafView().template end<0>();
  const LeafIndexSet& set = gv.leafView().indexSet();
    
  rows.resize(gv.logicalCartesianSize()[0]*gv.logicalCartesianSize()[1]);
  int cell=0, currdomain=0;
  for (LeafIterator it = gv.leafView().template begin<0>(); it != itend; ++it, ++cell) {
    if (cell / gv.logicalCartesianSize()[2] > 0 
        && cell % gv.logicalCartesianSize()[2] == 0)
      currdomain++;
    for (size_t i=0;i<8;++i) {
      int idx=set.subIndex(*it,i,dim);
      for (int d=0;d<3;++d) {
        MPC* mpc = A.getMPC(idx,d);
        if (mpc) {
          for (size_t n=0;n<mpc->getNoMaster();++n) {
            int row = A.getEquationForDof(mpc->getMaster(n).node,d);
            if (row > -1)
              rows[currdomain].insert(row);
          }
        } else {
          int row = A.getEquationForDof(idx,d);
          if (row > -1)
            rows[currdomain].insert(row);
        }
      }
    }
  }
  ovl = new OverlappingSchwarz(A.getOperator(),rows,1,false);
}

  template<class GridType>
void ElasticityUpscale<GridType>::setupSolvers(Solver solver)
{
  if (!B.empty() && A.getOperator().N() == A.getEqns()) { // LLM in use
    // The LLM linear system is of the form
    // [A    B1  B2    B3   B4    0   0] [   u]   [b]
    // [B1'   0   0     0    0   L1   0] [ l_1]   [0]
    // [B2'   0   0     0    0   L2   0] [ l_2]   [0]
    // [B3'   0   0     0    0    0  L3] [ l_3] = [0]
    // [B4'   0   0     0    0    0  L4] [ l_4]   [0]
    // [ 0   L1'  L2'   0    0    0   0] [ub_1]   [0]
    // [ 0    0   0   L3'  L4'    0   0] [ub_2]   [0]
    int r = A.getOperator().N();
    int c = A.getOperator().M();
    for (size_t i=0;i<B.size();++i) {
      A.getOperator() = MatrixOps::augment(A.getOperator(),B[i],0,c,true);
      c += B[i].M();
    }
    for (size_t i=0;i<L.size();++i) {
      A.getOperator() = MatrixOps::augment(A.getOperator(),L[i],r,c,true);
      r += L[i].N();
      if (i % 2 == 1)
        c += L[i].M();
    }
    int siz=A.getOperator().N();
    for (int i=0;i<6;++i)
      b[i].resize(A.getOperator().N());
  }
  // Mortar in use
  if (B.empty() && !L.empty() && A.getOperator().N() == A.getEqns()) { 
    // The Mortar linear system is of the form
    // [A   L1 L2] [  u]   [b]
    // [L1'  0  0] [l_1] = [0]
    // [L2'  0  0] [l_2]   [0]
    int c = A.getOperator().M();
    for (size_t i=0;i<L.size();++i) { 
      A.getOperator() = MatrixOps::augment(A.getOperator(),L[i],0,c,true);
      c += L[i].M();
    }
    int siz=A.getLoadVector().size();
    for (int i=0;i<6;++i) {
      b[i].resize(A.getOperator().N());
      for (int j=siz;j<b[i].size();++j)
        b[i][j] = 0;
    }
    if (L.empty() && B.empty()) { // MPC
      // overlapping schwarz is much more memory efficient
      // and usually more effective
      setupPreconditioner(); 
    }
  }
  if (solver == SLU) {
#if HAVE_SUPERLU
    if (!slu)
      slu = new Dune::SuperLU<Matrix>(A.getOperator(),false);
#else
    std::cerr << "SuperLU solver not enabled" << std::endl;
    exit(1);
#endif
  } else if (solver == CG) {
    if (!op)
      op = new Dune::MatrixAdapter<Matrix,Vector,Vector>(A.getOperator());
    if (!ovl && !ilu)
      ilu = new Dune::SeqILU0<Matrix,Vector,Vector>(A.getOperator(),0.92);
    if (!cgsolver) {
      if (ovl)
        cgsolver = new Dune::CGSolver<Vector>(*op, *ovl, tol, 5000, 1);
      else
        cgsolver = new Dune::CGSolver<Vector>(*op, *ilu, tol, 5000, 1);
    }
  }
}

  template<class GridType>
void ElasticityUpscale<GridType>::solve(Solver solver, double tol, int loadcase)
{
  Dune::InverseOperatorResult r;

  // initialize u to some arbitrary value
  u[loadcase].resize(A.getOperator().N(), false);
  u[loadcase] = 2.0;

#if HAVE_SUPERLU
  if (solver == SLU) {
    slu->apply(u[loadcase], b[loadcase], r);
  }
  else 
#endif
  if (solver == CG) {
    cgsolver->apply(u[loadcase], b[loadcase], r);
  }
  std::cout << "\t solution norm: " << u[loadcase].two_norm() << std::endl;
}
