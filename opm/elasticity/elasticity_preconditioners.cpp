//==============================================================================
//!
//! \file elasticity_preconditioners.cpp
//!
//! \date Aug 30 2013
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Preconditioners for elasticity upscaling
//!
//==============================================================================

#include "config.h"

#include "elasticity_preconditioners.hpp"

namespace Opm {
namespace Elasticity {

std::shared_ptr<FastAMG::type> FastAMG::setup(int pre, int post, int target,
                                              int zcells,
                                              std::shared_ptr<Operator>& op,
                                              const Dune::CpGrid& gv,
                                              ASMHandler<Dune::CpGrid>& A,
                                              bool& copy)
{
  AMG1<JACSmoother>::Criterion crit;
  crit.setCoarsenTarget(target);
  crit.setGamma(1);
  crit.setDefaultValuesIsotropic(3, zcells);

  std::cout << "\t collapsing 2x2x" << zcells << " cells per level" << std::endl;
  copy = true;
  return std::shared_ptr<type>(new type(*op, crit));
}

Schwarz::type* Schwarz::setup2(std::shared_ptr<Operator>& op,
                               const Dune::CpGrid& gv,
                               ASMHandler<Dune::CpGrid>& A, bool& copy)
{
  const int cps = 1;
  Schwarz::type::subdomain_vector rows;
  int nel1 = gv.logicalCartesianSize()[0];
  int nel2 = gv.logicalCartesianSize()[1];
  rows.resize(nel1/cps*nel2/cps);

  auto set = gv.leafView().indexSet();
  for (auto it  = gv.leafView().begin<0>(), e = gv.leafView().end<0>();
            it != e; ++it) {
    std::array<int, 3> ijk;
    gv.getIJK(set.index(*it), ijk);
    const int rowix = (ijk[0]/cps) + (nel1/cps)*(ijk[1]/cps);
    // loop over nodes
    for (int n=0;n<8;++n) {
      int idx = set.subIndex(*it, n, 3);
      for (int m=0;m<3;++m) {
        const MPC* mpc = A.getMPC(idx, m);
        if (mpc) {
          for (size_t q=0;q<mpc->getNoMaster();++q) {
            int idx2 = A.getEquationForDof(mpc->getMaster(q).node, m);
            if (idx2 > -1)
              rows[rowix].insert(idx2);
          }
        } else {
          if (A.getEquationForDof(idx, m) > -1)
            rows[rowix].insert(A.getEquationForDof(idx, m));
        }
      }
    }
  }

  copy = false;
  return new type(op->getmat(), rows, 1.0, false);
}

}
}
