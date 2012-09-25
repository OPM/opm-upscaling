//==============================================================================
//!
//! \file mpc.cpp
//!
//! \date Oct 1 2007
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Representation of multi-point constraint (MPC) equations.
//!
//==============================================================================

#include "mpc.hh"

bool MPC::Less::compareSlaveDofOnly = false;


//! \brief Global operator for comparing two MPC::DOF objects.

bool operator< (const MPC::DOF& lhs, const MPC::DOF& rhs)
{
  if (lhs.node < rhs.node) return true;
  if (lhs.node > rhs.node) return false;
  if (lhs.dof  < rhs.dof)  return true;
  if (lhs.dof  > rhs.dof)  return false;

  if (MPC::Less::compareSlaveDofOnly)
    return false; // ignore coefficient differences, if any
  else
    return lhs.coeff < rhs.coeff ? true : false;
}


bool MPC::Less::operator() (const MPC* lhs, const MPC* rhs) const
{
  if (!rhs) return false;
  if (!lhs) return true;

  if (lhs->getSlave() < rhs->getSlave()) return true;
  if (rhs->getSlave() < lhs->getSlave()) return false;

  if (compareSlaveDofOnly) return false; // ignore master differences, if any

  size_t lMaster = lhs->getNoMaster();
  size_t rMaster = rhs->getNoMaster();
  for (size_t i = 0; i < lMaster && i < rMaster; i++)
  {
    if (lhs->getMaster(i) < rhs->getMaster(i)) return true;
    if (rhs->getMaster(i) < lhs->getMaster(i)) return false;
  }

  return lMaster < rMaster ? true : false;
}


//! \brief Global stream operator printing a constraint equation.

std::ostream& operator<< (std::ostream& s, const MPC& mpc)
{
  s <<"Slave "<< mpc.slave.node <<","<< mpc.slave.dof;
  if (mpc.slave.coeff != double(0))
    s <<" = "<< mpc.slave.coeff;
  else if (mpc.master.empty())
    return s <<" = 0"<< std::endl;
  else
    s <<" = ";

  for (size_t i = 0; i < mpc.master.size(); i++)
  {
    if (i == 0 && mpc.slave.coeff == double(0))
      s << mpc.master[i].coeff;
    else if (mpc.master[i].coeff >= double(0))
      s <<" + "<< mpc.master[i].coeff;
    else
      s <<" - "<< -mpc.master[i].coeff;
    s <<"*("<< mpc.master[i].node <<","<< mpc.master[i].dof <<")";
  }
  return s << std::endl;
}
