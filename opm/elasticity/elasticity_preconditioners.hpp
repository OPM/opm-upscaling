//==============================================================================
//!
//! \file elasticity_preconditioners.hpp
//!
//! \date Aug 30 2013
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Preconditioners for elasticity upscaling
//!
//==============================================================================
#ifndef ELASTICITY_PRECONDITIONERS_HPP_
#define ELASTICITY_PRECONDITIONERS_HPP_

#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixmatrix.hh>
#include <dune/istl/ilu.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/grid/CpGrid.hpp>

#include <opm/elasticity/asmhandler.hpp>
#include <opm/elasticity/matrixops.hpp>

#include <dune/istl/superlu.hh>
#include <dune/istl/umfpack.hh>
#include <dune/istl/paamg/amg.hh>
#include <dune/istl/paamg/fastamg.hh>
#include <dune/istl/paamg/twolevelmethod.hh>
#include <dune/istl/overlappingschwarz.hh>

namespace Opm {
namespace Elasticity {

#if defined(HAVE_UMFPACK)
typedef Dune::UMFPack<Matrix> LUSolver;
#elif defined(HAVE_SUPERLU)
typedef Dune::SuperLU<Matrix> LUSolver;
#else
static_assert("Enable either SuperLU or UMFPACK");
#endif

//! \brief A linear operator
typedef Dune::MatrixAdapter<Matrix,Vector,Vector> Operator;

//! \brief SSOR AMG smoother
typedef Dune::SeqSSOR<Matrix, Vector, Vector> SSORSmoother;

//! \brief GJ AMG smoother
typedef Dune::SeqJac<Matrix, Vector, Vector> JACSmoother;

//! \brief ILU0 AMG smoother
typedef Dune::SeqILU0<Matrix, Vector, Vector> ILUSmoother;

//! \brief Schwarz + ILU0 AMG smoother
typedef Dune::SeqOverlappingSchwarz<Matrix,Vector,
                              Dune::SymmetricMultiplicativeSchwarzMode, LUSolver> SchwarzSmoother;

//! \brief Overlapping Schwarz preconditioner
struct Schwarz {
  typedef Dune::SeqOverlappingSchwarz<Matrix, Vector,
                                      Dune::SymmetricMultiplicativeSchwarzMode,
                                      LUSolver> type;
  //! \brief Setup preconditioner
  //! \param[in] pre The number of pre-smoothing steps
  //! \param[in] post The number of post-smoothing steps
  //! \param[in] target The coarsening target
  //! \param[in] zcells The wanted number of cells to collapse in z per level
  //! \param[in] op The linear operator
  //! \param[in] gv The cornerpoint grid
  //! \param[out] thread Whether or not to clone for threads
  static std::shared_ptr<type>
                setup(int pre, int post, int target, int zcells,
                      std::shared_ptr<Operator>& op, const Dune::CpGrid& gv,
                      ASMHandler<Dune::CpGrid>& A, bool& copy)
  {
    return std::shared_ptr<type>(setup2(op, gv, A, copy));
  }

  //! \brief Setup preconditioner
  //! \param[in] op The linear operator
  //! \param[in] gv The cornerpoint grid
  //! \param[in] A The ASMHandler for the elasticity operator(s)
  //! \param[out] copy Whether or not to clone for threads
  static type* setup2(std::shared_ptr<Operator>& op, const Dune::CpGrid& gv,
                      ASMHandler<Dune::CpGrid>& A, bool& copy);
};

//! \brief An AMG
template<class Smoother>
struct AMG1 {
  //! \brief The coupling metric used in the AMG
  typedef Dune::Amg::FirstDiagonal CouplingMetric;

  //! \brief The coupling criterion used in the AMG
  typedef Dune::Amg::SymmetricCriterion<Matrix, CouplingMetric> CritBase;

  //! \brief The coarsening criterion used in the AMG
  typedef Dune::Amg::CoarsenCriterion<CritBase> Criterion;

  typedef Dune::Amg::AMG<Operator, Vector, Smoother> type;

  //! \brief Setup preconditioner
  //! \param[in] pre The number of pre-smoothing steps
  //! \param[in] post The number of post-smoothing steps
  //! \param[in] target The coarsening target
  //! \param[in] zcells The wanted number of cells to collapse in z per level
  //! \param[in] op The linear operator
  //! \param[in] gv The cornerpoint grid
  //! \param[out] thread Whether or not to clone for threads
  static std::shared_ptr<type>
                setup(int pre, int post, int target, int zcells,
                      std::shared_ptr<Operator>& op, const Dune::CpGrid& gv,
                      ASMHandler<Dune::CpGrid>& A, bool& copy)
  {
    Criterion crit;
    typename AMG1<Smoother>::type::SmootherArgs args;
    args.relaxationFactor = 1.0;
    crit.setCoarsenTarget(target);
    crit.setGamma(1);
    crit.setNoPreSmoothSteps(pre);
    crit.setNoPostSmoothSteps(post);
    crit.setDefaultValuesIsotropic(3, zcells);

    std::cout << "\t collapsing 2x2x" << zcells << " cells per level" << std::endl;
    copy = true;
    return std::shared_ptr<type>(new type(*op, crit, args));
  }
};

//! \brief A FastAMG 
struct FastAMG {
  typedef Dune::Amg::FastAMG<Operator, Vector> type;

  //! \brief Setup preconditioner
  //! \param[in] pre The number of pre-smoothing steps
  //! \param[in] post The number of post-smoothing steps
  //! \param[in] target The coarsening target
  //! \param[in] zcells The wanted number of cells to collapse in z per level
  //! \param[in] op The linear operator
  //! \param[in] gv The cornerpoint grid
  //! \param[out] thread Whether or not to clone for threads
  static std::shared_ptr<type>
                setup(int pre, int post, int target, int zcells,
                      std::shared_ptr<Operator>& op, const Dune::CpGrid& gv,
                      ASMHandler<Dune::CpGrid>& A, bool& copy);
};


//! \brief A two-level method with a coarse AMG solver
  template<class Smoother>
struct AMG2Level {
  //! \brief AMG transfer policy
  typedef Dune::Amg::AggregationLevelTransferPolicy<Operator,
                             typename AMG1<Smoother>::Criterion> TransferPolicy;

  typedef Dune::Amg::LevelTransferPolicy<Operator, Operator> LevelTransferPolicy;

  typedef Dune::Amg::OneStepAMGCoarseSolverPolicy<Operator, Smoother,
                                             typename AMG1<Smoother>::Criterion> CoarsePolicy;

  typedef typename Dune::Amg::SmootherTraits<Smoother>::Arguments SmootherArgs;
  
  typedef Dune::Amg::TwoLevelMethod<Operator, CoarsePolicy, Schwarz::type> type;

  //! \brief Setup preconditioner
  //! \param[in] pre The number of pre-smoothing steps
  //! \param[in] post The number of post-smoothing steps
  //! \param[in] target The coarsening target
  //! \param[in] zcells The wanted number of cells to collapse in z per level
  //! \param[in] op The linear operator
  //! \param[in] gv The cornerpoint grid
  //! \param[out] thread Whether or not to clone for threads
  static std::shared_ptr<type>
                setup(int pre, int post, int target, int zcells,
                      std::shared_ptr<Operator>& op, const Dune::CpGrid& gv,
                      ASMHandler<Dune::CpGrid>& A, bool& copy)
  {
    typename AMG1<Smoother>::Criterion crit;
    SmootherArgs args;
    args.relaxationFactor = 1.0;
    crit.setCoarsenTarget(target);
    crit.setGamma(1);
    crit.setNoPreSmoothSteps(pre);
    crit.setNoPostSmoothSteps(post);
    crit.setDefaultValuesIsotropic(3, zcells);
    CoarsePolicy coarsePolicy(args, crit);
    TransferPolicy policy(crit);
    Dune::shared_ptr<Schwarz::type> fsp(Schwarz::setup2(op, gv, A, copy));
    copy = true;
    return std::shared_ptr<type>(new type(*op, fsp, policy, coarsePolicy, pre, post));
  }
};

}
}

#endif
