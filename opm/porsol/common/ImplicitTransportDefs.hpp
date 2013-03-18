
#ifndef OPENRS_IMPLICITTRANSPORTDEFS_HEADER
#define OPENRS_IMPLICITTRANSPORTDEFS_HEADER

#include <opm/porsol/common/ReservoirPropertyCapillary.hpp>
#include <opm/porsol/common/LinearSolverISTL.hpp>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <dune/grid/common/GridAdapter.hpp>
#include <opm/core/linalg/sparse_sys.h>
#include <opm/core/pressure/tpfa/ifs_tpfa.h>
#include <opm/core/pressure/tpfa/trans_tpfa.h>
#include <opm/core/transport/implicit/NormSupport.hpp>
#include <opm/core/transport/implicit/ImplicitTransport.hpp>

#if 0
template <class Ostream, class Collection>
Ostream&
operator<<(Ostream& os, const Collection& c)
{
    typedef typename Collection::value_type VT;

    os << "[ ";
    std::copy(c.begin(), c.end(), ::std::ostream_iterator<VT>(os, " "));
    os << "]";

    return os;
}
#endif
namespace Opm{
template <class Vector>
class MaxNormStl {
	public:
	    static double
	    norm(const Vector& v) {
		return  ImplicitTransportDefault::AccumulationNorm <Vector,  ImplicitTransportDefault::MaxAbs>::norm(v);
	    }
};

template <class Vector>
class MaxNormDune {
public:
    static double
    norm(const Vector& v) {
        return v.infinity_norm();
    }
};

template <int np = 2>
class ReservoirState { 
public:
    ReservoirState(const UnstructuredGrid* g)
        : press_ (g->number_of_cells),
          fpress_(g->number_of_faces),
          flux_  (g->number_of_faces),
          sat_   (np * g->number_of_cells)
    {}

    ::std::vector<double>& pressure    () { return press_ ; }
    ::std::vector<double>& facepressure() { return fpress_; }
    ::std::vector<double>& faceflux    () { return flux_  ; }
    ::std::vector<double>& saturation  () { return sat_   ; }

    const ::std::vector<double>& faceflux    () const { return flux_; }
    const ::std::vector<double>& saturation  () const { return sat_ ; }
 
private:
    ::std::vector<double> press_ ;
    ::std::vector<double> fpress_;
    ::std::vector<double> flux_  ;
    ::std::vector<double> sat_   ;
};
class Rock {
public:
    Rock(::std::size_t nc, ::std::size_t dim)
        : dim_ (dim           ),
          perm_(nc * dim * dim),
          poro_(nc            ) {}

    const ::std::vector<double>& perm() const { return perm_; }
    const ::std::vector<double>& poro() const { return poro_; }

    void
    perm_homogeneous(double k) {
        setVector(0.0, perm_);

        const ::std::size_t d2 = dim_ * dim_;

        for (::std::size_t c = 0, nc = poro_.size(); c < nc; ++c) {
            for (::std::size_t i = 0; i < dim_; ++i) {
                perm_[c*d2 + i*(dim_ + 1)] = k;
            }
        }
    }

    void
    poro_homogeneous(double phi) {
        setVector(phi, poro_);
    }

private:
    void
    setVector(double x, ::std::vector<double>& v) {
        ::std::fill(v.begin(), v.end(), x);
    }

    ::std::size_t         dim_ ;
    ::std::vector<double> perm_;
    ::std::vector<double> poro_;
};

//template <class P>
class TwophaseFluidWrapper {
public:
	TwophaseFluidWrapper(const Opm::ReservoirPropertyCapillary<3>& r)
	        : r_(r)
	{}
	//TwophaseFluidWrapper(){}
	//TwophaseFluidWrapper(TwophaseFluidWrapper){}

    /*
    void init(const Opm::ReservoirPropertyCapillary<3>& r)
    {
        r_ = r;
    }
    */

    template <class Sat ,
              class Mob ,
              class DMob>
    void
    mobility(int c, const Sat& s, Mob& mob, DMob& dmob) const {
        const double s1 = s[0];

        r_.phaseMobilities     (c, s1, mob );
        r_.phaseMobilitiesDeriv(c, s1, dmob);
    }
    template <class Sat ,
              class PC ,
              class DPC>
    void
    pc(int c, const Sat& s, PC& pc, DPC& dpc) const {
    	const double s1 = s[0];
    	pc = r_.capillaryPressure(c, s1);
    	dpc= r_.capillaryPressureDeriv(c, s1);
    }
    double density(int p) const {
        if (p == 0) {
            return r_.densityFirstPhase();
        } else {
            return r_.densitySecondPhase();
        }
    }

    double s_min(int c) const {
        return r_.s_min(c);
    }
    double s_max(int c) const {
    	return r_.s_max(c);
    }


private:
    Opm::ReservoirPropertyCapillary<3> r_;
};
    
class LinearSolverBICGSTAB {
public:
	typedef Dune::FieldVector<double, 1>    ScalarVectorBlockType;
	typedef Dune::FieldMatrix<double, 1, 1> ScalarMatrixBlockType;

	typedef Dune::BlockVector<ScalarVectorBlockType> ScalarBlockVector;
	typedef Dune::BCRSMatrix <ScalarMatrixBlockType> ScalarBCRSMatrix;
    void
    solve(const ScalarBCRSMatrix&  A,
          const ScalarBlockVector& b,
          ScalarBlockVector&       x) {

        Dune::MatrixAdapter<ScalarBCRSMatrix,
                            ScalarBlockVector,
                            ScalarBlockVector> opA(A);

        Dune::SeqILU0<ScalarBCRSMatrix,
                      ScalarBlockVector,
                      ScalarBlockVector> precond(A, 1.0);

        int maxit  = A.N();
        double tol = 5.0e-7;
        int verb   = 0;

        Dune::BiCGSTABSolver<ScalarBlockVector>
            solver(opA, precond, tol, maxit, verb);

        ScalarBlockVector           bcpy(b);
        Dune::InverseOperatorResult res;
        solver.apply(x, bcpy, res);
    }
};

class LinearSolverISTLAMG {
public:
    LinearSolverISTLAMG()
    {
        Opm::parameter::ParameterGroup params;

        params.insertParameter("linsolver_tolerance",
                               boost::lexical_cast<std::string>(5.0e-9));

        params.insertParameter("linsoler_verbosity",
                               boost::lexical_cast<std::string>(1));

        params.insertParameter("linsolver_type",
                               boost::lexical_cast<std::string>(1));

        ls_.init(params);
    }

    void
    solve(struct CSRMatrix* A,
          const double*     b,
          double*           x)
    {
        ls_.solve(A->m, A->nnz, A->ia, A->ja, A->sa, b, x);
    }

private:
    Opm::LinearSolverISTL ls_;
};

class TransportSource {
public:
    TransportSource() : nsrc(0),pf(0) {}

    int                   nsrc      ;
    int                   pf      ;
    ::std::vector< int  > cell      ;
    ::std::vector<double> pressure  ;
    ::std::vector<double> flux      ;
    ::std::vector<double> saturation;
    ::std::vector<double> periodic_cells;
    ::std::vector<double> periodic_faces;
};
template <class Arr>
void
append_transport_source(int c, double p, double v, const Arr& s,
                        TransportSource& src)
{
    src.cell      .push_back(c);
    src.pressure  .push_back(p);
    src.flux      .push_back(v);
    src.saturation.insert(src.saturation.end(),
                          s.begin(), s.end());
    ++src.nsrc;
}
void
compute_porevolume(const UnstructuredGrid*        g,
                   const Rock&          rock,
                   std::vector<double>& porevol)
{
    const ::std::vector<double>& poro = rock.poro();

    assert (poro.size() == (::std::size_t)(g->number_of_cells));

    porevol.resize(rock.poro().size());

    ::std::transform(poro.begin(), poro.end(),
                     g->cell_volumes,
                     porevol.begin(),
                     ::std::multiplies<double>());
}
}
#endif // /OPENRS_IMPLICITTRANSPORTDEFS_HEADER
