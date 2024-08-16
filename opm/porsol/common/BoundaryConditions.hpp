//===========================================================================
//
// File: BoundaryConditions.hpp
//
// Created: Mon Jun 29 15:19:42 2009
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//            Bård Skaflestad     <bard.skaflestad@sintef.no>
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2009, 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2009, 2010 Statoil ASA.

  This file is part of The Open Reservoir Simulator Project (OpenRS).

  OpenRS is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OpenRS is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OpenRS.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPENRS_BOUNDARYCONDITIONS_HEADER
#define OPENRS_BOUNDARYCONDITIONS_HEADER


#include <vector>
#include <ostream>
#include <type_traits>
#include <opm/common/ErrorMacros.hpp>
#include <dune/common/fvector.hh>

namespace Opm
{

    /// @brief A class for building boundary conditions in a uniform way.
    template <typename T>
    class BCBase
    {
    public:
	/// @brief Enum for the allowed boundary condition types.
	/// So far, we support Dirichlet, Neumann and Periodic conditions.
	/// In this class, these are just tags, it's up to the code using it
	/// to attach meaning to them.
        enum BCType { Dirichlet, Neumann, Periodic };

	/// @brief Write type and value to an output stream.
	/// @tparam traits character type.
	/// @tparam traits character traits.
	/// @param os output stream.
	template<typename charT, class traits>
	void write(std::basic_ostream<charT,traits>& os) const
	{
	    os << "Type: " << type_ << "   Value: " << value_;
	}

    protected: // methods
	/// @brief Default constructor, that makes a Neumann condition with value 0.0.
        BCBase()
            : type_(Neumann), value_(0.0)
        {
        }
	/// @brief Constructor taking a type and value.
	/// @param type the condition type.
	/// @param value the condition value.
        BCBase(BCType type, T value)
            : type_(type), value_(value)
        {
        }
	/// @brief Type query.
	/// @return true if the type is Dirichlet.
        bool isDirichlet() const
        {
            return type_ == Dirichlet;
        }
	/// @brief Type query.
	/// @return true if the type is Neumann.
        bool isNeumann() const
        {
            return type_ == Neumann;
        }
	/// @brief Type query.
	/// @return true if the type is Periodic.
        bool isPeriodic() const
        {
            return type_ == Periodic;
        }

    protected: // data
	BCType type_;
	T value_;
    };

    /// @brief Stream insertion for BCBase.
    template<typename charT, class traits, typename T>
    std::basic_ostream<charT,traits>&
    operator<<(std::basic_ostream<charT,traits>& os,
               const BCBase<T>& bc)
    {
        bc.write(os);
        return os;
    }




    /// @brief A class for representing a flow boundary condition.
    class FlowBC : public BCBase<double>
    {
    public:
	/// @brief Default constructor, that makes a noflow condition (Neumann, value 0.0).
        FlowBC()
            : BCBase<double>(Neumann, 0.0)
        {
        }
	/// @brief Constructor taking a type and value.
	/// @param type the condition type.
	/// @param value the condition value.
        FlowBC(BCType type, double value)
            : BCBase<double>(type, value)
        {
	    assert(isNeumann() || isDirichlet() || isPeriodic());
        }

	/// @brief Forwarding the relevant type queries.
	using BCBase<double>::isDirichlet;
	using BCBase<double>::isNeumann;
	using BCBase<double>::isPeriodic;

	/// @brief Query a Dirichlet condition.
	/// @return the pressure condition value
        double pressure() const
        {
#ifndef NDEBUG
            OPM_ERROR_IF(!isDirichlet(), "Pressure boundary conditions are only valid for Dirichlet boundaries");
#endif
            return value_;
        }
	/// @brief Query a Neumann condition.
	/// @return the outwards flux condition value.
        double outflux() const
        {
#ifndef NDEBUG
            OPM_ERROR_IF(!isNeumann(), "Outflux boundary conditions are only valid for Neumann boundaries");
#endif
            return value_;
        }
	/// @brief Query a Periodic condition.
	/// @return the pressure difference condition value.
        double pressureDifference() const
        {
#ifndef NDEBUG
            OPM_ERROR_IF(!isPeriodic(), "Pressure difference boundary conditions are only valid for periodic boundaries");
#endif
            return value_;
        }
    };



    /// @brief A class for representing a saturation boundary condition.
    class SatBC : public BCBase<double>
    {
    public:
	/// @brief Default constructor, that makes a Dirichlet condition with value 1.0.
        SatBC()
            : BCBase<double>(Dirichlet, 1.0)
        {
        }
	/// @brief Constructor taking a type and value.
	/// @param type the condition type.
	/// @param value the condition value.
        SatBC(BCType type, double value)
            : BCBase<double>(type, value)
        {
	    assert(isDirichlet() || isPeriodic());
        }
	/// @brief Forwarding the relevant type queries.
	using BCBase<double>::isDirichlet;
	using BCBase<double>::isPeriodic;

	/// @brief Query a Dirichlet condition.
	/// @return the boundary saturation value
        double saturation() const
        {
            assert (isDirichlet());
            return value_;
        }
	/// @brief Query a Periodic condition.
	/// @return the saturation difference value.
        double saturationDifference() const
        {
            assert (isPeriodic());
            return value_;
        }
    };


    /// @brief A class for representing a surface volume boundary condition.
    template <int numComponents>
    class SurfvolBC : public BCBase<Dune::FieldVector<double, numComponents> >
    {
    public:
        typedef Dune::FieldVector<double, numComponents> CompVec;
        typedef BCBase<CompVec> Base;
	/// @brief Default constructor, that makes a Dirichlet condition with all zero component values.
        SurfvolBC()
            : Base(Base::Dirichlet, CompVec(0.0))
        {
        }
	/// @brief Constructor taking a value.
	/// @param value the condition value.
        explicit SurfvolBC(Dune::FieldVector<double, numComponents> value)
            : Base(Base::Dirichlet, value)
        {
	    assert(isDirichlet());
        }
	/// @brief Forwarding the relevant type query.
	using Base::isDirichlet;

	/// @brief Query a Dirichlet condition.
	/// @return the boundary saturation value
        CompVec surfvol() const
        {
            return Base::value_;
        }
    };



    class PeriodicConditionHandler
    {
    public:
        PeriodicConditionHandler()
        {
        }

        explicit PeriodicConditionHandler(int num_different_boundary_ids)
	    : periodic_partner_bid_(num_different_boundary_ids, 0),
	      canonical_bid_(num_different_boundary_ids, 0)
        {
        }

        void resize(int new_size)
        {
            periodic_partner_bid_.resize(new_size, 0);
            canonical_bid_.resize(new_size, 0);
        }

        bool empty() const
        {
            return periodic_partner_bid_.empty() && canonical_bid_.empty();
        }

        void clear()
        {
            periodic_partner_bid_.clear();
	    canonical_bid_.clear();
        }

        int size() const
        {
            return periodic_partner_bid_.size();
        }

        void setPeriodicPartners(int boundary_id_1, int boundary_id_2)
        {
            assert(boundary_id_1 >= 0 && boundary_id_1 < int(periodic_partner_bid_.size()));
            assert(boundary_id_2 >= 0 && boundary_id_2 < int(periodic_partner_bid_.size()));
            assert(periodic_partner_bid_[boundary_id_1] == 0);
            assert(periodic_partner_bid_[boundary_id_2] == 0);
            periodic_partner_bid_[boundary_id_1] = boundary_id_2;
            periodic_partner_bid_[boundary_id_2] = boundary_id_1;
        }

        int getPeriodicPartner(int boundary_id) const
        {
            assert(boundary_id >= 0 && boundary_id < int(periodic_partner_bid_.size()));
            return periodic_partner_bid_[boundary_id];
        }

	void setCanonicalBoundaryId(int boundary_id, int canonical_bid)
	{
	    assert(boundary_id >= 0 && boundary_id < int(canonical_bid_.size()));
            canonical_bid_[boundary_id] = canonical_bid;
	}

        int getCanonicalBoundaryId(int boundary_id) const
        {
            assert(boundary_id >= 0 && boundary_id < int(canonical_bid_.size()));
            return canonical_bid_[boundary_id];
        }

        template<typename charT, class traits>
        void write(std::basic_ostream<charT,traits>& os) const
        {
            for (int i = 0;  i < int(periodic_partner_bid_.size()); ++i) {
                os << "Partner of bid " << i << " is " << periodic_partner_bid_[i]
		   << "  (canonical bid is " << canonical_bid_[i] << '\n';
            }
            os << std::endl;
        }
    private:
        std::vector<int> periodic_partner_bid_;
	std::vector<int> canonical_bid_;
    };

    template<typename charT, class traits>
    std::basic_ostream<charT,traits>&
    operator<<(std::basic_ostream<charT,traits>& os,
               const PeriodicConditionHandler& pch)
    {
        pch.write(os);
        return os;
    }


    template <typename T>
    class DummyVec
    {
    public:
	DummyVec() {}
    explicit DummyVec(int) {}
	void resize(int) {}
	void clear() {}
    };

    template <bool FC = false, bool SC = false, bool ZC = false, int numComponents = 3>
    class BasicBoundaryConditions : public PeriodicConditionHandler,
        private std::conditional<FC, std::vector<FlowBC>, DummyVec<FlowBC> >::type,
        private std::conditional<SC, std::vector<SatBC>,  DummyVec<SatBC>  >::type,
        private std::conditional<ZC, std::vector<SurfvolBC<numComponents> >, DummyVec<SurfvolBC<numComponents> > >::type
    {
    public:
	typedef typename std::conditional<FC, std::vector<FlowBC>, DummyVec<FlowBC> >::type FlowConds;
	typedef typename std::conditional<SC,  std::vector<SatBC>, DummyVec<SatBC>  >::type SatConds;
	typedef typename std::conditional<ZC, std::vector<SurfvolBC<numComponents> >, DummyVec<SurfvolBC<numComponents> > >::type SurfvolConds;
	const static bool HasFlowConds = FC;
	const static bool HasSatConds = SC;
	const static bool HasSurfvolConds = SC;

        BasicBoundaryConditions()
	{
        }

        explicit BasicBoundaryConditions(int num_different_boundary_ids)
	    : PeriodicConditionHandler(num_different_boundary_ids),
	      FlowConds(num_different_boundary_ids),
	      SatConds(num_different_boundary_ids),
 	      SurfvolConds(num_different_boundary_ids)
	{
        }

        void resize(int new_size)
        {
	    PeriodicConditionHandler::resize(new_size);
	    FlowConds::resize(new_size);
	    SatConds::resize(new_size);
            SurfvolConds::resize(new_size);
        }

        bool empty() const
        {
            return PeriodicConditionHandler::empty();
        }

        void clear()
        {
	    PeriodicConditionHandler::clear();
	    FlowConds::clear();
	    SatConds::clear();
	    SurfvolConds::clear();
        }

        int size() const
        {
            return PeriodicConditionHandler::size();
        }

	FlowBC& flowCond(int i)
	{
	    return FlowConds::operator[](i);
	}

	const FlowBC& flowCond(int i) const
	{
	    return FlowConds::operator[](i);
	}

        template <class BoundaryFace>
	const FlowBC& flowCond(const BoundaryFace& bf) const
	{
            assert(bf.boundary());
	    return FlowConds::operator[](bf.boundaryId());
	}

	SatBC& satCond(int i)
	{
	    return SatConds::operator[](i);
	}

	const SatBC& satCond(int i) const
	{
	    return SatConds::operator[](i);
	}

        template <class BoundaryFace>
	const SatBC& satCond(const BoundaryFace& bf) const
	{
            assert(bf.boundary());
	    return SatConds::operator[](bf.boundaryId());
	}

	SurfvolBC<numComponents>& surfvolCond(int i)
	{
	    return SurfvolConds::operator[](i);
	}

	const SurfvolBC<numComponents>& surfvolCond(int i) const
	{
	    return SurfvolConds::operator[](i);
	}

        template <class BoundaryFace>
	const SurfvolBC<numComponents>& surfvolCond(const BoundaryFace& bf) const
	{
            assert(bf.boundary());
	    return SurfvolConds::operator[](bf.boundaryId());
	}

        template<typename charT, class traits>
        void write(std::basic_ostream<charT,traits>& os) const
        {
	    PeriodicConditionHandler::write(os);
	}
    };

    template<typename charT, class traits, bool F, bool S> //, bool P>
    std::basic_ostream<charT,traits>&
    operator<<(std::basic_ostream<charT,traits>& os,
               const BasicBoundaryConditions<F,S>& bcs)
    {
        bcs.write(os);
        return os;
    }


} // namespace Opm


#endif // OPENRS_BOUNDARYCONDITIONS_HEADER
