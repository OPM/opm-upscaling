//==============================================================================
//!
//! \file material.hh
//!
//! \date Oct 1 2007
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Material interface.
//!
//==============================================================================
#ifndef MATERIAL_HH_
#define MATERIAL_HH_

#include <opm/common/utility/platform_dependent/disable_warnings.h>

#include <dune/common/fmatrix.hh>
#include <dune/common/dynvector.hh>

#include <opm/common/utility/platform_dependent/reenable_warnings.h>

namespace Opm {
namespace Elasticity {


/*!
  \brief This is a base class for linear elastic materials.
  \details It is an abstract class since some of the member functions are
  purely virtual.
*/

class Material
{
protected:
  //! \brief Default constructor creating an empty material.
  Material(int ID = 0, double density = 0.0)
    : id(ID), rho(density)
  {
  }

  //! \brief Prints the material properties to a stream.
  virtual std::ostream& write(std::ostream& os) const
  { 
    return os;
  }
public:
  //! \brief Empty virtual destructor.
  virtual ~Material() {}

  //! \brief Returns the external material id.
  int num() const
  {
    return id;
  }

  //! \brief Returns the number of parameters describing this material.
  virtual int numPar() const = 0;
  //! \brief Returns the \a ipar'th parameter describing this material.
  virtual double getPar(int ipar = 1) const
  {
    return double(0);
  }

  //! \brief Establishes the full constitutive matrix for this material.
  //! \param[out] C The constitutive matrix
  //! \param[in] invers If \e true, set up the inverse matrix instead
  virtual bool getConstitutiveMatrix(Dune::FieldMatrix<double,6,6>& C,
				     bool invers = false) const = 0;

  //! \brief Establishes the full constitutive matrix for this material.
  //! \param[out] C The constitutive matrix
  //! \param[in] invers If \e true, set up the inverse matrix instead
  virtual bool getConstitutiveMatrix(Dune::FieldMatrix<double,3,3>& C,
				     bool invers = false) const = 0;

  //! \brief Returns the mass density of this material.
  double getMassDensity() const
  {
    return rho;
  }

  //! \brief Global stream operator printing a material properties object.
  friend std::ostream& operator<<(std::ostream& os, const Material& m)
  {
    return m.write(os);
  }

  //! \brief Creates a material object of a given type.
  //! \details The material type depends on the number of parameters provided.
  //! \param[in] ID External number for this material
  //! \param[in] params Array of material parameters
  static Material* create(int ID, const Dune::DynamicVector<double>& params);

  //! \brief Creates a material object from a rocklist
  //! \param[in] ID ID of the material
  //! \param[in] file The URL to the rocklist
  static Material* create(int ID, const std::string& file);
private:
  int  id;  //!< External material number
  double rho; //!< Mass density
};

}
}

#endif
