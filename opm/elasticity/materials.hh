//==============================================================================
//!
//! \file materials.hh
//!
//! \date Oct 1 2007
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Material properties.
//!
//==============================================================================
#ifndef MATERIALS_HH_
#define MATERIALS_HH_

#include "material.hh"


namespace Opm {
namespace Elasticity {

/*!
  \brief Isotropic linear elastic material.
*/

class Isotropic : public Material
{
public:
  //! \brief Constructor creating a new isotropic material.
  //! \param[in] ID      External material number
  //! \param[in] Emod    Young's modulus
  //! \param[in] Poisson Poisson's ratio
  //! \param[in] rho_    Mass density
  Isotropic(int ID, double Emod, double Poisson, double rho_= 0.0)
    : Material(ID,rho_)
  {
    E = Emod;
    nu = Poisson;
  }

  //! \brief Empty virtual destructor.
  virtual ~Isotropic() {}

  //! \brief Returns the number of parameters describing this material.
  virtual int numPar() const
  {
    return 2;
  }

  //! \brief Returns the \a ipar'th parameter describing this material.
  virtual double getPar(int ipar = 1) const
  {
    return ipar == 1 ? E : nu;
  }

  //! \brief Set the E modulus of the material
  //! \param[in] E_ The E modulus to be set
  void setE(double E_)
  {
    E = E_;
  }

  //! \brief Returns the E modulus of the material
  double getE() const
  {
    return E;
  }

  //! \brief Establishes the full constitutive matrix for this material.
  //! \param[out] C The constitutive matrix
  //! \param[in] invers If \e true, set up the inverse matrix instead
  virtual bool getConstitutiveMatrix(Dune::FieldMatrix<double,6,6>& C,
                                     bool invers = false) const;

  //! \brief Establishes the full constitutive matrix for this material.
  //! \param[out] C The constitutive matrix
  //! \param[in] invers If \e true, set up the inverse matrix instead
  virtual bool getConstitutiveMatrix(Dune::FieldMatrix<double,3,3>& C,
                                     bool invers = false) const;
protected:
  //! \brief Prints the material properties to a stream.
  virtual std::ostream& write(std::ostream& os) const;

private:
  double E;  //!< Young's modulus
  double nu; //!< Poisson's ratio
};

//! \brief Orthotropic linear elastic material with diagonal constitutive matrix.
class OrthotropicD : public Material
{
public:
  //! \brief Constructor creating a new material.
  //! \param[in] ID External material number
  //! \param[in] Ex Elasticity modulus in local x-direction
  //! \param[in] Ey Elasticity modulus in local y-direction
  //! \param[in] Ez Elasticity modulus in local z-direction
  //! \param[in] Gxy Shear modulus in the local xy-plane
  //! \param[in] Gxz Shear modulus in the local xz-plane, default = Gxy
  //! \param[in] Gyz Shear modulus in the local yz-plane, default = Gxz
  OrthotropicD(int ID, double Ex, double Ey, double Ez,
               double Gxy, double Gxz = double(-1), double Gyz = double(-1));

  //! \brief Empty virtual destructor.
  virtual ~OrthotropicD() {}

  //! \brief Returns the number of parameters describing this material.
  virtual int numPar() const
  {
    return 6;
  }

  //! \brief Returns the \a ipar'th parameter describing this material.
  virtual double getPar(int ipar = 1) const;

  //! \brief Establishes the full constitutive matrix for this material.
  //! \param[out] C The constitutive matrix
  //! \param[in] invers If \e true, set up the inverse matrix instead
  virtual bool getConstitutiveMatrix(Dune::FieldMatrix<double,6,6>& C,
				     bool invers = false) const;

  //! \brief Establishes the full constitutive matrix for this material.
  //! \param[out] C The constitutive matrix
  //! \param[in] invers If \e true, set up the inverse matrix instead
  virtual bool getConstitutiveMatrix(Dune::FieldMatrix<double,3,3>& C,
				     bool invers = false) const;

protected:
  //! \brief Prints the material properties to a stream.
  virtual std::ostream& write(std::ostream& os) const;

private:
  double E[6]; //!< The diagonal of the constitutive matrix
};

//! \brief Orthotropic linear elastic material with symmetric constitutive matrix.
class OrthotropicSym : public Material
{
public:
  //! \brief Constructor creating a new material.
  //! \param[in] ID External material number
  //! \param[in] Cu Upper triangle of the symmetric constitutive matrix
  OrthotropicSym(int ID, const Dune::DynamicVector<double>& Cu);

  //! \brief Empty virtual destructor.
  virtual ~OrthotropicSym() {}

  //! \brief Returns the number of parameters describing this material.
  virtual int numPar() const
  {
    return 21;
  }

  //! \brief Returns the \a ipar'th parameter describing this material.
  virtual double getPar(int ipar = 1) const;

  //! \brief Establishes the full constitutive matrix for this material.
  //! \param[out] C The constitutive matrix
  //! \param[in] invers If \e true, set up the inverse matrix instead
  virtual bool getConstitutiveMatrix(Dune::FieldMatrix<double,6,6>& C,
				     bool invers = false) const;

  //! \brief Establishes the full constitutive matrix for this material.
  //! \param[out] C The constitutive matrix
  //! \param[in] invers If \e true, set up the inverse matrix instead
  virtual bool getConstitutiveMatrix(Dune::FieldMatrix<double,3,3>& C,
				     bool invers = false) const;
protected:
  //! \brief Prints the material properties to a stream.
  virtual std::ostream& write(std::ostream& os) const;

private:
  double Cupper[21]; //!< Upper triangle of the symmetric constitutive matrix
};

}
}

#endif
