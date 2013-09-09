//==============================================================================
//!
//! \file material.cpp
//!
//! \date Oct 1 2007
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Material interface.
//!
//==============================================================================

#include "config.h"
#include "material.hh"
#include "materials.hh"

#include <fstream>
#include <iostream>

namespace Opm {
namespace Elasticity {

Material* Material::create (int ID, const Dune::DynamicVector<double>& params)
{
  switch (params.size())
    {
      // Isotropic material, E and nu given
    case 1: return new Isotropic(ID,params[0],double(0));
    case 2: return new Isotropic(ID,params[0],params[1]);
    case 3: return new Isotropic(ID,params[0],params[1],params[2]);

      // Diagonal-orthotropic material, Ex, Ey, Ez, Gxy, Gxz, Gyz given
    case 4: return new OrthotropicD(ID,
                                    params[0],params[1],params[2],
                                    params[3]);
    case 5: return new OrthotropicD(ID,
                                    params[0],params[1],params[2],
                                    params[3],params[4]);
    case 6: return new OrthotropicD(ID,
                                    params[0],params[1],params[2],
                                    params[3],params[4],params[5]);

      // General symmetric orthotropic material, upper triangle of C given
    case 9:
      {
        Dune::DynamicVector<double>mat(21,double(0));
	mat[0]  = params[0]; mat[1] = params[1]; mat[2] = params[2];
	mat[6]  = params[3]; mat[7] = params[4];
	mat[11] = params[5];
	mat[15] = params[6];
	mat[18] = params[7];
	mat[20] = params[8];
	return new OrthotropicSym(ID,mat);
      }
      break;

    case 21: return new OrthotropicSym(ID,params);
    }

  std::cerr <<"Material::create: Invalid number of material parameters, "
            << params.size() << std::endl;
  return 0;
}

Material* Material::create(int ID, const std::string& file)
{
  std::ifstream f;
  f.open(file.c_str());
  if (!f.good()) {
    std::cerr << "Rock file " << file << " not found, using default isotropic material" << std::endl;
    return new Isotropic(ID,100,0.38,3);
  }
  std::string str,str2;
  f >> str;
  if (str == "ti") { // transverse isotropic material with unit axis in z-direction. 
      Dune::DynamicVector<double>mat(21,double(0));
      double a,b,c,d,e,rho;
      f >> a >> b >> c >> d >> e;
      f >> str2;
      if (str2 != "density") {
          std::cerr << "Rock file " << file << " has wrong format. when isotropycase is \"ti\", " << std::endl 
                    << "the next line need 5 a,b,c,d and e for C matrix with C11=C22=a, C33=c, " << std::endl
                    << "C44=C55=d, C66=e, C12=a-2e and C13=b. Then keyword density on next " << std::endl
                    << "line followed by material density (double) on last line." << std::endl;
          exit(1);
          assert(0);
      }
      f >> rho;
      mat[0] = a; mat[1] = a-2.f*e; mat[2] = b; 
      mat[6] = a; mat[7] = b; mat[11] = c;
      mat[15] = d; mat[18] = d; mat[20] = e;
      return new  OrthotropicSym(ID,mat);
  }
  else if (str == "anisotropic") { // full symmetric matrix with 21 elements (upper triangle, indexed c11 c12 ... c16 c22 c23 ... c66)
      Dune::DynamicVector<double>mat(21,double(0));
      double rho;
      for (int cuidx=0; cuidx<21; cuidx++) {
          f >> mat[cuidx];
      } 
      f >> str2;
      if (str2 != "density") {
          std::cerr << "Rock file " << file << " has wrong format. when isotropycase is \"anisotropic\", " << std::endl 
                    << "the next line need 21 entries for upper C matrix. Then keyword density on next " << std::endl
                    << "line followed by material density (double) on last line." << std::endl;
          exit(1);
          assert(0);
      }
      
      f >> rho;
      return new  OrthotropicSym(ID,mat);
  }
  else {
      double p1, p2, rho, E, nu;
      f >> p1 >> p2;
      f >> str2;
      if (str2 != "density") {
          std::cerr << "Rock file " << file << " has wrong format. when isotropycase is either \"km\", " << std::endl
                    << "\"lm\", \"en\" or \"vpvs\", the next line need 21 entries for upper C matrix. " << std::endl
                    << "Then keyword density on next line followed by material density (double) on last line." << std::endl;
          exit(1);
          assert(0);
      }
      f >> rho;
      if (str == "vpvs") {
          p1 = rho*p1*p1-4.f/3*rho*p2*p2;
          p2 = rho*p2*p2;
          str = "km";
      }
      if (str == "km") { // bulk modulus and shear modulus
          nu = (3*p1-2*p2)/(6*p1+2*p2);
          E = 2*p2*(1+nu);
      } else if (str == "lm") { // lame's lambda and shear modulus
          nu = p1/(2*(p2+p1));
          E = (1+nu)*(1-2*nu)*p1/nu;
      } else if (str == "en") { // young's modulus and poisson's ratio
          E = p1;
          nu = p2;
      } else {
        std::cerr << "Could not parse rock file " << file << ", bailing" << std::endl;
        exit(1);
      }

      if (nu < 0 || nu > 0.5) {
        std::cerr << "Material in " << file << " is not isotropic (nu=" << nu << "), bailing" << std::endl;
        exit(1);
      }
      
      return new Isotropic(ID,E,nu,rho);
  }
}

}
}
