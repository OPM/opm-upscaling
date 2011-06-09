/*
  Copyright 2011 - Magne Haveraaen, Helmer André Friis and Hans Munthe-Kaas.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/


//---------------------- Include Files ----------------------------------------
#include <iostream>
//#include <strstream>
#include <sstream>
#include <fstream>
#include <iomanip>

#include <cstdlib>
#include <iostream>

#include "MeshSFSD.h"

using namespace std;  // This leads to that all names in the namespace std may be 
                      // used without the qualifier "std::".




//---------------------  Constants --------------------------------------------


//---------------------- Types ------------------------------------------------






//---------------------- Public Functions -------------------------------------


Boolean MeshSFSD::u_DI_0u() const 
{
  Boolean is_OK = ((m_shape . u_domainDimension_0u() == m_domainDimension) && (m_p != NULL));
  if (((m_domainDimension < 0) || (m_domainDimension > 2)))
    {
      (is_OK = false);
    }
  else
    if ((m_domainDimension == 0))
      {
        (is_OK = false);
      }
    else
      if ((m_domainDimension == 1))
        {
          {
            int a_0;
            m_shape . u_getElement_0g(0, a_0);
            (is_OK = (((m_nRows == 1) && (m_nCols == m_size)) && (m_size == a_0)));
          }
        }
      else
        if ((m_domainDimension == 2))
          {
            {
              int b_0;
              m_shape . u_getElement_0g(0, b_0);
              {
                int c_0;
                m_shape . u_getElement_0g(1, c_0);
                (is_OK = ((((m_nRows == m_size) && (m_nCols == m_size)) && (m_size == b_0)) && (m_size == c_0)));
              }
            }
          }
  return is_OK;
}

MeshSFSD::MeshSFSD()  : m_shape(MeshShape()), m_domainDimension(0), m_size(0), m_nRows(0), m_nCols(0), m_p(NULL), m_isMatBanded(false), m_numbSubDiag(0), m_numbSuperDiag(0) 
{
}

MeshSFSD::MeshSFSD(const MeshSFSD &V)  : m_shape(V . m_shape), m_domainDimension(V . m_domainDimension), m_size(V . m_size),m_isMatBanded(V.m_isMatBanded), m_numbSubDiag(V.m_numbSubDiag), m_numbSuperDiag(V.m_numbSuperDiag) 
{
  int r;
  int c;
  double *smp;
  double *dmp;
  (m_nRows = ((V . m_domainDimension == 1) ? 1 : V . m_size));
  (m_nCols = V . m_nCols);
  (m_p = new double[(m_nRows * m_nCols)]);
  //TEST_ASSERT(m_p != NULL,"Memory allocation failed",FATAL);
  for ((r = 0);(r < m_nRows); (r ++))
    {
      for ((c = 0);(c < m_nCols); (c ++))
        {
          (dmp = ul_GetItemPtr(m_p, r, c));
          (smp = ul_GetItemPtr(V . m_p, r, c));
          ((* dmp) = (* smp));
        }
    }
}

MeshSFSD::MeshSFSD(const MeshShape &S)  : m_shape(S), m_domainDimension(S.u_domainDimension_0u())
{
  int r;
  int c;
  double *dmp;

  m_isMatBanded = false; 
  m_numbSubDiag = 0; 
  m_numbSuperDiag = 0;
  
  S.u_getElement_0g(0, m_size);
  
  (m_nRows = ((m_domainDimension == 1) ? 1 : m_size));
  (m_nCols = m_size);
  (m_p = new double[(m_nRows * m_nCols)]);
  //TEST_ASSERT(m_p != NULL,"Memory allocation failed",FATAL);
  for ((r = 0);(r < m_nRows); (r ++))
    {
      for ((c = 0);(c < m_nCols); (c ++))
        {
          (dmp = ul_GetItemPtr(m_p, r, c));
          ((* dmp) = 0.0);
        }
    }
}

MeshSFSD::MeshSFSD(const MeshShape &S, const double &v)  : m_shape(S), m_domainDimension(S.u_domainDimension_0u())
{
  int r;
  int c;
  double *dmp;

  m_isMatBanded = false; 
  m_numbSubDiag = 0; 
  m_numbSuperDiag = 0;
  
  S.u_getElement_0g(0, m_size);

  (m_nRows = ((m_domainDimension == 1) ? 1 : m_size));
  (m_nCols = m_size);
  (m_p = new double[(m_nRows * m_nCols)]);
  //TEST_ASSERT(m_p != NULL,"Memory allocation failed",FATAL);
  for ((r = 0);(r < m_nRows); (r ++))
    {
      for ((c = 0);(c < m_nCols); (c ++))
        {
          (dmp = ul_GetItemPtr(m_p, r, c));
          ((* dmp) = v);
        }
    }
}

void MeshSFSD::u_copy_1u() 
{
  if (m_p)
    {
      delete[] m_p;
      (m_p = NULL);
    }
  m_shape . u_copy_1u();
  (m_domainDimension = 0);
  (m_size = 0);
  (m_nRows = 0);
  (m_nCols = 0);

  m_isMatBanded = false; 
  m_numbSubDiag = 0; 
  m_numbSuperDiag = 0;

}

void MeshSFSD::u_copy_1u(const MeshShape &S) 
{
  int oldDomDim = m_domainDimension;
  int oldSize = m_size;
  (m_shape = S);
  (m_domainDimension = S . u_domainDimension_0u());
  S . u_getElement_0g(0, m_size);

  (m_nRows = ((m_domainDimension == 1) ? 1 : m_size));
  (m_nCols = m_size);
  if (((m_size != oldSize) || (m_domainDimension != oldDomDim)))
    {
      ul_Allocate();
    }
  ul_Preset(0.0);

  m_isMatBanded = false; 
  m_numbSubDiag = 0; 
  m_numbSuperDiag = 0;

}

void MeshSFSD::u_copy_1u(const MeshShape &S, const double &v) 
{
  int oldDomDim = m_domainDimension;
  int oldSize = m_size;
  (m_shape = S);
  (m_domainDimension = S . u_domainDimension_0u());
  S . u_getElement_0g(0, m_size);

  (m_nRows = ((m_domainDimension == 1) ? 1 : m_size));
  (m_nCols = m_size);
  if (((m_size != oldSize) || (m_domainDimension != oldDomDim)))
    {
      ul_Allocate();
    }
  ul_Preset(v);

  m_isMatBanded = false; 
  m_numbSubDiag = 0; 
  m_numbSuperDiag = 0;

}

MeshShape MeshSFSD::u_getShape_0u() const 
{
  return m_shape;
}


int MeshSFSD::u_getNumberOfRows_0u() const 
{
  return m_nRows;
}


int MeshSFSD::u_getNumberOfColumns_0u() const 
{
  return m_nCols;
}


Boolean MeshSFSD::u_sameShape_0u(const MeshSFSD &CP) const 
{
  return (m_shape == CP . m_shape);
}

Boolean MeshSFSD::u_legalPoint_0u(const MeshPoint &P) const 
{
  return (m_shape == P . u_getShape_0u());
}

double MeshSFSD::u_getData_0u(const MeshPoint &P) const 
{
  int r;
  int c;
  double *dmp;
  {
    int d_0;
    P . u_getElement_0g(0, d_0);
    (r = ((m_domainDimension == 1) ? 0 : d_0));
    ;
    P . u_getElement_0g(0, d_0);
    {
      int g_0;
      P . u_getElement_0g(1, g_0);
      (c = ((m_domainDimension == 1) ? d_0 : g_0));
    }
    (dmp = ul_GetItemPtr(m_p, r, c));
    //    cout << endl;
    //cout << "r= " << r << " c= " << c << endl;
    return double((* dmp));
  }
}



void MeshSFSD::u_print_0u() const 
{

  int r;
  int c;
  // double *dmp;

  if (m_domainDimension == 1)
  {
    //We print the vector
    cout << "Vector:" << endl;
    for (int i=0; i < m_size; i++)
    {
      cout << "  " << m_p[i] << endl;
    }
  }
  else if (m_domainDimension == 2)
  {
    //We print the matrix
    cout << "Matrix:" << endl;
    for ((r = 0);(r < m_size); (r ++))
    {
      for ((c = 0);(c < m_size); (c ++))
      {
	//dmp = ul_GetItemPtr(m_p, r, c);
	//	cout << "  " << *m_p;
	cout << "  " << m_p[(r * m_nCols) + c];
      }
      cout << endl;
    }
  }

 
}



int MeshSFSD::u_domainDimension_0u() const 
{
  return m_domainDimension;
}




double MeshSFSD::u_getElementFromMatMultOfNonQuadraticMatrices_0u(const int& indRowObj, const MeshSFSD& MatB, const int& indColB) const 
{

  //A precondition:
  if (m_nCols != MatB.m_nRows)
  {
    cout << "Wrong input in matrix multiplication of (possibly) non-quadratic matrices.  Program is stopped." << endl;
    exit(0);
  }

  double ret = 0.0;

  double *dmp1, *dmp2;

  //  cout << "m_nRows= " << m_nRows << " m_nCols= " << m_nCols << " MatB.m_nRows= " << MatB.m_nRows << " MatB.m_nCols= " << MatB.m_nCols << endl;

  for (int k=0; k < m_nCols; k++)
  {
    dmp1 = ul_GetItemPtr(m_p, indRowObj, k);
    dmp2 = ul_GetItemPtrGeneral(MatB.m_p, k, indColB, MatB.m_nCols);//Test 11.12.2006
    //dmp2 = ul_GetItemPtr(MatB.m_p, indColB, k);//Test 11.12.2006
    //cout << "dmp1= " << (*dmp1) << " dmp2= " << (*dmp2) << endl;
    ret += (*dmp1)*(*dmp2);
  }

  return ret;
}




void MeshSFSD::u_setupForNonQuadraticMatrices_1u(const int& numbOfRows, const int& numbOfColumns, const double &v) 
{

  //******************************************************************************
  //NB: This routine is only for matrices; NOT for vectors !!!
  //
  //NOTE: A given matrix is also allowed to be quadratic here, and a matrix with
  //only 1 column may be represented as well also in this case with
  //m_domainDimension = 2.
  //******************************************************************************
  m_isMatBanded = false;
  m_nRows = numbOfRows;
  m_nCols = numbOfColumns;
  m_domainDimension = 2;
  m_size = m_nRows;

  int vHelp[MAXDIRECTIONSIZE];
  vHelp[0] = m_nRows;
  vHelp[1] = m_nCols;
  m_shape.u_setShape_1u(2, vHelp);


  if (m_p)
    {
      delete[] m_p;
      (m_p = NULL);
    }
  m_p = new double[(m_nRows * m_nCols)];
  //TEST_ASSERT(m_p != NULL,"Memory allocation failed",FATAL);

  ul_Preset(v);

  //  cout << "MARTA!!" << "m_nRows= " << m_nRows << " m_nCols= " << m_nCols << endl;

}




void MeshSFSD::u_setupForBandedMatrices_1u(const MeshShape &S, const double &v, const int& numbSubDiag, const int& numbSuperDiag) 
{

  //******************************************************************************
  //NB: This routine is only for matrices; NOT for vectors !!!
  //The storage pattern follows the LAPACK convention (but it is transposed since LAPACK is in FORTRAN)
  //******************************************************************************
  m_isMatBanded = true;
  m_numbSubDiag = numbSubDiag;
  m_numbSuperDiag = numbSuperDiag;
  m_domainDimension = 2;
  m_shape = S;
  S.u_getElement_0g(0, m_size);

  //JUST TO BE SURE:
  if (m_shape.u_domainDimension_0u() != m_domainDimension)
  {
    cout << "Input MeshShape is WRONG. This is only for matrices !!! Program is STOPPED." << endl;
    exit(0);
  }

  m_nRows = m_size;
  m_nCols = 2*m_numbSubDiag + m_numbSuperDiag + 1;//NB: This is due to the LAPACK convention. WRONG ???

  //  m_nCols = m_size;
  //m_nRows = 2*m_numbSubDiag + m_numbSuperDiag + 1;//NB: This is due to the LAPACK convention. WRONG ???

  if (m_p)
    {
      delete[] m_p;
      (m_p = NULL);
    }
  m_p = new double[(m_nRows * m_nCols)];
  //TEST_ASSERT(m_p != NULL,"Memory allocation failed",FATAL);

  ul_Preset(v);

}



void MeshSFSD::u_setData_1u(const MeshPoint &P, const double &v) 
{

  double *dmp;
  if (!m_isMatBanded)
  {
    int r;
    int c;
    int i_0;
    P . u_getElement_0g(0, i_0);
    (r = ((m_domainDimension == 1) ? 0 : i_0));
    ;
    //    P . u_getElement_0g(0, i_0);
    {
      int l_0;
      P . u_getElement_0g(1, l_0);
      (c = ((m_domainDimension == 1) ? i_0 : l_0));
    }
    (dmp = ul_GetItemPtr(m_p, r, c));
    ((* dmp) = v);
  }
  else
  {
    //The following code is only for banded matrices and exploits the LAPACK convention.
    int i = P[0];
    int j = P[1];
    //if ((i==1) && (j==6))
    // {
      //      cout << "val= " << v << " pos= " << ((m_numbSubDiag + m_numbSuperDiag + i - j)*m_nCols + j) << endl;
      //cout << "val= " << v << " pos= " << (j*m_nCols) + (m_numbSubDiag + m_numbSuperDiag + i - j) << endl;
    // }
    //    dmp = ul_GetItemPtr(m_p, m_numbSubDiag + m_numbSuperDiag + 1 + i - j, j); //OK, eller skal det transponeres ???
    //    dmp = ul_GetItemPtr(m_p, j, m_numbSubDiag + m_numbSuperDiag + 1 + i - j); //OK ???
    //    dmp = ul_GetItemPtrBanded(m_p, j, i); //OK ???
    //    dmp = ul_GetItemPtrBanded(m_p, i, j); //OK ???

    //    cout << "i= " << i << " j= " << j << " val= " << v << "--------------" << " dr= " << j << " dc= " << m_numbSubDiag + m_numbSuperDiag + i - j << " 1DPos= " << (j*m_nCols) + (m_numbSubDiag + m_numbSuperDiag + i - j) << endl;

    dmp = ul_GetItemPtrBanded(m_p, i, j); //OK ???
    (* dmp) = v;
  }
}


void MeshSFSD::u_mMatMult_2u(MeshSFSD &M) const 
{
  int r;
  int c;
  int k;
  double *dmp;
  double *smp1;
  double *smp2;
  double *temp = new double[(m_nRows * m_nCols)];
  //TEST_ASSERT(temp != NULL,"Memory allocation failed",FATAL);

  //****************************************************
  //NB: It seems like this routine assumes that we are
  //multiplying two quadratic matrices of the SAME size !!!
  //****************************************************

  for ((r = 0);(r < m_size); (r ++))
    {
      for ((c = 0);(c < m_size); (c ++))
        {
          (dmp = ul_GetItemPtr(temp, r, c));
          ((* dmp) = 0.0);
          for ((k = 0);(k < m_size); (k ++))
            {
              (dmp = ul_GetItemPtr(temp, r, c));
              (smp1 = ul_GetItemPtr(m_p, r, k));
              (smp2 = ul_GetItemPtr(M . m_p, k, c));
              ((* dmp) = ((* dmp) + ((* smp1) * (* smp2))));
            }
        }
    }
  for ((r = 0);(r < m_size); (r ++))
    {
      for ((c = 0);(c < m_size); (c ++))
        {
          (dmp = ul_GetItemPtr(M . m_p, r, c));
          (smp1 = ul_GetItemPtr(temp, r, c));
          ((* dmp) = (* smp1));
        }
    }
  delete[] temp;
}

void MeshSFSD::u_mVecMult_2u(MeshSFSD &V) const 
{
  int i;
  int k;
  double *smp1;
  double *temp = new double[m_size];
  //TEST_ASSERT(temp != NULL,"Memory allocation failed",FATAL);

  //****************************************************
  //NB: It seems like this routine assumes that we are
  //multiplying a quadratic matrix with a vector
  //having the SAME number of rows !!!
  //****************************************************

  for ((i = 0);(i < m_size); (i ++))
    {
      (temp[i] = 0.0);
    }
  for ((i = 0);(i < m_size); (i ++))
    {
      for ((k = 0);(k < m_size); (k ++))
        {
          (smp1 = ul_GetItemPtr(m_p, i, k));
          (temp[i] = (temp[i] + ((* smp1) * V . m_p[k])));
        }
    }
  for ((i = 0);(i < m_size); (i ++))
    {
      (V . m_p[i] = temp[i]);
    }
  delete[] temp;
}

void MeshSFSD::u_mdivide_1u(const MeshSFSD &V) 
{
  int nors = 1;
  int r;
  int c;
  double *dmp;
  double *smp;
  double *mat_trans = new double[(V . m_nRows * V . m_nCols)];
  //TEST_ASSERT(mat_trans != NULL,"Memory allocation failed",FATAL);
  for ((r = 0);(r < V . m_nRows); (r ++))
    {
      for ((c = 0);(c < V . m_nCols); (c ++))
        {
          (smp = ul_GetItemPtr(V . m_p, r, c));
          (dmp = ul_GetItemPtr(mat_trans, c, r));
          ((* dmp) = (* smp));
        }
    }
  eqsolve_(mat_trans, m_p, (& m_size), (& nors));
  delete[] mat_trans;
}



void MeshSFSD::u_mdivideBanded_1u(MeshSFSD &V) 
{
  int nors = 1;
//   int r;
//   int c;
//   double *dmp;
//   double *smp;

  int nCols = V.m_nCols;
  int numbSubDiag = V.m_numbSubDiag;
  int numbSuperDiag = V.m_numbSuperDiag;
    eqsolvb_(V.m_p, m_p, (& m_size), (& nors), (& numbSubDiag), (& numbSuperDiag), (& nCols));//Just a test...
}



void MeshSFSD::u_invert_1u() 
{
  matinv_(m_p, (& m_size));
}



void MeshSFSD::operator+=(const MeshSFSD &CP) 
{

  //A precondition:
  if ((m_nRows != CP.m_nRows) || (m_nCols != CP.m_nCols))
  {
    cout << "m_nRows= " << m_nRows << " m_nCols= " << m_nCols << endl;
    cout << "CP.m_nRows= " << CP.m_nRows << " CP.m_nCols= " << CP.m_nCols << endl;
    cout << "Wrong input in matrix addition of (possibly) non-quadratic matrices.  Program is stopped." << endl;
    exit(0);
  }

  double *dmp1, *dmp2;

  for (int i=0; i < m_nRows; i++)
  {
    for (int j=0; j < m_nCols; j++)
    {
      dmp1 = ul_GetItemPtr(m_p, i, j);
      dmp2 = ul_GetItemPtr(CP.m_p, i, j);
      //cout << "dmp1= " << (*dmp1) << " dmp2= " << (*dmp2) << endl;
      (*dmp1) = (*dmp1) + (*dmp2);
    }
  }
  

}




void MeshSFSD::operator*=(const double &v) 
{

  double *dmp1;

  for (int i=0; i < m_nRows; i++)
  {
    for (int j=0; j < m_nCols; j++)
    {
      dmp1 = ul_GetItemPtr(m_p, i, j);
      (*dmp1) = (*dmp1)*v;
    }
  }
  

}


Boolean MeshSFSD::operator==(const MeshSFSD &CP) const 
{
  int i;
  int nItems = (m_nRows * m_nCols);
  Boolean t_OK = true;
  (t_OK = (((((m_shape == CP . m_shape) && (m_domainDimension == CP . m_domainDimension)) && (m_size == CP . m_size)) && (m_nRows == CP . m_nRows)) && (m_nCols == CP . m_nCols)));
  for ((i = 0);((i < nItems) && t_OK); (i ++))
    {
      (t_OK = (m_p[i] == CP . m_p[i]));
    }
  return t_OK;
}

void MeshSFSD::operator=(const MeshSFSD &CP) 
{
  int oldDomDim = m_domainDimension;
  int oldSize = m_size;
  int nItems = (CP . m_nRows * CP . m_nCols);
  (m_shape = CP . m_shape);
  (m_domainDimension = CP . m_domainDimension);
  (m_size = CP . m_size);

  m_nRows = ((m_domainDimension == 1) ? 1 : m_size);
  m_nCols = CP.m_nCols;
  if (((m_size != oldSize) || (m_domainDimension != oldDomDim)))
    {
      ul_Allocate();
    }
  for (int i = 0;(i < nItems); (i ++))
    {
      (m_p[i] = CP . m_p[i]);
    }

  m_isMatBanded = CP.m_isMatBanded;
  m_numbSubDiag = CP.m_numbSubDiag;
  m_numbSuperDiag = CP.m_numbSuperDiag;

}

MeshSFSD::~MeshSFSD() 
{
  if (m_p)
    {
      delete[] m_p;
    }
}









//#ifdef COMPLETE_SOPHUS
//---------------------- Removed from -----------------------------------------







//---------------------- Removed to -------------------------------------------
//#endif //COMPLETE_SOPHUS







double *MeshSFSD::ul_GetItemPtr(double *p, const int r, const int c) const 
{
  return ((p + (r * m_nCols)) + c);
}


double *MeshSFSD::ul_GetItemPtrGeneral(double *p, const int r, const int c, const int StepLength) const 
{
  //

  return ((p + (r * StepLength)) + c);
}



double *MeshSFSD::ul_GetItemPtrBanded(double *p, const int r, const int c) const 
{

  int minPos = 0;
  int maxPos = m_nRows*m_nCols - 1;
  int pos = (c*m_nCols) + (m_numbSubDiag + m_numbSuperDiag + r - c);
  if ((pos < minPos) || (pos > maxPos))
  {
    cout << "WARNING: TROUBLE !!!" << endl;
    cout << "maxPos= " << maxPos << " m_nRows= " << m_nRows << " m_nCols= " << m_nCols << endl;
    cout << "r= " << r << " c= " << c << " pos= " << pos << endl;
  }

  return ((p + (c*m_nCols)) + (m_numbSubDiag + m_numbSuperDiag + r - c));

}



void MeshSFSD::ul_Allocate() 
{
  if (m_p)
    {
      delete[] m_p;
      (m_p = NULL);
    }

  (m_p = new double[(m_nRows * m_nCols)]);
  //TEST_ASSERT(m_p != NULL,"Memory allocation failed",FATAL);
}

void MeshSFSD::ul_Preset(const double initVal) 
{
  int r;
  int c;
  double *dmp;
  for ((r = 0);(r < m_nRows); (r ++))
    {
      for ((c = 0);(c < m_nCols); (c ++))
        {
          (dmp = ul_GetItemPtr(m_p, r, c));
          ((* dmp) = initVal);
        }
    }
}










//---------------------- Protected Functions ----------------------------------



//---------------------- Private Functions ------------------------------------








/*
===============================================================================
    ------------------------ CLASS DESCRIPTION --------------------------
===============================================================================

  DESCRIPTION:
  ------------------
    <Description applicable regarding implementation, modifications and
     corrections etc>



  ENHANCEMENTS:
  -------------

    
    
  MISCELLANEOUS:
  --------------
  

===============================================================================
    --------------------- END OF CLASS DESCRIPTION ----------------------
===============================================================================
*/

