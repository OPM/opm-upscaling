/*
  Copyright 2011 IRIS - International Research Institute of Stavanger.

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


#ifndef IRISDuneGridInterfaceH
#define IRISDuneGridInterfaceH

#include <iostream>
#include <iomanip>

//#include<dune/grid/sgrid.hh>
#include<dune/grid/common/gridinfo.hh>

#include"dune/common/mpihelper.hh" // An initializer of MPI
#include"opm/porsol/blackoil/co2fluid/opm/common/exceptions.hh" // We use exceptions
#include<dune/grid/common/mcmgmapper.hh>
#include<dune/grid/common/entity.hh>
#include<dune/grid/common/entitypointer.hh>

#include<opm/porsol/twophase2/OPMKvasiSophusCode/RnPoint.h>
//#include<NestedArray.hpp>
//#include<DoublyNestedArray.hpp>



//§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§

//---------------------- Include Files ----------------------------------------




//---------------------  Constants --------------------------------------------




//---------------------- typedef's and enumerations ----------------------------



//---------------------- Directives -------------------------------------------
//using namespace std;  // This leads to that all names in the namespace std may be
                      // used without the qualifier "std::".


//§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§



template <class DuneGrid>
class IRISDuneGridInterface
{
  //§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§
  //NOTE: The types VertexPointer, EdgePointer, ElementPointer etc. MUST be defined.
  //Together with their corresponding codimensions (integer types) they essentially
  //represent the way we relate to the actual dune-grid through the interface
  //of this class.
  //Obviously this means that we would like it to be so that the user of this class
  //simply needs to apply (relate to) these few concepts above from the dune-grid;
  //The "real dune-grid coding" is performed in the implementation of this class.
  //§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§

  //NBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNB
  //MEN HVORDAN SKAL TYPENE VertexPointer, EdgePointer, ElementPointer etc. "BEKJENTGJØRES"
  //OVER DE DELENE AV PROGRAMMET HVOR DE MÅ VÆRE KJENT??? typedef kommandoen må selvsagt
  //brukes, men hvordan, og dessuten, må den gjentas mange steder etc.???????????
  //DETTE ER SPESIELLT ET AKTUELLT SPØRSMÅL SIDEN DuneGrid I SEG SELV ER EN template her...
  //
  //An example of such a typedef:
  //          typedef typename DuneGrid::LevelGridView LevelGridView;
  //        typedef typename LevelGridView::template Codim<0>::Iterator ElementLevelIterator;
  //          typedef ElementLevelIterator ElementPointer;
  //        typedef typename LevelGridView::template Codim<1>::Iterator EdgeLevelIterator;
  //          typedef EdgeLevelIterator EdgePointer;
  //        typedef typename LevelGridView::template Codim<2>::Iterator VertexLevelIterator;
  //          typedef VertexLevelIterator VertexPointer;
  //
  // MEN ER NOE SLIKT BRUKBART???, OG HVORDAN KAN DEN GJØRES "ALLMENT KJENT" I PROGRAMMET?????
  //
  //NOTE: Should discuss this with Ove! If it is a problem to operate directly with the types
  //VertexPointer, EdgePointer and ElementPointer below we could of course change the relevant
  // functions to be appropriate function templates instead, if that should make things
  //easier???????????????????????????????????+
  //
  //NBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNB

  public:

  typedef typename DuneGrid::ctype ctype;
  enum { dim = DuneGrid::dimension };

  typedef typename DuneGrid::LevelGridView::template Codim<DuneGrid::dimension>::EntityPointer VertexPointer;
  typedef typename DuneGrid::LevelGridView::template Codim<DuneGrid::dimension-1>::EntityPointer EdgePointer;
  typedef typename DuneGrid::LevelGridView::template Codim<0>::EntityPointer ElementPointer;
  typedef typename DuneGrid::LevelGridView::template Codim<1>::EntityPointer FacePointer;

  typedef typename DuneGrid::LevelGridView::template Codim<0>::Iterator ElementIterator;
  typedef typename DuneGrid::LevelGridView::template Codim<1>::Iterator FaceIterator;
  typedef typename DuneGrid::LevelGridView::template Codim<DuneGrid::dimension-1>::Iterator EdgeIterator;
  typedef typename DuneGrid::LevelGridView::template Codim<DuneGrid::dimension>::Iterator VertexIterator;

  typedef typename DuneGrid::LevelGridView::IntersectionIterator IntersectionIterator;

  //Her er GridType en enum f.eks. enum GridType{Onedimensional,Structured2D,Triangular,Cornerpoint2D,.......};
  enum GridType{_somethingRotten, _Structured2D, _Triangular, _Structured3D};

  // Inntil videre gjør vi slik:
  typedef typename Dune::FieldVector<typename DuneGrid::ctype, DuneGrid::dimensionworld> RnVector;


  struct compareVector
  {
    bool operator() (const RnVector& v0, const RnVector& v1)
    {
      return (v0[0]*v1[1]-v0[1]*v1[0] > 0.0);
    }
  };

  template <class EntityPointer>
  struct compareEntityCCWise
  {
    compareEntityCCWise(const RnVector& origo, const RnVector& interiorPoint, IRISDuneGridInterface* intFace)
    :intFace_(intFace)
    ,pi_(std::acos(-1.0))
    ,origo_(origo)
    {
      RnVector interiorDirection = interiorPoint - origo_;
      thetaBranchCut_ = getTheta(interiorDirection) + pi_; //Assuming regular bnd: branch outside domain for bnd-vtx
      if (thetaBranchCut_ >= 2*pi_) thetaBranchCut_ -= 2*pi_;
    }
    bool operator() (const EntityPointer& entP0, const EntityPointer& entP1)
    {
      RnVector ent0Dir = intFace_->getCentroid_V(entP0) - origo_;
      RnVector ent1Dir = intFace_->getCentroid_V(entP1) - origo_;

      double th0 = getTheta(ent0Dir);
      if (th0 > thetaBranchCut_) th0 -= 2*pi_;
      double th1 = getTheta(ent1Dir);
      if (th1 > thetaBranchCut_) th1 -= 2*pi_;

      return (th0 < th1);
    }
    double getTheta(RnVector v)
    {
      Dune::FieldVector<typename DuneGrid::ctype, 2> v2;
      v2[0] = v[0];
      v2[1] = v[1];
      v2 /= v2.two_norm();
      double theta = std::acos(v2[0]);
      if (v2[1] < 0.0)
      {
        theta += 2*(pi_ - theta);
      }
      return theta;
    }

    IRISDuneGridInterface* intFace_;
    const double pi_;
    RnVector origo_;
    double thetaBranchCut_;
  };



//--------------------------------------------------------------------------
//---------------------- Constructors --------------------------------------
//--------------------------------------------------------------------------
  IRISDuneGridInterface(DuneGrid & grid);

//--------------------------------------------------------------------------
//---------------------- Generators ----------------------------------------
//--------------------------------------------------------------------------
//NOTE: It is assumed that the generators can be tailor made for different problems,
//such as to use no more memory than needed from the data structures below.
//--------------------------------------------------------------------------
  void setupForMPFAFPS_FVM();

//--------------------------------------------------------------------------
//---------------------- Iterator operations--------------------------------
//--------------------------------------------------------------------------
  template <class EntityIterator>
  EntityIterator setEntityPointerToFirst(const int& gridLevel) const;
  template <class EntityIterator>
  EntityIterator entityPointerIsAtEnd(const int& gridLevel) const;

//--------------------------------------------------------------------------
//---------------------- Observers -----------------------------------------
//--------------------------------------------------------------------------

  //"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  //NOTE: Visse operasjoner nedenfor (kanskje de fleste?) burde muligens være overlastet med indekser som argumenter...?
  //Det meste av implementasjonen derimot kunne i såfall være felles og kan selvsagt skrives i lokale funksjoner.
  //"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

  //-----------------------------------------------------------------------------
  //Some more or less "standard" type of operations:
  //-----------------------------------------------------------------------------

  int domainDimension() const;
  int cellDimension() const;
  int cellCount(const int& gridLevel, const int& codimEntity) const;

  template <class EntityPointer>
  RnVector getCentroid_V(const EntityPointer entP) const;

  template <class EntityPointer>
  RnPoint getCentroid(const EntityPointer entP) const {return convertToRnPoint(getCentroid_V(entP));}

  //Ny operasjon i IRISDuneGridInterface. (Bør i første omgang kun (innvendig) implementeres for ElementPointer og codimEntity2==1)
  template <int codimEntity2>
  RnVector getCentroid_V(const ElementPointer ElmP, const int& Ent2RelIdx) const;

  template <int codimEntity2>
  RnPoint getCentroid(const ElementPointer ElmP, const int& Ent2RelIdx) const {return convertToRnPoint(getCentroid_V<codimEntity2>(ElmP,Ent2RelIdx));}


  double getVolume(const ElementPointer entP) const;

  template <class EntityPointer>
  int getGlobalIndexFromEntityPointer(const EntityPointer EntP) const;

  template <class EntityPointer>
  int getNumberOfNeighbours(const EntityPointer EntP, const int& codimNeighbour) const;

  void getNeighbours(const VertexPointer Vtx, std::vector<ElementPointer>& neighbours) const;
  void getNeighbours(const VertexPointer Vtx, std::vector<FacePointer>& neighbours) const;
  void getNeighbours(const ElementPointer Vtx, std::vector<ElementPointer>& neighbours) const;
  void getNeighbours(const ElementPointer Vtx, std::vector<FacePointer>& neighbours) const;

  bool isNeighbour(const VertexPointer Vtx, const FacePointer FaceP) const;
  bool isNeighbour(const ElementPointer ElmP, const FacePointer FaceP) const;

  RnVector getNormalVector_V(const ElementPointer ihatP, const int& EdgeRelIdx) const;
  RnPoint getNormalVector(const ElementPointer ihatP, const int& EdgeRelIdx) const {return convertToRnPoint(getNormalVector_V(ihatP,EdgeRelIdx));}

  GridType getGridType() const {return gridType_; }

  //Denne burde være rett frem for et strukturert grid. For et triangulært grid kan jeg selv bidra med en lite kodesnutt...
  double getMinGridSize() const;


  //-----------------------------------------------------------------------------
  //A group of more topological type of operations:
  //-----------------------------------------------------------------------------

  //A relative index function which delivers the corresponding (global) entity pointer:
  template <class EntityPointer>
  EntityPointer getGlobalEntityPointerFromRelativeIndex(const ElementPointer & ElmP, const int& EntRelIdx) const;

  template <class EntityPointer>
  EntityPointer getGlobalEntityPointerFromRelativeIndex(const VertexPointer & VtxP, const int& EntRelIdx) const;

  //An overloaded version of the above function:
  ElementPointer getGlobalEntityPointerFromRelativeIndex(const VertexPointer & ElmP, const int& FaceRelIdx, const int& codimFace, const int& ElemRelIdx) const;

  template <class EntityPointer>
  EntityPointer getGlobalEntityPointerFromRelativeIndex(const ElementPointer & VtxP, const int& Ent2RelIdx, const int& codimEntity2, const int& Ent3RelIdx) const;

  //Another type of relative index function which delivers a relative index:
  template <class EntityPointer>
  int getDecreasedRelativeIndexFromRelativeIndex(const EntityPointer EntP1, const int& Ent2RelIdx, const int& codimEntity2, const int& Ent3RelIdx, const int& codimEntity3) const;

  int getDecreasedRelativeIndexFromRelativeIndex(const VertexPointer & VtxP, const int&FaceRelVtxIdx, const int& ElemRelFaceIdx) const;

  //This will replace the "SLK-set operation".
  void getIntersectionRelativeEdgeIndexSetRightHandCSSorted(const VertexPointer qhatP, const int& ElementGlobIdx, std::vector<int> & out) const;

/*
  // Vente med denne: HAF 30/11-09
  //This is meant to replace parts of the "findTopPath... operation" (the rest can be found from getIntersectionRelativeEdgeIndexSetRightHandCSSorted).
  template <int SIZE>
  void getIntersectionRelativeElementIndexSetRightHandCSSorted(const VertexPointer qhatP, const int& ElementGlobIdx, int[SIZE]& out) const;//NOTE: SIZE==2 for triangular and quadrilateral grids in 2D.
*/

  //-----------------------------------------------------------------------------
  //A few boundary related operations:
  //-----------------------------------------------------------------------------

  int getNumberOfSubBoundaries() const {return facesAtSubBoundary_.size(); }
  int getNumberOfEdgesAtSubBoundary(const int& subBoundaryNumber) const {return facesAtSubBoundary_[subBoundaryNumber].size(); }

  //Decides whether an entity resides on the boundary (probably does NOT need an impl. for elements)
  bool boundaryIndex(const VertexPointer VtxP) const {return vertexAtBoundary_[mapperV_.map(*VtxP)]; }
  bool boundaryIndex(const FacePointer FaceP) const {return faceAtBoundary_[mapperFa_.map(*FaceP)]; }
  bool boundaryIndex(const ElementPointer ElemP) const {return elementAtBoundary_[mapperEl_.map(*ElemP)]; }

  // Ny operasjon i IRISDuneGridInterface. (Bør i første omgang kun (innvendig) implementeres for ElementPointer og codimEntity2==1)
  template <int codimEntity2>
  bool boundaryIndex(const ElementPointer ElmP, const int& Ent2RelIdx) const { return boundaryIndex(ElmP->template subEntity<codimEntity2>(Ent2RelIdx));}


  //A global edge index's actual position at the boundary (NOTE: The outer boundary is here considered as a connsecutive Dune::array...)
  int findBoundaryPosition(const FacePointer khatP) const;
  //-----------------------------------------------------------------------------------

  private :

  //-----------------Useful local functions--------------------------------------------

  // We don't need copy constructor and assignment operator ...
  IRISDuneGridInterface& operator=(const IRISDuneGridInterface& IGI); // Not implemented ...
  IRISDuneGridInterface(const IRISDuneGridInterface& IGI); // Not implemented ...
  RnPoint convertToRnPoint(const RnVector& vec) const;

  ElementPointer getEntP(const ElementPointer & ElmP, const VertexPointer & VtxP, const int& EntRelIdx) const;
  FacePointer getEntP(const FacePointer & FaceP, const VertexPointer & VtxP, const int& EntRelIdx) const;

  //-----------------Data structure----------------------------------------------------

  const DuneGrid & Dgrid_;

  GridType gridType_;

  //Extra index stuff:
  std::vector<std::vector<ElementPointer> > ElementNeighboursOfVertex_;
  std::vector<std::vector<FacePointer> > FaceNeighboursOfVertex_;
  std::vector<std::vector<std::vector<ElementPointer> > > ElementNeighboursOfVertexThroughFaces_;

  /*
  NestedArray<ElementPointer> ElementNeighboursOfVertex_;
  NestedArray<EdgePointer> EdgeNeighboursOfVertex_;
  DoublyNestedArray<ElementPointer> ElementNeighboursOfVertexThroughEdges_;
  */

  //Extra boundary stuff:
  /*
  int numbSubBoundaries_;
  std::vector<int> numbOfEdgesAtSubBoundary_;//NOTE: The length of this Dune::array is obviously: numbSubBoundaries_;
  std::vector<EdgePointer> EdgesAtSubBoundary0_;//NOTE: The length of this Dune::array is obviously: numbOfEdgesAtSubBoundary_[0];
  std::vector<EdgePointer> EdgesAtSubBoundary1_;//NOTE: The length of this Dune::array is obviously: numbOfEdgesAtSubBoundary_[1];
  std::vector<EdgePointer> EdgesAtSubBoundary2_;//NOTE: The length of this Dune::array is obviously: numbOfEdgesAtSubBoundary_[2];
  std::vector<EdgePointer> EdgesAtSubBoundary3_;//NOTE: The length of this Dune::array is obviously: numbOfEdgesAtSubBoundary_[3];
  std::vector<EdgePointer> EdgesAtSubBoundary4_;//NOTE: Length=0 in 2D, else length of this Dune::array is obviously: numbOfEdgesAtSubBoundary_[4];
  std::vector<EdgePointer> EdgesAtSubBoundary5_;//NOTE: Length=0 in 2D, else length of this Dune::array is obviously: numbOfEdgesAtSubBoundary_[5];
  */
  std::vector<std::vector<FacePointer> > facesAtSubBoundary_;
  std::vector<bool> vertexAtBoundary_;
  std::vector<bool> faceAtBoundary_;
  std::vector<bool> elementAtBoundary_;

  //! Parameter for mapper class
  template<int dim>
  struct PElLayout
  {
    bool contains (Dune::GeometryType gt)
    {
     if (gt.dim()==dim) return true;
     return false;
    }
  };
  template<int dim>
  struct PFaLayout
  {
    bool contains (Dune::GeometryType gt)
    {
     if (gt.dim()==dim-1) return true;
     return false;
    }
  };
  template<int dim>
  struct PEdLayout
  {
    bool contains (Dune::GeometryType gt)
    {
     if (gt.dim()==1) return true;
     return false;
    }
  };

  template<int dim>
  struct PVLayout
  {
    bool contains (Dune::GeometryType gt)
    {
      if (gt.dim()==0) return true;
     return false;
    }
  };

  Dune::LevelMultipleCodimMultipleGeomTypeMapper<DuneGrid,PElLayout> mapperEl_;
  Dune::LevelMultipleCodimMultipleGeomTypeMapper<DuneGrid,PFaLayout> mapperFa_;
  Dune::LevelMultipleCodimMultipleGeomTypeMapper<DuneGrid,PEdLayout> mapperEd_;
  Dune::LevelMultipleCodimMultipleGeomTypeMapper<DuneGrid,PVLayout> mapperV_;

};




//-----------------------------------------------------
//              Constructors
//-----------------------------------------------------

template <class DuneGrid>
IRISDuneGridInterface<DuneGrid>::IRISDuneGridInterface(DuneGrid & grid)
:Dgrid_(grid)
,ElementNeighboursOfVertex_(Dgrid_.size(VertexPointer::codim))
,FaceNeighboursOfVertex_(Dgrid_.size(VertexPointer::codim))
,ElementNeighboursOfVertexThroughFaces_(Dgrid_.size(VertexPointer::codim))
,facesAtSubBoundary_(2*DuneGrid::dimension)
,vertexAtBoundary_(Dgrid_.size(VertexPointer::codim),false)
,faceAtBoundary_(Dgrid_.size(FacePointer::codim),false)
,mapperEl_(Dgrid_,0)
,mapperFa_(Dgrid_,0)
,mapperEd_(Dgrid_,0)
,mapperV_(Dgrid_,0)

{
  if (2 == DuneGrid::dimension)
    gridType_ = _Structured2D;
  else if (3 == DuneGrid::dimension)
    gridType_ = _Structured3D;
  else
    gridType_ = _somethingRotten;
}


//-----------------------------------------------------
//              Generators
//-----------------------------------------------------
template <class DuneGrid>
void IRISDuneGridInterface<DuneGrid>::setupForMPFAFPS_FVM()
{
  std::cout << "\n ------------------  Begin - setupForMPFAFPS_FVM() -------------------------" << std::endl;

  ////////////////////////////////////////////////
  // - ElementNeighboursOfVertex_
  // - FaceNeighboursOfVertex_
  // - ElementNeighboursOfVertexThroughEdges_
  ////////////////////////////////////////////////

  int NFace = Dgrid_.size(FacePointer::codim);
  int NVertex = Dgrid_.size(VertexPointer::codim);

  std::vector<std::vector<ElementPointer> > ElementNeighboursOfFace(NFace);
  std::vector<bool> faceDone(NFace, false);

  std::vector<RnVector> vtxCoord(NVertex);
  std::vector<bool> vtxCoordDone(NVertex, false);

  for (ElementIterator it = setEntityPointerToFirst<ElementIterator>(0); it != entityPointerIsAtEnd<ElementIterator>(0); ++it)
  {
    ElementPointer ep = it;

   Dune::GeometryType gt = it->type();
   typedef typename DuneGrid::ctype ct;
   const int dim = DuneGrid::dimension;

    int nFace = ep->template count<FacePointer::codim>();
    int nVertex = ep->template count<VertexPointer::codim>();

    //ElementNeighboursOfVertex:
    for (int v=0; v < nVertex; v++)
    {
      VertexPointer vp = ep->template subEntity<VertexPointer::codim>(v);
      ElementNeighboursOfVertex_[mapperV_.map(*vp)].push_back(ep);

      if (!vtxCoordDone[mapperV_.map(*vp)])
      {
        vtxCoord[mapperV_.map(*vp)] = getCentroid_V(vp);
      }
    }

   //FaceNeighboursOfVertex:
   for (int f=0; f<nFace; ++f)
   {
     FacePointer fp = ep->template subEntity<FacePointer::codim>(f);
     if (! faceDone[mapperFa_.map(*fp)])
     {
       Dune::GeometryType gt = fp->type();
       typedef typename DuneGrid::ctype ct;
       const int faceDim = DuneGrid::dimension - FacePointer::codim;

       int nVertexFace = Dune::GenericReferenceElements<ct,faceDim>::general(gt).size(faceDim);
       for (int v=0; v < nVertexFace; v++)
        {
          VertexPointer vp = ep->template subEntity<VertexPointer::codim>(
           Dune::GenericReferenceElements<ct,dim>::general(ep->type()).subEntity(f, FacePointer::codim, v, VertexPointer::codim));

          FaceNeighboursOfVertex_[mapperV_.map(*vp)].push_back(fp);
        }

       faceDone[mapperFa_.map(*fp)] = true;
     }
     ElementNeighboursOfFace[mapperFa_.map(*fp)].push_back(ep);
   }
  }

  // Sortering cc-wise...
  for (int v=0; v<NVertex; ++v)
  {
    RnVector interiorPoint = getCentroid_V(ElementNeighboursOfVertex_[v][0]);

    compareEntityCCWise<FacePointer> compareFace(vtxCoord[v],interiorPoint,this);
    std::stable_sort(FaceNeighboursOfVertex_[v].begin(),FaceNeighboursOfVertex_[v].end(),compareFace);

    compareEntityCCWise<ElementPointer> compareElement(vtxCoord[v],interiorPoint,this);
    std::stable_sort(ElementNeighboursOfVertex_[v].begin(),ElementNeighboursOfVertex_[v].end(),compareElement);
  }

  //ElementNeighboursOfVertexThroughFaces:
  for (int v=0; v<NVertex; ++v)
  {
    for(unsigned int f=0; f<FaceNeighboursOfVertex_[v].size(); ++f)
    {
      ElementNeighboursOfVertexThroughFaces_[v].push_back(ElementNeighboursOfFace[mapperFa_.map(*(FaceNeighboursOfVertex_[v][f]))]);

      if (ElementNeighboursOfVertexThroughFaces_[v][f].size() > 1)
      {
        // Sortering cc-wise...
        RnVector interiorPoint = getCentroid_V(FaceNeighboursOfVertex_[v][f]);
        compareEntityCCWise<ElementPointer> compareElement(vtxCoord[v],interiorPoint,this);
        std::stable_sort(ElementNeighboursOfVertexThroughFaces_[v][f].begin(),ElementNeighboursOfVertexThroughFaces_[v][f].end(),compareElement);
      }
    }
  }

  /*
  // Output
  std::cout << "Vertex->Element:" << std::endl;
  typename std::vector<std::vector<ElementPointer> >::const_iterator itV=ElementNeighboursOfVertex_.begin();
  while(itV !=  ElementNeighboursOfVertex_.end())
  {
    typename std::vector<ElementPointer>::const_iterator itElem=itV->begin();
    int idV = int(itV-ElementNeighboursOfVertex_.begin());
    std::cout << std::setw(2) << idV << ":";
    while(itElem != itV->end())
    {
      std::cout << "  " << mapperEl_.map(**itElem);
      ++itElem;
    }
    std::cout << std::endl;
    ++itV;
  }

  std::cout << "Vertex->Face:" << std::endl;
  typename std::vector<std::vector<FacePointer> >::const_iterator itVF=FaceNeighboursOfVertex_.begin();
  while(itVF !=  FaceNeighboursOfVertex_.end())
  {
    typename std::vector<FacePointer>::const_iterator itFace=itVF->begin();
    int idVF = int(itVF-FaceNeighboursOfVertex_.begin());
    std::cout << std::setw(2) << idVF << ":";
    while(itFace != itVF->end())
    {
      std::cout << "  " << mapperFa_.map(**itFace);
      ++itFace;
    }
    std::cout << std::endl;
    ++itVF;
  }

  std::cout << "Vertex->Face->Element:" << std::endl;
  typename std::vector<std::vector<std::vector<ElementPointer> > >::const_iterator itVFE=ElementNeighboursOfVertexThroughFaces_.begin();
  while(itVFE !=  ElementNeighboursOfVertexThroughFaces_.end())
  {
    typename std::vector<std::vector<ElementPointer> >::const_iterator itFace=itVFE->begin();
    int idVFE = int(itVFE-ElementNeighboursOfVertexThroughFaces_.begin());
    std::cout << std::setw(2) << idVFE << ":";
    while(itFace != itVFE->end())
    {
      int iFace = int(itFace-itVFE->begin());
      if (iFace>0) std::cout << "   ";
      std::cout << "  " << iFace << ":";
      typename std::vector<ElementPointer>::const_iterator itElem=itFace->begin();
      while (itElem != itFace->end())
      {
        std::cout << "  " << mapperEl_.map(**itElem);
        ++itElem;
      }
      std::cout << std::endl;
      ++itFace;
    }
    ++itVFE;
  }
  */

  ////////////////////////////////////////
  // - facesAtSubBoundary_
  // - vertexAtBoundary_
  // - faceAtBoundary_
  ////////////////////////////////////////

  // For now: Quick 'n dirty, assuming box-shaped domain and all elements of same orientation ...

  for (ElementIterator it = setEntityPointerToFirst<ElementIterator>(0); it != entityPointerIsAtEnd<ElementIterator>(0); ++it)
  {
    ElementPointer ep = it;
    if (ep->hasBoundaryIntersections())
    {
      //for (IntersectionIterator itFace=ep->ileafbegin(); itFace!=ep->ileafend(); ++itFace)
      for (IntersectionIterator itFace=ep->ilevelbegin(); itFace!=ep->ilevelend(); ++itFace)
      {
        if (itFace->boundary())
        {
          int iFace = itFace->indexInInside();
          FacePointer fp = ep->template subEntity<FacePointer::codim>(iFace);
          int externalFaceOrdering = iFace;
          if (externalFaceOrdering == 1)
          {
            externalFaceOrdering = 2;
          }
          else if (externalFaceOrdering == 2)
          {
            externalFaceOrdering = 1;
          }
          facesAtSubBoundary_[externalFaceOrdering].push_back(fp);
          faceAtBoundary_[mapperFa_.map(*fp)] = true;
          int numberOfVrtxOnFace = Dune::GenericReferenceElements<typename DuneGrid::ctype,DuneGrid::dimension>::general(ep->type()).size(iFace,1,VertexPointer::codim);
          for (int i=0; i<numberOfVrtxOnFace; ++i)
          {
            int iVrtx = Dune::GenericReferenceElements<typename DuneGrid::ctype,DuneGrid::dimension>::general(ep->type()).subEntity(iFace, 1, i, VertexPointer::codim);
            VertexPointer vP= ep->template subEntity<VertexPointer::codim>(iVrtx);
            vertexAtBoundary_[mapperV_.map(*vP)] = true;
          }
        }
      }
    }
  }

  /*
  // Output
  for (unsigned int i=0; i<facesAtSubBoundary_.size(); ++i)
  {
    std::cout << "\nBnd-" << i << ":";
    for (unsigned int j=0; j<facesAtSubBoundary_[i].size(); ++j)
    {
      std::cout << " " << mapperFa_.map(*facesAtSubBoundary_[i][j]);
    }
    std::cout << std::endl;
  }

  std::cout << "Faces:";
  for (unsigned int i=0; i<faceAtBoundary_.size(); ++i)
    std::cout << faceAtBoundary_[i];
  std::cout << std::endl;


  std::cout << "Vrtx:";
  for (unsigned int i=0; i<vertexAtBoundary_.size(); ++i)
    std::cout << vertexAtBoundary_[i];
  std::cout << std::endl;
  */

  std::cout << "\n ------------------  End - setupForMPFAFPS_FVM() -------------------------" << std::endl;

}

//-----------------------------------------------------
//              Iterator operations
//-----------------------------------------------------
template <class DuneGrid>
 template <class EntityIterator>
EntityIterator IRISDuneGridInterface<DuneGrid>::setEntityPointerToFirst(const int& gridLevel) const
{
  return Dgrid_.levelView(gridLevel).template begin<EntityIterator::codim>();
}

template <class DuneGrid>
 template <class EntityIterator>
EntityIterator IRISDuneGridInterface<DuneGrid>::entityPointerIsAtEnd(const int& gridLevel) const
{
  return Dgrid_.levelView(gridLevel).template end<EntityIterator::codim>();
}


//-----------------------------------------------------
//              Observers
//-----------------------------------------------------

template <class DuneGrid>
int IRISDuneGridInterface<DuneGrid>::domainDimension() const
{
  return DuneGrid::dimensionworld;
}

template <class DuneGrid>
int
IRISDuneGridInterface<DuneGrid>::cellDimension() const
{
  return DuneGrid::dimension;
}

template <class DuneGrid>
int
IRISDuneGridInterface<DuneGrid>::cellCount(const int& gridLevel, const int& codimEntity) const
{
  return Dgrid_.size(gridLevel,codimEntity);
}

template <class DuneGrid>
template <class EntityPointer>
typename IRISDuneGridInterface<DuneGrid>::RnVector
IRISDuneGridInterface<DuneGrid>::getCentroid_V(const EntityPointer entP) const
{
  Dune::FieldVector<typename DuneGrid::ctype,DuneGrid::dimension-EntityPointer::codim> cogLocal =
    Dune::GenericReferenceElements<typename DuneGrid::ctype,DuneGrid::dimension-EntityPointer::codim>::general(entP->type()).position(0,0);
  RnVector cogGlobal = (*entP).geometry().global(cogLocal);
  //std::cout << "cogLocal= (" << cogLocal << ") - " << "cogGlobal= (" << cogGlobal << ") - " << entP->type() << " - " << EntityPointer::codim << std::endl;
  return cogGlobal;
}

template <class DuneGrid>
template <int codimEntity2>
typename IRISDuneGridInterface<DuneGrid>::RnVector
IRISDuneGridInterface<DuneGrid>::getCentroid_V(const IRISDuneGridInterface<DuneGrid>::ElementPointer ElmP, const int& Ent2RelIdx) const
{
  return getCentroid_V(ElmP->template subEntity<codimEntity2>(Ent2RelIdx));
}

template <class DuneGrid>
double
IRISDuneGridInterface<DuneGrid>::getVolume(const IRISDuneGridInterface<DuneGrid>::ElementPointer entP) const
{
  return (*entP).geometry().volume();
}

template <class DuneGrid>
template <class EntityPointer>
int
IRISDuneGridInterface<DuneGrid>::getGlobalIndexFromEntityPointer(const EntityPointer EntP) const
{
  if (EntityPointer::codim == 0)
    return mapperEl_.map(*EntP);
  else if (EntityPointer::codim == 1)
    return mapperFa_.map(*EntP);
  else if (EntityPointer::codim == DuneGrid::dimension-1)
    return mapperEd_.map(*EntP);
  else if (EntityPointer::codim == int(DuneGrid::dimension))
    return mapperV_.map(*EntP);
}

template <class DuneGrid>
template <class EntityPointer>
int IRISDuneGridInterface<DuneGrid>::getNumberOfNeighbours(const EntityPointer EntP, const int& codimNeighbour) const
{
  assert(EntityPointer::codim == int(DuneGrid::dimension));

  if (codimNeighbour == 0)
    return ElementNeighboursOfVertex_[mapperV_.map(*EntP)].size();
  else if (codimNeighbour == 1)
     return FaceNeighboursOfVertex_[mapperV_.map(*EntP)].size();
  else
    assert(false);
}

template <class DuneGrid>
void
IRISDuneGridInterface<DuneGrid>::getNeighbours(const IRISDuneGridInterface<DuneGrid>::VertexPointer VtxP, std::vector<typename IRISDuneGridInterface<DuneGrid>::ElementPointer>& neighbours) const
{
  neighbours = ElementNeighboursOfVertex_[mapperV_.map(*VtxP)];
}

template <class DuneGrid>
void
IRISDuneGridInterface<DuneGrid>::getNeighbours(const IRISDuneGridInterface<DuneGrid>::VertexPointer VtxP, std::vector<typename IRISDuneGridInterface<DuneGrid>::FacePointer>& neighbours) const
{
  neighbours = FaceNeighboursOfVertex_[mapperV_.map(*VtxP)];
}


template <class DuneGrid>
void
IRISDuneGridInterface<DuneGrid>::getNeighbours(const IRISDuneGridInterface<DuneGrid>::ElementPointer ElmP, std::vector<typename IRISDuneGridInterface<DuneGrid>::ElementPointer>& neighbours) const
{
  assert(neighbours.empty());
  /*
  const int codimVtx = IRISDuneGridInterface<DuneGrid>::VertexPointer::codim;
  int nVtx = Dune::ReferenceElements<typename DuneGrid::ctype,DuneGrid::dimension>::general(ElmP->type()).size(codimVtx);
  for (int i=0; i<nVtx; ++i)
  {
    typename IRISDuneGridInterface<DuneGrid>::VertexPointer vtxP = getGlobalEntityPointerFromRelativeIndex<typename IRISDuneGridInterface<DuneGrid>::VertexPointer>(ElmP,i);
    for (unsigned int ii=0; ii<ElementNeighboursOfVertex_[mapperV_.map(*vtxP)].size(); ++ii)
    {
      if (ElementNeighboursOfVertex_[mapperV_.map(*vtxP)][ii] != ElmP)
      {
        if (neighbours.end() == find(neighbours.begin(),neighbours.end(),ElementNeighboursOfVertex_[mapperV_.map(*vtxP)][ii]) )
        {
          neighbours.push_back(ElementNeighboursOfVertex_[mapperV_.map(*vtxP)][ii]);
        }
      }
    }
  }
  */
    //for (IntersectionIterator itFace=ElmP->ileafbegin(); itFace!=ElmP->ileafend(); ++itFace)
  for (IntersectionIterator itFace=ElmP->ilevelbegin(); itFace!=ElmP->ilevelend(); ++itFace)
  {
    if (!itFace->boundary())
    {
      neighbours.push_back(itFace->outside());
    }
  }
}

template <class DuneGrid>
void
IRISDuneGridInterface<DuneGrid>::getNeighbours(const IRISDuneGridInterface<DuneGrid>::ElementPointer ElmP, std::vector<typename IRISDuneGridInterface<DuneGrid>::FacePointer>& neighbours) const
{
  assert(neighbours.empty());
  const int codimFace = IRISDuneGridInterface<DuneGrid>::FacePointer::codim;
  int nFace = Dune::GenericReferenceElements<typename DuneGrid::ctype,DuneGrid::dimension>::general(ElmP->type()).size(codimFace);
  for (int i=0; i<nFace; ++i)
  {
    neighbours.push_back(getGlobalEntityPointerFromRelativeIndex<typename IRISDuneGridInterface<DuneGrid>::FacePointer>(ElmP,i));
  }
}

template <class DuneGrid>
bool
IRISDuneGridInterface<DuneGrid>::isNeighbour(const IRISDuneGridInterface<DuneGrid>::VertexPointer VtxP, const IRISDuneGridInterface<DuneGrid>::FacePointer FaceP) const
{
  int iVtx = mapperV_.map(*VtxP);
  for (unsigned int i=0; i<FaceNeighboursOfVertex_[iVtx].size(); ++i)
  {
    if (FaceNeighboursOfVertex_[iVtx][i] == FaceP)
      return true;
  }
  return false;
}

template <class DuneGrid>
bool
IRISDuneGridInterface<DuneGrid>::isNeighbour(const IRISDuneGridInterface<DuneGrid>::ElementPointer ElmP, const IRISDuneGridInterface<DuneGrid>::FacePointer FaceP) const
{
  const int codimFace = IRISDuneGridInterface<DuneGrid>::FacePointer::codim;
  int nFace = Dune::GenericReferenceElements<typename DuneGrid::ctype,DuneGrid::dimension>::general(ElmP->type()).size(codimFace);
  for (int i=0; i<nFace; ++i)
  {
    if (ElmP->template entity<codimFace>(i) == FaceP)
      return true;
  }
  return false;
}

template <class DuneGrid>
typename IRISDuneGridInterface<DuneGrid>::RnVector
IRISDuneGridInterface<DuneGrid>::getNormalVector_V(const ElementPointer ihatP, const int& EdgeRelIdx) const
{
  //NOTE: An important assumtion must be that EdgeRelIdx is completely compatible with the local ordering that the IntersectionIterator employs.
  //      Obviously the IntersectionIterator class must be used in the implementation of this function.  RnVector ...

  //IntersectionIterator itFace = ihatP->ileafbegin();
  IntersectionIterator itFace = ihatP->ilevelbegin();
  for (int i=0; i<EdgeRelIdx;  ++i) ++itFace; // No operator+ for IntersectionIterator ...
  Dune::FieldVector<typename DuneGrid::ctype, DuneGrid::dimension-1> coordCogFace(0);
  RnVector nGlobal = itFace->unitOuterNormal(coordCogFace); // Velger å tro at dette er globale coord ...
  //std::cout << "Outer unit normal:     coordCogFace= (" << coordCogFace << ") - " << "nGlobal= (" << nGlobal << ")" << std::endl;
  return nGlobal;
}

//NOTE: This code will NOT work for a triangular grid!!!
template <class DuneGrid>
double
IRISDuneGridInterface<DuneGrid>::getMinGridSize() const
{
  double diamMin = 1.0e6;
  for (ElementIterator it = setEntityPointerToFirst<ElementIterator>(0); it != entityPointerIsAtEnd<ElementIterator>(0); ++it)
  {
    ElementPointer ep = it;
    double dM = (dim==2) ? std::pow(getVolume(ep),0.5) : std::pow(getVolume(ep),0.33333333);
    if (diamMin > dM)
      diamMin = dM;
  }
  return diamMin;
}

//-----------------------------------------------------------------------------
//A group of more topological type of operations:
//-----------------------------------------------------------------------------

//A relative index function which delivers the corresponding (global) entity pointer:
template <class DuneGrid>
template <class EntityPointer>
EntityPointer
IRISDuneGridInterface<DuneGrid>::getGlobalEntityPointerFromRelativeIndex(const IRISDuneGridInterface<DuneGrid>::ElementPointer & ElmP, const int& EntRelIdx) const
{
  return ElmP->template subEntity<EntityPointer::codim>(EntRelIdx);
}

template <class DuneGrid>
template <class EntityPointer>
EntityPointer
IRISDuneGridInterface<DuneGrid>::getGlobalEntityPointerFromRelativeIndex(const IRISDuneGridInterface<DuneGrid>::VertexPointer & VtxP, const int& EntRelIdx) const
{
  EntityPointer entP = getEntP(entP,VtxP,EntRelIdx);
  return entP;
}

template <class DuneGrid>
typename IRISDuneGridInterface<DuneGrid>::ElementPointer
IRISDuneGridInterface<DuneGrid>::getEntP(const IRISDuneGridInterface<DuneGrid>::ElementPointer & ElmP, const IRISDuneGridInterface<DuneGrid>::VertexPointer & VtxP, const int& EntRelIdx) const
{
  assert(EntRelIdx < int(ElementNeighboursOfVertex_[mapperV_.map(*VtxP)].size()));
  return ElementNeighboursOfVertex_[mapperV_.map(*VtxP)][EntRelIdx];
}

template <class DuneGrid>
typename IRISDuneGridInterface<DuneGrid>::FacePointer
IRISDuneGridInterface<DuneGrid>::getEntP(const IRISDuneGridInterface<DuneGrid>::FacePointer & FaceP, const IRISDuneGridInterface<DuneGrid>::VertexPointer & VtxP, const int& EntRelIdx) const
{
  assert(EntRelIdx < int(FaceNeighboursOfVertex_[mapperV_.map(*VtxP)].size()));
  return FaceNeighboursOfVertex_[mapperV_.map(*VtxP)][EntRelIdx];
}

//An overloaded version of the above function:
template <class DuneGrid>
template <class EntityPointer>
EntityPointer
IRISDuneGridInterface<DuneGrid>::getGlobalEntityPointerFromRelativeIndex(const IRISDuneGridInterface<DuneGrid>::ElementPointer & ElmP, const int& Ent2RelIdx, const int& codimEntity2, const int& Ent3RelIdx) const
{
  return ElmP->template entity<EntityPointer::codim>(
    Dune::GenericReferenceElements<typename DuneGrid::ctype,DuneGrid::dimension>::general(ElmP->type()).subEntity(Ent2RelIdx, codimEntity2, Ent3RelIdx, EntityPointer::codim));
}

template <class DuneGrid>
typename IRISDuneGridInterface<DuneGrid>::ElementPointer
IRISDuneGridInterface<DuneGrid>::getGlobalEntityPointerFromRelativeIndex(const IRISDuneGridInterface<DuneGrid>::VertexPointer & VtxP, const int& FaceRelIdx, const int& codimFace, const int& ElemRelIdx) const
{
  assert(FaceRelIdx < int(ElementNeighboursOfVertexThroughFaces_[mapperV_.map(*VtxP)].size()) && codimFace==1
         && ElemRelIdx < int(ElementNeighboursOfVertexThroughFaces_[mapperV_.map(*VtxP)][FaceRelIdx].size()));
  return ElementNeighboursOfVertexThroughFaces_[mapperV_.map(*VtxP)][FaceRelIdx][ElemRelIdx];
}

//Another type of relative index function which delivers a relative index:
template <class DuneGrid>
template <class EntityPointer>
int
IRISDuneGridInterface<DuneGrid>::getDecreasedRelativeIndexFromRelativeIndex(const EntityPointer EntP, const int& Ent2RelIdx, const int& codimEntity2, const int& Ent3RelIdx, const int& codimEntity3) const
{
  return Dune::GenericReferenceElements<typename DuneGrid::ctype,DuneGrid::dimension-EntityPointer::codim>::general(EntP->type()).subEntity(Ent2RelIdx, codimEntity2, Ent3RelIdx, codimEntity3);
}

template <class DuneGrid>
int
IRISDuneGridInterface<DuneGrid>::getDecreasedRelativeIndexFromRelativeIndex(const IRISDuneGridInterface<DuneGrid>::VertexPointer & VtxP, const int&FaceRelVtxIdx, const int& ElemRelFaceIdx) const
{
  int iVtx = mapperV_.map(*VtxP);
  for (unsigned int iElem=0; iElem<ElementNeighboursOfVertex_[iVtx].size(); ++iElem)
  {
    if (ElementNeighboursOfVertex_[iVtx][iElem] == ElementNeighboursOfVertexThroughFaces_[iVtx][FaceRelVtxIdx][ElemRelFaceIdx])
      return iElem;
  }
  return -1;
}

template <class DuneGrid>
void
IRISDuneGridInterface<DuneGrid>::getIntersectionRelativeEdgeIndexSetRightHandCSSorted(const VertexPointer qhatP, const int& ElementGlobIdx, std::vector<int> & out) const
{
  //Just to be sure:
  out.clear();
  // Identify element:
  int iVtx = mapperV_.map(*qhatP);
  typename std::vector<ElementPointer>::const_iterator itElem=ElementNeighboursOfVertex_[iVtx].begin();
  while (itElem != ElementNeighboursOfVertex_[iVtx].end())
  {
    if ( mapperEl_.map(**itElem) == ElementGlobIdx)
      break;
    ++itElem;
  }
  assert(itElem != ElementNeighboursOfVertex_[iVtx].end());
  ElementPointer ep=*itElem;

  for (unsigned int iFace=0; iFace<ElementNeighboursOfVertexThroughFaces_[iVtx].size(); ++iFace)
  {
    for (unsigned int iElem=0; iElem<ElementNeighboursOfVertexThroughFaces_[iVtx][iFace].size(); ++iElem)
    {
      if (ElementNeighboursOfVertexThroughFaces_[iVtx][iFace][iElem] == ep)
      {
        out.push_back(iFace);
      }
    }
  }

  if (!vertexAtBoundary_[iVtx])
  {
    int iCut=0;
    for (unsigned int i=1; i<out.size(); ++i)
    {
      if (out[i]-out[i-1] > 1)
      {
        iCut = i;
        break;
      }
    }
    if (iCut > 0)
    {
      if (out.size() == 2)
      {
        std::swap(out[0],out[1]);
      }
      else
      {
        std::vector<int> outTmp(out.size());
        for (unsigned int i=0; i<out.size(); ++i)
          outTmp[i] = out[(i+iCut)%(out.size()+1)];
        out = outTmp;
      }
    }
  }
  //--------------------------HAF: 12.02.2010----------------------------------
  //Change to a clockwise rotation:
  if (out.size() == 2)
  {
    int help = out[0];
    out[0] = out[1];
    out[1] = help;
  }
  //---------------------------------------------------------------------------
}


// Boundaries ...

template <class DuneGrid>
int
IRISDuneGridInterface<DuneGrid>::findBoundaryPosition(const IRISDuneGridInterface<DuneGrid>::FacePointer khatP) const
{
  int pos=-1;
  typename std::vector<std::vector<FacePointer> >::const_iterator itSub = facesAtSubBoundary_.begin();
  while(itSub != facesAtSubBoundary_.end())
  {
    typename std::vector<FacePointer>::const_iterator itFace = itSub->begin();
    while (itFace != itSub->end())
    {
      ++pos;
      if (khatP == *itFace)
        return pos;
      ++itFace;
    }
    ++itSub;
  }
  return -1;
}

template <class DuneGrid>
RnPoint
IRISDuneGridInterface<DuneGrid>::convertToRnPoint(const RnVector& vec) const
{
  RnShape sh(DuneGrid::dimensionworld);
  RnPoint rp;
  rp.u_setPoint_1u(sh, 0.0); //initialiserer koords. til 0.0.
  for (int i=0; i<DuneGrid::dimensionworld; ++i)
  {
    rp.u_changeElement_1u(i, vec[i]);
  }
  return rp;
}


/*
===============================================================================
    ------------------------ CLASS DESCRIPTION --------------------------
===============================================================================

  DESCRIPTION:
  ------------------

    (NOTE: Should also decribe Boundary convensions etc......)



  USAGE/EXAMPLES:
  ---------------
    <Examples how to use the class>


  MISCELLANEOUS:
  --------------






===============================================================================
    --------------------- END OF CLASS DESCRIPTION ----------------------
===============================================================================
*/

#endif


