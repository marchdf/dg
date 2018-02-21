/*!
  \file simpleMesh.h  
  \brief Class deals with the mesh.
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \license This project is released under the GNU Public License. See LICENSE.
  \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
  \section Description
  Classes that set up the mesh (read/write), build the interfaces,
  builds the neighbors, etc
*/
#ifndef _SIMPLE_MESH_H_
#define _SIMPLE_MESH_H_
#include <vector>
#include "fullMatrix.h"
#include "scalar_def.h"
#include <string>
#include <map>
#include "macros.h"
#include <math.h>
#include "philmesh.h"

/*!
  \class simpleElement simpleMesh.h
  \brief Class that creates a simple element
*/  
class simpleElement {
  std::vector<int> _nodes;
  int _ptag, _id, _partition; // partition is numbered from 0 to match with MPI processes
  public:
  /*!
    \brief Constructor of a simple element
    \param[in] id  its id
    \param[in] ptag its physical tag
    \param[in] partition its partition
    \param[in] nodes the nodes of the element
  */
  simpleElement (int id, int ptag, int partition, std::vector<int> nodes) : _id (id), _ptag (ptag), _partition(partition), _nodes (nodes) {}
  inline int getId () const {/*!Return the element id*/return _id; }
  inline int getPhysicalTag () const {/*!Return the physical tag of an element*/return _ptag; }
  inline int getPartition () const {/*!Return the partition the element belongs to*/ return _partition; }
  inline int getNbNodes () const {/*!Return the number of nodes of the element*/return _nodes.size() ;}
  inline int getNode (int i) const {/*!Return nodes of an element*/return _nodes[i];}
};

class simpleMesh;

/*!
  \class simpleInterface simpleMesh.h
  \brief Class that creates a simple interface
*/  
class simpleInterface {
  const simpleElement *_elements[2];
  int _closureId[2];
  int _physicalTag;
 public:
  inline int getPhysicalTag() const {/*!return physical tag of interface*/return _physicalTag;}
  inline const simpleElement *getElement(int i) const {/*!return the element of a side of the interface*/return _elements[i];}
  inline int getClosureId(int i) const {/*!return the closure id of a side of the interface*/return _closureId[i];}
  static void BuildInterfaces(simpleMesh &mesh, std::vector<simpleInterface> &interfaces, int typeInterface, int typeElement, int nsides, int PeriCo2D);
  static void BuildInterfaces_PERIFIX(simpleMesh &mesh, std::vector<simpleInterface> &interfaces, int typeInterface, int typeElement, int nsides);
  static void FixClosures(simpleMesh &mesh, std::vector<simpleInterface> &interfaces, int typeInterface, int typeElement, int nsides, int PeriCo2D, scalar DomainLength);
  simpleInterface(int physicalTag = 0); //I think this is how we declare interfaces
};

class simpleMesh {
  int _myid; // CPU identity for MPI
  int _numprocs; // total number of processors in use
  std::vector<std::vector<simpleElement> > _elements;      // holds elements in my partition
  std::vector<std::vector<simpleElement> > _otherElements; // holds elements in other partitions
  std::vector<simpleInterface> _interfaces;
  fullMatrix<double> _nodes;
  fullMatrix<scalar> _normals;
  std::map<int,int> _elementMap;
  std::map<int,int> _ghostElementMap;

  //PEJ 12/01/2017: Messing with the communication routines today,
  //going to make _ghostElementMap obsolete.
  std::vector<int> _ghostElement_Key;
  std::vector<int> _ghostElement_Val;
  std::vector<int> _ghostElement_HNe;

  bool _cartesian;
  int* _ghostElementSend;
  int* _ghostElementRecv;
  int _N_B; // number of boundaries
  int _N_ghosts;
  int* _boundary;
  int* _boundaryIdx;
  int* _neighbors; // N_N x N_E : element | neighbor1 | neighbor2 | ...
  fullMatrix<scalar> _shifts;
  scalar _Dx;
 public:
  /*!
    \brief Constructor for simpleMesh
    \param[in] myid my processor id (default=0)
    \param[in] numprocs number of processors (default=1)
  */
  simpleMesh(int myid=0, int numprocs=1) : _myid (myid), _numprocs(numprocs){}
  inline const std::vector<simpleInterface> & getInterfaces () const {/*! Return interfaces*/return  _interfaces;}
  inline const std::vector<simpleElement> & getElements (int type) const {/*! Return elements*/return  _elements[type];}
  inline const std::vector<simpleElement> & getOtherElements (int type) const {/*! Return elements not in my partition*/return  _otherElements[type];}
  inline const fullMatrix<double> & getNodes () const {/*! Return nodes*/return _nodes;}
  inline const fullMatrix<scalar> & getNormals () const {/*! Return normals*/return _normals;}
  inline const std::map<int,int> & getElementMap () const {/*! Return elements mappings*/return _elementMap;}
  inline const std::map<int,int> & getGhostElementMap () const {/*! Return ghost elements mappings*/return _ghostElementMap;}
  inline const std::vector<int> & getGhostElement_Key () const {/*! Return ghost element key*/return _ghostElement_Key;}
  inline const std::vector<int> & getGhostElement_Val () const {/*! Return ghost element values*/return _ghostElement_Val;}
  inline const std::vector<int> & getGhostElement_HNe () const {/*! Return ghost element values*/return _ghostElement_HNe;}
  void load (const char *fileName);
  void writeSolution (const scalar* solution, const int N_s, const int N_E, int type, std::vector<std::string> fnames, std::vector<std::string> names, int step, double time) const;
  void readSolution (const int N_s, const int N_E, int type, std::vector<std::string> fnames, std::vector<std::string> names, const int step, double &time, scalar* solution);
  void readSolutionNoTime (const int N_s, const int N_E, int type, std::vector<std::string> fnames, std::vector<std::string> names, const int step, scalar* solution);
  inline void buildInterfaces(int typeInterface, int typeElement, int nsides, int PeriCo2D) {
    /*! Calls BuildInterfaces*/
    simpleInterface::BuildInterfaces(*this, _interfaces, typeInterface, typeElement, nsides, PeriCo2D);
  }
  inline void buildInterfaces_PERIFIX(int typeInterface, int typeElement, int nsides) {
    /*! Calls BuildInterfaces*/
    simpleInterface::BuildInterfaces_PERIFIX(*this, _interfaces, typeInterface, typeElement, nsides);
  }
  
  inline void fixclosures(int typeInterface, int typeElement, int nsides, int PeriCo2D, scalar DomainLength) {
    simpleInterface::FixClosures(*this, _interfaces, typeInterface, typeElement, nsides, PeriCo2D, DomainLength);
  }

  void buildElementMap(int elem_type);
  void buildCommunicators(int elem_type);
  void buildNormals(int typeInterface, int typeElement, int PeriCo2D, std::vector<std::vector<int> > elemNodes, int order, int N_N);

  int getNbGhosts() const {/*! Return number of ghost elements*/return _N_ghosts;}
  int* getGhostElementSend()const {/*! Return ghost elements to send*/return _ghostElementSend;}
  int* getGhostElementRecv()const {/*! Return ghost elements to send*/return _ghostElementRecv;}
  
  void buildBoundary();
  int  getBoundarySize()const {/*! Return boundary size*/return _N_B;}
  int* getBoundaryMap()const {/*! Return boundary map*/return _boundary;}
  int* getBoundaryIdx()const {/*! Return boundary index mesh*/return _boundaryIdx;}

  int SeekGhost(/*const simpleElement* elHome, int om_on, */int GmshId_on, int GmshId_off);
  void buildNeighbors(int N_N, int N_E);
  int* getNeighbors()const {/*! Return element neighbors mesh*/return _neighbors;}
  void buildBoundaryElementShift1D(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes);
  void buildBoundaryElementShift2D(int order, const fullMatrix<scalar> &XYZNodesF, std::map<int,int> &ElementMap);
  inline const fullMatrix<scalar> & getShifts() const {/*! Return shifts*/return _shifts;}

  fullMatrix<scalar> getElementCentroids(const int N_E, const int ncorners, const fullMatrix<scalar> XYZNodes);
  bool iscartesian(std::string typeElement, const int elem_type);
  void setDx(const int N_N, const int N_E, const fullMatrix<scalar> &XYZCen, const fullMatrix<scalar> &XYZNodes);
  scalar getDx()const {/*! Return smallest Dx in mesh*/ return _Dx;} 
  inline scalar distance_to_edge(scalar x0,scalar y0,scalar x1,scalar y1,scalar x2,scalar y2);
    
  // Destructor
  ~simpleMesh(){
    delete[] _ghostElementSend;
    delete[] _ghostElementRecv;
    delete[] _boundary;
    delete[] _boundaryIdx;
    delete[] _neighbors; 
  };

};
#endif
