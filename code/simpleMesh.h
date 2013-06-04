#ifndef _SIMPLE_MESH_H_
#define _SIMPLE_MESH_H_
#include <vector>
#include "fullMatrix.h"
#include <scalar_def.h>
#include <string>
#include <map>
#include <macros.h>
#include <math.h>

class simpleElement {
  std::vector<int> _nodes;
  int _ptag, _id, _partition; // partition is numbered from 0 to match with MPI processes
  public:
  simpleElement (int id, int ptag, int partition, std::vector<int> nodes) : _id (id), _ptag (ptag), _partition(partition), _nodes (nodes) {}
  inline int getId () const { return _id; }
  inline int getPhysicalTag () const { return _ptag; }
  inline int getPartition () const { return _partition; }
  inline int getNbNodes () const { return _nodes.size() ;}
  inline int getNode (int i) const {return _nodes[i];}
};

class simpleMesh;
class simpleInterface {
  const simpleElement *_elements[2];
  int _closureId[2];
  int _physicalTag;
 public:
  inline int getPhysicalTag() const {return _physicalTag;}
  inline const simpleElement *getElement(int i) const {return _elements[i];}
  inline int getClosureId(int i) const {return _closureId[i];}
  static void BuildInterfaces(simpleMesh &mesh, std::vector<simpleInterface> &interfaces, int tagInterface, int tagElement, int nsides);
  //bool pairPeriodic(const fullMatrix<double> & meshNodes, std::vector<int> nodes1, std::vector<int> nodes2);
  simpleInterface(int physicalTag = 0);
};

class simpleMesh {
  int _myid; // CPU identity for MPI
  std::vector<std::vector<simpleElement> > _elements;      // holds elements in my partition
  std::vector<std::vector<simpleElement> > _otherElements; // holds elements in other partitions
  std::vector<simpleInterface> _interfaces;
  fullMatrix<double> _nodes;
  fullMatrix<scalar> _normals;
  std::map<int,int> _elementMap;
  std::map<int,int> _ghostElementMap;
  bool _cartesian;
  int* _ghostInterfaces;
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
  simpleMesh(int myid=0) : _myid (myid) {}
  inline const std::vector<simpleInterface> & getInterfaces () const {return  _interfaces;}
  inline const std::vector<simpleElement> & getElements (int type) const {return  _elements[type];}
  inline const std::vector<simpleElement> & getOtherElements (int type) const {return  _otherElements[type];}
  inline const fullMatrix<double> & getNodes () const {return _nodes;}
  inline const fullMatrix<scalar> & getNormals () const {return _normals;}
  inline const std::map<int,int> & getElementMap () const {return _elementMap;}
  inline const std::map<int,int> & getGhostElementMap () const {return _ghostElementMap;}
  void load (const char *fileName);
  void writeSolution (const fullMatrix<scalar> &solution, int type, std::string filename, const char *name, int step, double time, int append) const;
  inline void buildInterfaces(int typeInterface, int typeElement, int nsides) {
    simpleInterface::BuildInterfaces(*this, _interfaces, typeInterface, typeElement, nsides);
  }
  void buildElementMap(int elem_type);
  void buildCommunicators(int elem_type);
  void buildNormals(int typeInterface, int typeElement, const int D);

  int getNbGhostInterfaces() const {return _N_ghosts;}
  int* getGhostInterfaces()const {return _ghostInterfaces;}
  int* getGhostElementSend()const {return _ghostElementSend;}
  int* getGhostElementRecv()const {return _ghostElementRecv;}
  
  void buildLineBoundary(int boundaryType);
  void buildSquareBoundary(int M_s, const fullMatrix<scalar> &XYZNodesF, const int D);
  int  getBoundarySize()  const {return _N_B;}
  int* getBoundaryMap()const {return _boundary;}
  int* getBoundaryIdx()const {return _boundaryIdx;}

  void buildNeighbors(int N_N, int N_E);
  void sortNeighbors(const int N_E, const int N_N, const fullMatrix<scalar> XYZCen);
  int* getNeighbors()const {return _neighbors;}
  void buildBoundaryElementShift1D(const int N_s, const int D, const int N_E, const fullMatrix<scalar> &XYZNodes);
  void buildBoundaryElementShift2D(int order, const fullMatrix<scalar> &XYZNodesF, const int D, std::map<int,int> &ElementMap);
  inline const fullMatrix<scalar> & getShifts() const {return _shifts;}

  fullMatrix<scalar> getElementCentroids(const int N_E, const int D, const int ncorners, const fullMatrix<scalar> XYZNodes);
  bool iscartesian(std::string typeElement, const int elem_type);
  void setDx(const int N_N, const int N_E, const int D, const fullMatrix<scalar> &XYZCen, const fullMatrix<scalar> &XYZNodes);
  scalar getDx()const {return _Dx;}
  inline scalar distance_to_edge(scalar x0,scalar y0,scalar x1,scalar y1,scalar x2,scalar y2);
    
  // Destructor
  ~simpleMesh();

};
#endif
