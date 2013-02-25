#ifndef _SIMPLE_MESH_H_
#define _SIMPLE_MESH_H_
#include <vector>
#include "fullMatrix.h"
#include <scalar_def.h>
#include <string>
#include <map>

class simpleElement {
  std::vector<int> _nodes;
  int _tag, _id;
  int _physicalTag;
  public:
  simpleElement (int id, int tag, std::vector<int> nodes) : _id (id), _tag (tag), _nodes (nodes) {}
  inline int getId () const { return _id; }
  inline int getPhysicalTag () const { return _tag; }
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
  simpleInterface(int physicalTag = 0);
};

class simpleMesh {
  std::vector<std::vector<simpleElement> > _elements;
  std::vector<simpleInterface> _interfaces;
  fullMatrix<double> _nodes;
  fullMatrix<scalar> _normals;
  int _N_B; // number of boundaries
  int* _boundary;
  int* _neighbors; // N_N x N_E : element | neighbor1 | neighbor2 | ...
  fullMatrix<scalar> _shifts;
  fullMatrix<scalar> _XYZCen;
  public:
  inline const std::vector<simpleInterface> & getInterfaces () const {return  _interfaces;}
  inline const std::vector<simpleElement> & getElements (int type) const {return  _elements[type];}
  inline const fullMatrix<double> & getNodes () const {return _nodes;}
  inline const fullMatrix<scalar> & getNormals () const {return _normals;}
  void load (const char *fileName);
  void writeSolution (const fullMatrix<scalar> &solution, int type, const char *filename, const char *name, int step, double time, int append) const;
  inline void buildInterfaces(int typeInterface, int typeElement, int nsides) {
    simpleInterface::BuildInterfaces(*this, _interfaces, typeInterface, typeElement, nsides);
  }
  void buildNormals(int typeInterface, int typeElement, const int D);
  void buildPeriodicSquare(int order, const fullMatrix<scalar> &XYZNodesF, const int D);
  void buildPeriodicLine();
  void buildFarfield();
  int getBoundaryNB()  const {return _N_B;}
  int* getBoundaryMap()const {return _boundary;}
  void buildNeighbors(int N_N, int N_E, std::map<int,int> &ElementMap);
  void sortNeighbors(const int N_E, const int N_N, const fullMatrix<scalar> XYZCen);
  int* getNeighbors()const {return _neighbors;}
  void buildBoundaryElementShift(int order, const fullMatrix<scalar> &XYZNodesF, const int D, std::map<int,int> &ElementMap);
  inline const fullMatrix<scalar> & getShifts() const {return _shifts;}
  fullMatrix<scalar> getElementCentroids(const int N_E, const int D, const int ncorners, const fullMatrix<scalar> XYZNodes);
  bool iscartesian(std::string typeElement, const int elem_type);
};
#endif
