
#ifndef XFEM_CIRCLE_CUT_H
#define XFEM_CIRCLE_CUT_H

#include "XFEM_geometric_cut.h"

class XFEM_circle_cut : public XFEM_geometric_cut
{
public:

  XFEM_circle_cut(std::vector<Real> square_nodes);
  ~XFEM_circle_cut();

  virtual bool cut_elem_by_geometry(const Elem* elem, std::vector<cutEdge> & cutEdges, Real time);
  virtual bool cut_elem_by_geometry(const Elem* elem, std::vector<cutFace> & cutFaces, Real time);

  virtual bool cut_frag_by_geometry(std::vector<std::vector<Point> > & frag_edges, std::vector<cutEdge> & cutEdges, Real time);
  virtual bool cut_frag_by_geometry(std::vector<std::vector<Point> > & frag_faces, std::vector<cutFace> & cutFaces, Real time);

private:

  std::vector<Point> _vertices;
  Point _center;
  Point _normal;
  Real _radius;
  Real _angle;

private:

  bool intersect_with_edge(Point p1, Point p2, Point &pint);
  bool isInsideCutPlane(Point p);
};

#endif
