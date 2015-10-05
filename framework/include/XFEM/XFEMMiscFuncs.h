#ifndef XFEMINTEGRATIONFUNCS_H
#define XFEMINTEGRATIONFUNCS_H

#include "libmesh/libmesh_common.h"
#include "libmesh/libmesh.h" // libMesh::invalid_uint
#include <iostream>
#include <vector>
#include <cmath>
#include "MooseError.h"
#include "Conversion.h"

using namespace libMesh;

void static dunavant_rule2(const Real* wts, const Real* a, const Real* b, const unsigned int* permutation_ids,
                           unsigned int n_wts, std::vector<Point> &points, std::vector<Real> &weights)
{
  // see libmesh/src/quadrature/quadrature_gauss.C
  // Figure out how many total points by summing up the entries
  // in the permutation_ids array, and resize the _points and _weights
  // vectors appropriately.
  unsigned int total_pts = 0;
  for (unsigned int p=0; p<n_wts; ++p)
    total_pts += permutation_ids[p];

  // Resize point and weight vectors appropriately.
  points.resize(total_pts);
  weights.resize(total_pts);

  // Always insert into the points & weights vector relative to the offset
  unsigned int offset=0;
  for (unsigned int p=0; p<n_wts; ++p)
  {
    switch (permutation_ids[p])
    {
      case 1:
      {
        // The point has only a single permutation (the centroid!)
        // So we don't even need to look in the a or b arrays.
        points [offset  + 0] = Point(1.0L/3.0L, 1.0L/3.0L);
        weights[offset + 0] = wts[p];

        offset += 1;
        break;
      }

      case 3:
      {
        // For this type of rule, don't need to look in the b array.
        points[offset + 0] = Point(a[p],         a[p]);         // (a,a)
        points[offset + 1] = Point(a[p],         1.L-2.L*a[p]); // (a,1-2a)
        points[offset + 2] = Point(1.L-2.L*a[p], a[p]);         // (1-2a,a)

        for (unsigned int j=0; j<3; ++j)
          weights[offset + j] = wts[p];

        offset += 3;
        break;
      }

      case 6:
      {
        // This type of point uses all 3 arrays...
        points[offset + 0] = Point(a[p], b[p]);
        points[offset + 1] = Point(b[p], a[p]);
        points[offset + 2] = Point(a[p], 1.L-a[p]-b[p]);
        points[offset + 3] = Point(1.L-a[p]-b[p], a[p]);
        points[offset + 4] = Point(b[p], 1.L-a[p]-b[p]);
        points[offset + 5] = Point(1.L-a[p]-b[p], b[p]);

        for (unsigned int j=0; j<6; ++j)
          weights[offset + j] = wts[p];

        offset += 6;
        break;
      }

      default:
        mooseError("Unknown permutation id: " << permutation_ids[p] << "!");
    }
  }
}

void static stdQuadr2D(unsigned int nen, unsigned int iord, std::vector<std::vector<Real> > &sg2) // ZZY code
{
  // Purpose: get Guass integration points for 2D quad and trig elems
  // N.B. only works for n_qp <= 6

  Real lr4[4] = {-1.0,1.0,-1.0,1.0}; // libmesh order
  Real lz4[4] = {-1.0,-1.0,1.0,1.0};
  Real lr9[9] = {-1.0,0.0,1.0,-1.0,0.0,1.0,-1.0,0.0,1.0}; // libmesh order
  Real lz9[9] = {-1.0,-1.0,-1.0,0.0,0.0,0.0,1.0,1.0,1.0};
  Real lw9[9] = {25.0,40.0,25.0,40.0,64.0,40.0,25.0,40.0,25.0};

  if (nen == 4) // 2d quad element
  {
    if (iord == 1) // 1-point Gauss
    {
      sg2.resize(1);
      sg2[0].resize(3);
      sg2[0][0] = 0.0;
      sg2[0][1] = 0.0;
      sg2[0][2] = 4.0;
    }
    else if (iord == 2) // 2x2-point Gauss
    {
      sg2.resize(4);
      for (unsigned int i = 0; i < 4; ++i) sg2[i].resize(3);
      for (unsigned int i = 0; i < 4; ++i)
      {
        sg2[i][0] = (1/sqrt(3))*lr4[i];
        sg2[i][1] = (1/sqrt(3))*lz4[i];
        sg2[i][2] = 1.0;
      }
    }
    else if (iord == 3) // 3x3-point Gauss
    {
      sg2.resize(9);
      for (unsigned int i = 0; i < 9; ++i) sg2[i].resize(3);
      for (unsigned int i = 0; i < 9; ++i)
      {
        sg2[i][0] = sqrt(0.6)*lr9[i];
        sg2[i][1] = sqrt(0.6)*lz9[i];
        sg2[i][2] = (1.0/81.0)*lw9[i];
      }
    }
    else
      mooseError("Invalid quadrature order = " + Moose::stringify(iord) + " for quad elements");
  }
  else if (nen == 3) // triangle
  {
    if (iord == 1) // one-point Gauss
    {
      sg2.resize(1);
      sg2[0].resize(4);
      sg2[0][0] = 1.0/3.0;
      sg2[0][1] = 1.0/3.0;
      sg2[0][2] = 1.0/3.0;
      sg2[0][3] = 0.5;
    }
    else if (iord == 2) // three-point Gauss
    {
      sg2.resize(3);
      for (unsigned int i = 0; i < 3; ++i) sg2[i].resize(4);
      sg2[0][0] = 2.0/3.0; sg2[0][1] = 1.0/6.0; sg2[0][2] = 1.0/6.0; sg2[0][3] = 1.0/6.0;
      sg2[1][0] = 1.0/6.0; sg2[1][1] = 2.0/3.0; sg2[1][2] = 1.0/6.0; sg2[1][3] = 1.0/6.0;
      sg2[2][0] = 1.0/6.0; sg2[2][1] = 1.0/6.0; sg2[2][2] = 2.0/3.0; sg2[2][3] = 1.0/6.0;
    }
    else if (iord == 3) // four-point Gauss
    {
      sg2.resize(4);
      for (unsigned int i = 0; i < 4; ++i) sg2[i].resize(4);
      sg2[0][0] = 1.5505102572168219018027159252941e-01;
      sg2[0][1] = 1.7855872826361642311703513337422e-01;
      sg2[0][2] = 1.0 - sg2[0][0] - sg2[0][1];
      sg2[0][3] = 1.5902069087198858469718450103758e-01;

      sg2[1][0] = 6.4494897427831780981972840747059e-01;
      sg2[1][1] = 7.5031110222608118177475598324603e-02;
      sg2[1][2] = 1.0 - sg2[1][0] - sg2[1][1];
      sg2[1][3] = 9.0979309128011415302815498962418e-02;

      sg2[2][0] = 1.5505102572168219018027159252941e-01;
      sg2[2][1] = 6.6639024601470138670269327409637e-01;
      sg2[2][2] = 1.0 - sg2[2][0] - sg2[2][1];
      sg2[2][3] = 1.5902069087198858469718450103758e-01;

      sg2[3][0] = 6.4494897427831780981972840747059e-01;
      sg2[3][1] = 2.8001991549907407200279599420481e-01;
      sg2[3][2] = 1.0 - sg2[3][0] - sg2[3][1];
      sg2[3][3] = 9.0979309128011415302815498962418e-02;
    }
    else if (iord == 4) // six-point Guass
    {
      const unsigned int n_wts = 2;
      const Real wts[n_wts] =
        {
          1.1169079483900573284750350421656140e-01L,
          5.4975871827660933819163162450105264e-02L
        };

      const Real a[n_wts] =
        {
          4.4594849091596488631832925388305199e-01L,
          9.1576213509770743459571463402201508e-02L
        };

      const Real b[n_wts] = {0., 0.}; // not used
      const unsigned int permutation_ids[n_wts] = {3, 3};

      std::vector<Point> points; 
      std::vector<Real> weights;
      dunavant_rule2(wts, a, b, permutation_ids, n_wts, points, weights); // 6 total points 

      sg2.resize(6);
      for (unsigned int i = 0; i < 6; ++i)
        sg2[i].resize(4);
      for (unsigned int i = 0; i < 6; ++i)
      {
        sg2[i][0] = points[i](0);
        sg2[i][1] = points[i](1);
        sg2[i][2] = 1.0 - points[i](0) - points[i](1);
        sg2[i][3] = weights[i];
      }
    }
    else
      mooseError("Invalid quadrature order = " + Moose::stringify(iord) + " for triangle elements");
  }
  else
    mooseError("Invalid 2D element type");
}

void static wissmannPoints(unsigned int nqp, std::vector<std::vector<Real> > &wss)
{
  if (nqp == 6)
  {
    wss.resize(6);
    for (unsigned int i = 0; i < 6; ++i) wss[i].resize(3);
    wss[0][0] = 0.0;
    wss[0][1] = 0.0;
    wss[0][2] = 1.1428571428571428;

    wss[1][0] = 0.0;
    wss[1][1] = 9.6609178307929590e-01;
    wss[1][2] = 4.3956043956043956e-01;

    wss[2][0] = 8.5191465330460049e-01;
    wss[2][1] = 4.5560372783619284e-01;
    wss[2][2] = 5.6607220700753210e-01;

    wss[3][0] = -wss[2][0];
    wss[3][1] = wss[2][1];
    wss[3][2] = wss[2][2];

    wss[4][0] = 6.3091278897675402e-01;
    wss[4][1] = -7.3162995157313452e-01;
    wss[4][2] = 6.4271900178367668e-01;

    wss[5][0] = -wss[4][0];
    wss[5][1] = wss[4][1];
    wss[5][2] = wss[4][2];
  }
  else
    mooseError("Unknown Wissmann quadrature type");
}

void static shapeFunc2D(unsigned int nen, std::vector<Real> &ss, std::vector<Point> &xl,
                        std::vector<std::vector<Real> > &shp, Real &xsj, bool natl_flg)
{
  // Get shape functions and derivatives
  Real s[4] = {-0.5, 0.5, 0.5, -0.5};
  Real t[4] = {-0.5, -0.5, 0.5, 0.5};

  if (nen == 4) // quad element
  {
    Real xs[2][2] = {{0.0,0.0},{0.0,0.0}};
    Real sx[2][2] = {{0.0,0.0},{0.0,0.0}};
    for (unsigned int i = 0; i < 4; ++i)
    {
      shp[i][2] = (0.5 + s[i]*ss[0])*(0.5 + t[i]*ss[1]);
      shp[i][0] = s[i]*(0.5 + t[i]*ss[1]);
      shp[i][1] = t[i]*(0.5 + s[i]*ss[0]);
    }
    for (unsigned int i = 0; i < 2; ++i) //x, y
    {
      for (unsigned int j = 0; j < 2; ++j) // xi, eta
      {
        xs[i][j] = 0.0;
        for (unsigned int k = 0; k < nen; ++k)
          xs[i][j] += xl[k](i)*shp[k][j];
      }
    }
    xsj = xs[0][0]*xs[1][1] - xs[0][1]*xs[1][0]; // det(j)
    if (natl_flg == false) // get global derivatives
    {
      Real temp = 1.0/xsj;
      sx[0][0] =  xs[1][1]*temp; // inv(j)
      sx[1][1] =  xs[0][0]*temp;
      sx[0][1] = -xs[0][1]*temp;
      sx[1][0] = -xs[1][0]*temp;
      for (unsigned int i = 0; i < nen; ++i)
      {
        temp      = shp[i][0]*sx[0][0] + shp[i][1]*sx[1][0];
        shp[i][1] = shp[i][0]*sx[0][1] + shp[i][1]*sx[1][1];
        shp[i][0] = temp;
      }
    }
  }
  else if (nen == 3) // triangle element
  {
    // x1*(y2 - y3) + x2*(y3 - y1) + x3*(y1 - y2)
    xsj  = xl[0](0)*(xl[1](1)-xl[2](1)) + xl[1](0)*(xl[2](1)-xl[0](1)) + xl[2](0)*(xl[0](1)-xl[1](1));
    Real xsjr = 1.0;
    if (xsj != 0.0)
        xsjr = 1.0/xsj;
    // xsj  *= 0.5; // we do not have this 0.5 here because in stdQuad2D the sum of all weights in trig is 0.5
    shp[0][2] = ss[0];
    shp[1][2] = ss[1];
    shp[2][2] = ss[2];
    if (natl_flg == false) // need global drivatives
    {
      shp[0][0] = (xl[1](1) - xl[2](1))*xsjr;
      shp[0][1] = (xl[2](0) - xl[1](0))*xsjr;
      shp[1][0] = (xl[2](1) - xl[0](1))*xsjr;
      shp[1][1] = (xl[0](0) - xl[2](0))*xsjr;
      shp[2][0] = (xl[0](1) - xl[1](1))*xsjr;
      shp[2][1] = (xl[1](0) - xl[0](0))*xsjr;
    }
    else
    {
      shp[0][0] = 1.0;
      shp[0][1] = 0.0;
      shp[1][0] = 0.0;
      shp[1][1] = 1.0;
      shp[2][0] = -1.0;
      shp[2][1] = -1.0;
    }
  }
  else
    mooseError("ShapeFunc2D only works for linear quads and trigs!");
}

double static inline r8vec_norm(int n, double a[])
{
//    John Burkardt geometry.cpp
  double v = 0.0;
  for (unsigned int i = 0; i < n; ++i)
    v = v + a[i] * a[i];
  v = std::sqrt(v);
  return v;
}

void static inline r8vec_copy(int n, double a1[], double a2[])
{
//    John Burkardt geometry.cpp
  for (unsigned int i = 0; i < n; ++i)
    a2[i] = a1[i];
  return;
}

bool static inline r8vec_eq (int n, double a1[], double a2[])
{
//    John Burkardt geometry.cpp
  for (unsigned int i = 0; i < n; ++i)
    if ( a1[i] != a2[i] )
      return false;
  return true;
}

double static inline r8vec_dot_product(int n, double a1[], double a2[])
{
//    John Burkardt geometry.cpp
  double value = 0.0;
  for (unsigned int i = 0; i < n; ++i)
    value += a1[i] * a2[i];
  return value;
}

bool static inline isInsideEdge(Point p1, Point p2, Point p)
{
  Point p1_to_p2 = p2 - p1;
  Point p_to_p1 = p1 - p;
  Point p_to_p2 = p2 - p;
  Real dotp1 = p_to_p1*p1_to_p2;
  Real dotp2 = p_to_p2*p1_to_p2;
  if (dotp1*dotp2 <= 0.0)
    return true;
  else
    return false;
}

Real static inline getRelativePosition(Point p1, Point p2, Point p)
{
  // get the relative position of p from p1
  Real full_len = std::sqrt((p2 - p1).size_sq());
  Real len_p1_p = std::sqrt((p - p1).size_sq());
  return len_p1_p/full_len;
}

bool static inline line_exp_is_degenerate_nd(int dim_num, double p1[], double p2[])
{
//    John Burkardt geometry.cpp
  bool value;
  value = r8vec_eq(dim_num, p1, p2);
  return value;
}


int static plane_normal_line_exp_int_3d(double pp[3], double normal[3], double p1[3], double p2[3], double pint[3])
{
//    John Burkardt geometry.cpp
//  Parameters:
//
//    Input, double PP[3], a point on the plane.
//
//    Input, double NORMAL[3], a normal vector to the plane.
//
//    Input, double P1[3], P2[3], two distinct points on the line.
//
//    Output, double PINT[3], the coordinates of a
//    common point of the plane and line, when IVAL is 1 or 2.
//
//    Output, integer PLANE_NORMAL_LINE_EXP_INT_3D, the kind of intersection;
//    0, the line and plane seem to be parallel and separate;
//    1, the line and plane intersect at a single point;
//    2, the line and plane seem to be parallel and joined.
# define DIM_NUM 3

  double direction[DIM_NUM];
  int ival;
  double temp;
  double temp2;
//
//  Make sure the line is not degenerate.
  if (line_exp_is_degenerate_nd(DIM_NUM, p1, p2))
  {
    std::cerr << "\n";
    std::cerr << "PLANE_NORMAL_LINE_EXP_INT_3D - Fatal error!\n";
    std::cerr << "  The line is degenerate.\n";
    exit ( 1 );
  }
//
//  Make sure the plane normal vector is a unit vector.
  temp = r8vec_norm(DIM_NUM, normal);
  if (temp == 0.0)
  {
    std::cerr << "\n";
    std::cerr << "PLANE_NORMAL_LINE_EXP_INT_3D - Fatal error!\n";
    std::cerr << "  The normal vector of the plane is degenerate.\n";
    exit ( 1 );
  }

  for (unsigned int i = 0; i < DIM_NUM; ++i)
    normal[i] = normal[i]/temp;
//
//  Determine the unit direction vector of the line.
  for (unsigned int i = 0; i < DIM_NUM; ++i)
    direction[i] = p2[i] - p1[i];
  temp = r8vec_norm(DIM_NUM, direction);

  for (unsigned int i = 0; i < DIM_NUM; ++i)
    direction[i] = direction[i]/temp;
//
//  If the normal and direction vectors are orthogonal, then
//  we have a special case to deal with.
  if (r8vec_dot_product(DIM_NUM, normal, direction) == 0.0)
  {
    temp = 0.0;
    for (unsigned int i = 0; i < DIM_NUM; ++i)
      temp = temp + normal[i]*(p1[i] - pp[i]);

    if ( temp == 0.0 )
    {
      ival = 2;
      r8vec_copy(DIM_NUM, p1, pint);
    }
    else
    {
      ival = 0;
      for (unsigned int i = 0; i < DIM_NUM; ++i)
        pint[i] = 1.0e20; // dummy huge value
    }
    return ival;
  }
//
//  Determine the distance along the direction vector to the intersection point.
  temp = 0.0;
  for (unsigned int i = 0; i < DIM_NUM; ++i)
    temp = temp + normal[i]*(pp[i] - p1[i]);
  temp2 = 0.0;
  for (unsigned int i = 0; i < DIM_NUM; ++i)
    temp2 = temp2 + normal[i]*direction[i];

  ival = 1;
  for (unsigned int i = 0; i < DIM_NUM; ++i)
    pint[i] = p1[i] + temp*direction[i]/temp2;

  return ival;
# undef DIM_NUM
}



#endif
