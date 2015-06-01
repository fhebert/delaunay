
#ifndef DELAUNAY_VECTOR_OPS_H
#define DELAUNAY_VECTOR_OPS_H

#include "Point.h"

#include <cmath>


double dotProduct(const Point& a, const Point& b)
{
  return a[0]*b[0] + a[1]*b[1];
}


double norm(const Point& a)
{
  return sqrt(dotProduct(a, a));
}


double distance(const Point& a, const Point& b)
{
  const Point diff = {{a[0] - b[0], a[1] - b[1]}};
  return norm(diff);
}


// area of triangle abc (points going around ccw) as computed by
// 1/2 * ((a->b) cross (a->c))
double orientedArea(const Point& a, const Point& b, const Point& c)
{
  const Point ab = {{b[0] - a[0], b[1] - a[1]}};
  const Point ac = {{c[0] - a[0], c[1] - a[1]}};
  return 0.5 * (ab[0]*ac[1] - ab[1]*ac[0]);
}


#endif // DELAUNAY_VECTOR_OPS_H
