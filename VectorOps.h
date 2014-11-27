
#ifndef DELAUNAY_VECTOR_OPS_H
#define DELAUNAY_VECTOR_OPS_H


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


// are p1 and p2 on the same side of line a-b ?
bool sameSide(const Point& p1, const Point& p2, const Point& a, const Point& b)
{
  const Point ab = {{b[0] - a[0], b[1] - a[1]}};
  const Point a1 = {{p1[0] - a[0], p1[1] - a[1]}};
  const Point a2 = {{p2[0] - a[0], p2[1] - a[1]}};

  // we're working in 2d, so cross product can be simplified:
  const double c1 = ab[0]*a1[1] - ab[1]*a1[0];
  const double c2 = ab[0]*a2[1] - ab[1]*a2[0];

  return (c1*c2 >= 0);
}


#endif // DELAUNAY_VECTOR_OPS_H
