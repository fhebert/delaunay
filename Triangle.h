
#ifndef INTERP_TRIANGLE_H
#define INTERP_TRIANGLE_H


#include "Point.h"
#include "VectorOps.h"

#include <string>


struct Triangle {

  Triangle(const Point& a, const Point& b, const Point& c)
    : a_(a), b_(b), c_(c) {}

  // TODO: make const, but this breaks copy c'tors used in the process of
  // deleting triangles from the std::vector (which i want to remove anyway)
  Point a_, b_, c_;
};


bool pointInTriangle(const Point& p, const Triangle& t)
{
  return (sameSide(p,t.a_, t.b_,t.c_)
      and sameSide(p,t.b_, t.a_,t.c_)
      and sameSide(p,t.c_, t.a_,t.b_)) ? true : false;
}


std::string ToString(const Triangle& tri)
{
  using std::to_string;
  return to_string(tri.a_[0]) + "  " + to_string(tri.a_[1])
    + "\n" + to_string(tri.b_[0]) + "  " + to_string(tri.b_[1])
    + "\n" + to_string(tri.c_[0]) + "  " + to_string(tri.c_[1])
    + "\n" + to_string(tri.a_[0]) + "  " + to_string(tri.a_[1]) + "\n";
}



#endif // INTERP_TRIANGLE_H
