

#include "DelaunayTri.h"
#include "Point.h"

#include <cassert>
#include <vector>


int main() {

  // initialize triangulation on (0,1) x (0,1)
  DelaunayTri delaunay(0, 1, 0, 1);

  std::vector<Point> points;
  points.push_back({{0.2, 0.9}});
  points.push_back({{0.1, 0.3}});
  points.push_back({{0.7, 0.4}});
  points.push_back({{0.8, 0.2}});
  points.push_back({{0.5, 0.6}});

  bool status = delaunay.addPoints(points);
  assert(status and "couldn't add points?");

  delaunay.WriteToFile("triangulation.data");
}
