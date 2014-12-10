
#include "DelaunayTri.h"
#include "Point.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <random>
#include <vector>


int main() {

  // parameters of the triangulation
  const double xmin = 0;
  const double xmax = 1;
  const double ymin = 0;
  const double ymax = 1;
  const double resolution = 0.15;

  // initialize triangulation
  DelaunayTri delaunay(xmin, xmax, ymin, ymax);

  // create initial point distribution
  std::vector<Point> points;
  // edge points
  const int nx = fabs(ceil( (xmax-xmin) / resolution));
  const int ny = fabs(ceil( (ymax-ymin) / resolution));
  for (int ix = 1; ix < nx-1; ++ix) {
    const double dx = (xmax-xmin) * ix/(nx-1);
    points.push_back({{xmin + dx, ymin}});
    points.push_back({{xmin + dx, ymax}});
  }
  for (int iy = 1; iy < ny-1; ++iy) {
    const double dy = (ymax-ymin) * iy/(ny-1);
    points.push_back({{xmin, ymin + dy}});
    points.push_back({{xmax, ymin + dy}});
  }
  // spiral on the inside
  for (int i=0; i<1000; ++i) {
    const double t = i;
    const double x = (xmin+xmax)/2.0 + 0.2*resolution * exp(0.01*t) * cos(0.4*t);
    const double y = (ymin+ymax)/2.0 + 0.2*resolution * exp(0.01*t) * sin(0.4*t);
    if (x < xmin or x > xmax or y < ymin or y > ymax) continue;
    //if (4*(x*x + y*y) > ((xmax-xmin)*
    points.push_back({{x, y}});
  }
  // randomize
  {
    std::random_device rd;
    std::mt19937 mt(rd());
    std::shuffle(points.begin(), points.end(), mt);
  }

  bool status = delaunay.addPoints(points);
  assert(status and "couldn't add points?");

  delaunay.writeToFile("triangulation.data");
}
