
#include "DelaunayTri.h"
#include "Interpolator.h"
#include "Point.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>


namespace {
  std::vector<Point> pointsToAdd(
      const double xmin, const double xmax,
      const double ymin, const double ymax)
  {
    std::vector<Point> points;

    // edge points
    const int resolution = 16;
    const int nx = fabs(ceil( (xmax-xmin) * resolution));
    const int ny = fabs(ceil( (ymax-ymin) * resolution));
    for (int ix = 1; ix < nx-1; ++ix) {
      const double dx = (xmax-xmin) * ix/(nx-1.0);
      points.push_back({{xmin + dx, ymin}});
      points.push_back({{xmin + dx, ymax}});
    }
    for (int iy = 1; iy < ny-1; ++iy) {
      const double dy = (ymax-ymin) * iy/(ny-1.0);
      points.push_back({{xmin, ymin + dy}});
      points.push_back({{xmax, ymin + dy}});
    }

    // interior points
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<> dist(0,1);
    for (int i = 0; i < nx*ny; ++i) {
      const double x = xmin + dist(mt) * (xmax-xmin);
      const double y = ymin + dist(mt) * (ymax-ymin);
      points.push_back({{x, y}});
    }
    return points;
  }

  double f(const double x, const double y) {
    return x*x*y;
  }

  double g(const double x, const double y) {
    return exp(-5 * ((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5)));
  }

  double h(const double x, const double y) {
    return sin(6*x) * cos(9*y);
  }
}



int main() {
  // parameters of the triangulation
  const double xmin = 0;
  const double xmax = 1;
  const double ymin = 0;
  const double ymax = 1;

  // construct triangulation
  DelaunayTri delaunay(xmin, xmax, ymin, ymax);
  {
    std::vector<Point> points = pointsToAdd(xmin, xmax, ymin, ymax);
    std::random_device rd;
    std::mt19937 mt(rd());
    std::shuffle(points.begin(), points.end(), mt);

    bool status = delaunay.addPoints(points);
    assert(status and "couldn't add points?");
    delaunay.writeToFile("triangulation.data");
  }

  // three functions to be interpolated
  // 1) x^2 * y
  // 2) exp(-5 (x-0.5)^2 + (y-0.5)^2)
  // 3) sin(6x) * cos(9y)
  const int numfns = 3;
  const size_t numpts = delaunay.getNumPoints(); // not the same as points.size() !!
  std::vector<std::vector<double>> data(numfns, std::vector<double>(numpts, 0));
  for (size_t i=0; i<numpts; ++i) {
    const double x = delaunay.point(i)[0];
    const double y = delaunay.point(i)[1];
    data[0][i] = f(x,y);
    data[1][i] = g(x,y);
    data[2][i] = h(x,y);
  }

  // grid to interpolate the data onto
  const int intres = 50;
  std::vector<Point> intgrid;
  for (int iy=0; iy<intres; ++iy) {
    for (int ix=0; ix<intres; ++ix) {
      const double x = xmin + (ix/(intres-1.0)) * (xmax-xmin);
      const double y = ymin + (iy/(intres-1.0)) * (ymax-ymin);
      intgrid.push_back({{x,y}});
    }
  }

  // interpolate
  const Interpolator interp(delaunay, data);
  std::vector<std::vector<double>> intdata(numfns, std::vector<double>(intres*intres, 0));
  for (int i=0; i<intres*intres; ++i) {
    const std::vector<double> intpt = interp.evaluateAt(intgrid[i]);
    intdata[0][i] = intpt[0];
    intdata[1][i] = intpt[1];
    intdata[2][i] = intpt[2];
  }

  // print to file
  std::ofstream outfile("interpolated.data");
  for (int i=0; i<intres*intres; ++i) {
    const double x = intgrid[i][0];
    const double y = intgrid[i][1];
    outfile << x << "  " << y << "  "
            << intdata[0][i] << "  " << intdata[1][i] << "  " << intdata[2][i] << "  "
            << f(x,y) << "  " << g(x,y) << "  " << h(x,y) << "\n";
  }
  outfile.close();

}
