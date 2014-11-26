
#ifndef INTERP_DELAUNAYTRI_H
#define INTERP_DELAUNAYTRI_H


#include "Point.h"
#include "Triangle.h"
#include "VectorOps.h"

#include <cassert>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>



namespace {

  size_t enclosingTriangle(const Point& p, const std::vector<Triangle>& tris)
  {
    // for now, brute-force by trying all triangles
    // TODO: implement real method that tries a triangle and then moves towards
    //       the known p for the next guess...
    for (size_t i=0; i<tris.size(); ++i) {
      if (pointInTriangle(p, tris[i]))
        return i;
    }
    assert(false and "should have found the enclosing triangle");
    return 0;
  }

  void insertPoint(const Point& p, std::vector<Triangle>& tris)
  {
    // start by identifying the triangle that contains point
    const int index = enclosingTriangle(p, tris);

    // add new sub-triangles
    tris.push_back({tris[index].a_, tris[index].b_, p});
    tris.push_back({tris[index].b_, tris[index].c_, p});
    tris.push_back({tris[index].c_, tris[index].a_, p});

    // erase the parent triangle from the list
    // TODO: avoid this inneficient array deletion
    // TODO: deal with triangle neighbors
    std::cout << "erasing triangle at index " << index << "\n"
      << ToString(tris[index]) << "\n";
    tris.erase(tris.begin()+index);

    // TODO: flip triangles for delaunay condition
  }

}



// following Lee & Schachter 1980
// "Two algorithms for constructing a delaunay trianguation"
// specifically, uses the iterative method in a rectangular region
class DelaunayTri {

  public:
    // initialize a mostly-empty triangulation
    DelaunayTri(
        const double xmin, const double xmax,
        const double ymin, const double ymax)
    : points_(), triangles_(),
      xmin_(xmin), xmax_(xmax),
      ymin_(ymin), ymax_(ymax)
    {
      // initialize the corner and center points
      points_.push_back({{xmin_, ymin_}});
      points_.push_back({{xmax_, ymin_}});
      points_.push_back({{xmax_, ymax_}});
      points_.push_back({{xmin_, ymax_}});
      points_.push_back({{(xmin_+xmax_)/2.0, (ymin_+ymax)/2.0}});

      // triangulate this small list
      // TODO: triangles need to know about their neighbors
      triangles_.push_back({points_[0], points_[1], points_[4]});
      triangles_.push_back({points_[1], points_[2], points_[4]});
      triangles_.push_back({points_[2], points_[3], points_[4]});
      triangles_.push_back({points_[3], points_[0], points_[4]});
    }


    bool addPoints(const std::vector<Point>& MorePoints)
    {
      // check there are no invalid points
      for (auto newpoint : MorePoints) {
        // check it's in the rectangle
        if (newpoint[0] < xmin_ or newpoint[0] > xmax_
            or newpoint[1] < ymin_ or newpoint[1] > ymax_) {
          return false;
        }
        // check it's not a duplicate
        for (auto oldpoint : points_) {
          if (distance(newpoint, oldpoint) < 1e-6) {
            return false;
          }
        }
      }

      // iteratively add each point
      for (auto newpoint : MorePoints) {
        insertPoint(newpoint, triangles_);
      }
      return true;
    }


    void WriteToFile(const std::string& filename) const
    {
      std::ofstream outfile(filename);
      for (auto tri : triangles_) {
        outfile << ToString(tri) << "\n";
      }
      outfile.close();
    }

  private:
    std::vector<Point> points_;
    std::vector<Triangle> triangles_;
    const double xmin_, xmax_;
    const double ymin_, ymax_;
};

#endif // INTERP_DELAUNAYTRI_H
