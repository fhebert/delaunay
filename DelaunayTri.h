
#ifndef DELAUNAY_DELAUNAYTRI_H
#define DELAUNAY_DELAUNAYTRI_H


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
      if (tris[i].isPointInside(p))
        return i;
    }
    assert(false and "should have found the enclosing triangle");
    return 0;
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
      triangles_.push_back({points_, 0, 1, 4, 1, 3, -1});
      triangles_.push_back({points_, 1, 2, 4, 2, 0, -1});
      triangles_.push_back({points_, 2, 3, 4, 3, 1, -1});
      triangles_.push_back({points_, 3, 0, 4, 0, 2, -1});
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
        insertPoint(newpoint);
      }
      return true;
    }


    void WriteToFile(const std::string& filename) const
    {
      std::ofstream outfile(filename);
      for (auto tri : triangles_) {
        outfile << tri.toString() << "\n";
      }
      outfile.close();
    }


  private:
    void insertPoint(const Point& p)
    {
      // add point to list
      const int newPoint = points_.size();
      points_.push_back(p);

      // start by identifying the triangle that contains point
      const int index = enclosingTriangle(p, triangles_);
      const Triangle& tri = triangles_[index];

      // TODO: special cases when new point is on internal or external edge

      // add new sub-triangles
      const int size = triangles_.size();
      triangles_.push_back({points_, tri.vertex(0), tri.vertex(1), newPoint, size+1, size+2, tri.neighbor(2)});
      triangles_.push_back({points_, tri.vertex(1), tri.vertex(2), newPoint, size+2, size, tri.neighbor(0)});
      triangles_.push_back({points_, tri.vertex(2), tri.vertex(0), newPoint, size, size+1, tri.neighbor(1)});

      // erase the parent triangle from the list
      // TODO: avoid this inneficient array deletion
      // TODO: deal with triangle neighbors
      //tris.erase(tris.begin()+index);

      // TODO: flip triangles for delaunay condition
    }


  private:
    std::vector<Point> points_;
    std::vector<Triangle> triangles_;
    const double xmin_, xmax_;
    const double ymin_, ymax_;
};

#endif // DELAUNAY_DELAUNAYTRI_H
