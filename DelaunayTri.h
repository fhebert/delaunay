
#ifndef DELAUNAY_DELAUNAYTRI_H
#define DELAUNAY_DELAUNAYTRI_H


#include "Connectivity.h"
#include "Point.h"
#include "Triangle.h"
#include "VectorOps.h"

#include <cassert>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>



// following Lee & Schachter 1980
// "Two algorithms for constructing a delaunay trianguation"
// specifically, uses the iterative method in a rectangular region
class DelaunayTri {

  public:
    // initialize a mostly-empty triangulation
    DelaunayTri(const double xmin, const double xmax, const double ymin, const double ymax)
    : xmin_(xmin), xmax_(xmax), ymin_(ymin), ymax_(ymax)
    {
      // initialize the corner and center points
      points_.push_back({{xmin_, ymin_}});
      points_.push_back({{xmax_, ymin_}});
      points_.push_back({{xmax_, ymax_}});
      points_.push_back({{xmin_, ymax_}});
      points_.push_back({{(xmin_+xmax_)/2.0, (ymin_+ymax)/2.0}});

      // explicitly connect triangles in this small list
      connections_.emplace_back(0, 1, 4, 1, 3, -1);
      connections_.emplace_back(1, 2, 4, 2, 0, -1);
      connections_.emplace_back(2, 3, 4, 3, 1, -1);
      connections_.emplace_back(3, 0, 4, 0, 2, -1);
    }


    bool addPoints(const std::vector<Point>& MorePoints)
    {
      // check there are no invalid points
      for (const auto& newpoint : MorePoints) {
        // check it's in the rectangle
        if (newpoint[0] < xmin_ or newpoint[0] > xmax_
            or newpoint[1] < ymin_ or newpoint[1] > ymax_) {
          return false;
        }
        // check it's not a duplicate
        for (const auto& oldpoint : points_) {
          if (distance(newpoint, oldpoint) < 1e-6) {
            return false;
          }
        }
      }

      // iteratively add each point
      for (const auto& newpoint : MorePoints) {
        insertPointAndRetriangulate(newpoint);
      }
      return true;
    }


    void WriteToFile(const std::string& filename) const
    {
      std::ofstream outfile(filename);
      for (const auto& c : connections_) {
        const Triangle tri(c, points_);
        outfile << tri.toString() << "\n";
      }
      outfile.close();
    }


  private:
    void insertPointAndRetriangulate(const Point& p)
    {
      // add point to list
      const int newVertex = points_.size();
      points_.push_back(p);

      // start by finding the indices of the new triangles
      // the split will generate new triangles. we store one by overwriting the
      // original triangles, and append the others at the end of the list
      const int triIndex0 = findEnclosingTriangleIndex(p);
      const int triIndex1 = connections_.size();

      // make a COPY not a reference because original triangle is destroyed in this function
      const Connectivity c(connections_[triIndex0]);

      // TODO: special cases when new point is on internal or external edge

      // add new sub-triangles
      const int triIndex2 = triIndex1 + 1;
      connections_[triIndex0] = {c.vertex(0), c.vertex(1), newVertex, triIndex1, triIndex2, c.neighbor(2)};
      connections_.emplace_back(c.vertex(1), c.vertex(2), newVertex, triIndex2, triIndex0, c.neighbor(0));
      connections_.emplace_back(c.vertex(2), c.vertex(0), newVertex, triIndex0, triIndex1, c.neighbor(1));

      // fix pre-existing neighbor triangles to point to new triangles
      connections_[c.neighbor(0)].updateNeighbor(triIndex0, triIndex1);
      connections_[c.neighbor(1)].updateNeighbor(triIndex0, triIndex2);
      // neighbor(2) still points to same index by 'replacement construction'

      // TODO: flip triangles for delaunay condition
    }

    size_t findEnclosingTriangleIndex(const Point& p)
    {
      // for now, brute-force by trying all leaf triangles
      // TODO: implement real method that tries a triangle and then moves towards
      //       the known p for the next guess...
      for (size_t i=0; i<connections_.size(); ++i) {
        const Triangle tri(connections_[i], points_);
        if (tri.isPointInside(p))
          return i;
      }
      assert(false and "should have found the enclosing triangle");
      return 0;
    }


  private:
    std::vector<Point> points_;
    std::vector<Connectivity> connections_;
    const double xmin_, xmax_;
    const double ymin_, ymax_;
};

#endif // DELAUNAY_DELAUNAYTRI_H
