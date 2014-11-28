
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

      // find and split the enclosing triangle
      const int index = findEnclosingTriangleIndex(p);
      splitTriangleAndReconnect(index, newVertex);

      // TODO: flip triangles for delaunay condition
    }

    size_t findEnclosingTriangleIndex(const Point& p) const
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

    // the split will generate new triangles. we store one by overwriting the
    // original triangles, and append the others at the end of the list
    void splitTriangleAndReconnect(const int index, const int newVertex)
    {
      // make a COPY not a reference because original triangle is destroyed in this function
      const Connectivity connection(connections_[index]);

      // shortcuts for sanity
      const int v0 = connection.vertex(0);
      const int v1 = connection.vertex(1);
      const int v2 = connection.vertex(2);
      const int n0 = connection.neighbor(0);
      const int n1 = connection.neighbor(1);
      const int n2 = connection.neighbor(2);

      // TODO: special cases when new point is on internal or external edge

      // add new sub-triangles
      const int subTriIndex1 = connections_.size();
      const int subTriIndex2 = subTriIndex1 + 1;
      connections_[index] = {v0, v1, newVertex, subTriIndex1, subTriIndex2, n2};
      connections_.emplace_back(v1, v2, newVertex, subTriIndex2, index, n0);
      connections_.emplace_back(v2, v0, newVertex, index, subTriIndex1, n1);

      // fix pre-existing neighbor triangles to point to new triangles
      connections_[n0].updateNeighbor(index, subTriIndex1);
      connections_[n1].updateNeighbor(index, subTriIndex2);
      // neighbor(2) still points to same index by 'replacement construction'

    }


  private:
    std::vector<Point> points_;
    std::vector<Connectivity> connections_;
    const double xmin_, xmax_;
    const double ymin_, ymax_;
};

#endif // DELAUNAY_DELAUNAYTRI_H
