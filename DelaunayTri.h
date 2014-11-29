
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
      triangles_.emplace_back(points_, 0, 1, 4, 1, 3, -1);
      triangles_.emplace_back(points_, 1, 2, 4, 2, 0, -1);
      triangles_.emplace_back(points_, 2, 3, 4, 3, 1, -1);
      triangles_.emplace_back(points_, 3, 0, 4, 0, 2, -1);
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
        writeToFile("temp.out");
      }
      return true;
    }


    void writeToFile(const std::string& filename) const
    {
      std::ofstream outfile(filename);
      for (const auto& tri : triangles_) {
        if (tri.isLeaf())
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
      splitTriangle(index, newVertex);

      // TODO: flip triangles for delaunay condition
    }


    size_t findEnclosingTriangleIndex(const Point& p) const
    {
      int tri = -1;
      for (size_t i=0; i<4; ++i) {
        if (triangles_[i].isPointInside(p)) {
          tri = i;
          break;
        }
      }

      while (not triangles_[tri].isLeaf()) {
        for (const auto& sub : triangles_[tri].subTriangles()) {
          if (triangles_[sub].isPointInside(p)) {
            tri = sub;
            break;
          }
        }
      }
      return tri;
    }


    void splitTriangle(const int index, const int newVertex)
    {
      // shortcuts for sanity
      const int v0 = triangles_[index].vertex(0);
      const int v1 = triangles_[index].vertex(1);
      const int v2 = triangles_[index].vertex(2);
      const int n0 = triangles_[index].neighbor(0);
      const int n1 = triangles_[index].neighbor(1);
      const int n2 = triangles_[index].neighbor(2);

      // TODO: special cases when new point is on internal or external edge

      // add new sub-triangles
      // i'th sub triangle is adjacent to i'th neighbor, opposite from i'th vertex
      const int subTri0 = triangles_.size();
      const int subTri1 = subTri0 + 1;
      const int subTri2 = subTri0 + 2;
      triangles_.emplace_back(points_, v1, v2, newVertex, subTri1, subTri2, n0);
      triangles_.emplace_back(points_, v2, v0, newVertex, subTri2, subTri0, n1);
      triangles_.emplace_back(points_, v0, v1, newVertex, subTri0, subTri1, n2);

      // fix pre-existing neighbor triangles to point to new triangles
      if (n0 >= 0) {
        triangles_[n0].updateNeighbor(index, subTri0);
      }
      if (n1 >= 0) {
        triangles_[n1].updateNeighbor(index, subTri1);
      }
      if (n2 >= 0) {
        triangles_[n2].updateNeighbor(index, subTri2);
      }

      // point triangle to children
      triangles_[index].setSubTriangles({subTri0, subTri1, subTri2});
    }



  private:
    std::vector<Point> points_;
    std::vector<Triangle> triangles_;
    const double xmin_, xmax_;
    const double ymin_, ymax_;
};

#endif // DELAUNAY_DELAUNAYTRI_H
