
#ifndef DELAUNAY_DELAUNAYTRI_H
#define DELAUNAY_DELAUNAYTRI_H

#include "Point.h"
#include "Triangle.h"
#include "VectorOps.h"
#include "Utils.h"

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

      // for each point, add its triangles to the connection list
      trianglesWithPoint_.push_back({{0,3}});
      trianglesWithPoint_.push_back({{0,1}});
      trianglesWithPoint_.push_back({{1,2}});
      trianglesWithPoint_.push_back({{2,3}});
      trianglesWithPoint_.push_back({{0,1,2,3}});
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


    void writeToFile(const std::string& filename) const
    {
      std::ofstream outfile(filename);
      for (const auto& tri : triangles_) {
        if (tri.isLeaf())
          outfile << tri.toString() << "\n";
      }
      outfile.close();
    }


    std::vector<int> getConnectedPoints(const int pt) const
    {
      std::vector<int> result;
      for (const int tri : trianglesWithPoint_[pt]) {
        for (const int vert : triangles_[tri].vertices()) {
          if (vert != pt) {
            result.push_back(vert);
          }
        }
      }
      erase_dupes(result);
      return result;
    }


    std::tuple<int,int> findEnclosingTriangleIndex(const Point& p) const
    {
      int tri = -1;
      int edge = -1;
      // find root-level triangle that encloses point
      for (size_t i=0; i<4; ++i) {
        bool inside;
        std::tie(inside, edge) = triangles_[i].isPointInside(p);
        if (inside) {
          tri = i;
          break;
        }
      }
      assert(tri >= 0 and "point was not enclosed by root-level triangles");

      // iterate through children to find leaf triangle enclosing point
      while (not triangles_[tri].isLeaf()) {
        for (const auto& child : triangles_[tri].children()) {
          bool inside;
          std::tie(inside, edge) = triangles_[child].isPointInside(p);
          if (inside) {
            tri = child;
            break;
          }
        }
      }
      return std::make_tuple(tri, edge);
    }


    const Triangle& triangle(const size_t i) const {return triangles_[i];}
    const Point& point(const size_t i) const {return points_[i];}
    size_t getNumPoints() const {return points_.size();}


  private:
    void insertPointAndRetriangulate(const Point& p)
    {
      // add point to list
      const int newVertex = points_.size();
      points_.push_back(p);
      trianglesWithPoint_.push_back({{}});

      // find and split the enclosing triangle
      int index, edge;
      std::tie(index, edge) = findEnclosingTriangleIndex(p);
      splitTriangle(index, edge, newVertex);
    }


    void splitTriangle(const int index, const int edge, const int newVertex)
    {
      // shortcuts for sanity
      assert(index >= 0 and "can't index an invalid triangle");
      const int v0 = triangles_[index].vertex(0);
      const int v1 = triangles_[index].vertex(1);
      const int v2 = triangles_[index].vertex(2);
      const int n0 = triangles_[index].neighbor(0);
      const int n1 = triangles_[index].neighbor(1);
      const int n2 = triangles_[index].neighbor(2);

      // TODO: special case for internal edge
      // for now, assert that only external edges are allowed
      assert(edge == -1 or triangles_[index].neighbor(edge) == -1);

      // special case when new point is on external edge
      if (edge >= 0) { // then newVertex lies on edge
        // going around the triangle in the positive (CCW) direction, starting
        // with the edge 'edge' (which contains newVertex), permutation holds
        // the correct ordering of the triangle's local vertex indices {0,1,2}
        const std::array<int,3> permutation = {{(edge+2)%3, edge, (edge+1)%3}};

        // from the permutation, get the global indices of points and nghbrs
        const std::array<int,3> pts = {{triangles_[index].vertex(permutation[0]),
                                        triangles_[index].vertex(permutation[1]),
                                        triangles_[index].vertex(permutation[2])}};
        const std::array<int,3> nbs = {{triangles_[index].neighbor(permutation[0]),
                                        -1, // by construction, this is the ext. edge
                                        triangles_[index].neighbor(permutation[2])}};
        const int child0 = triangles_.size();
        const int child1 = child0 + 1;

        triangles_[index].setChildren(child0, child1);
        triangles_.emplace_back(points_, pts[0], pts[1], newVertex, child1, nbs[1], nbs[2]);
        triangles_.emplace_back(points_, pts[1], pts[2], newVertex, nbs[1], child0, nbs[0]);

        trianglesWithPoint_[pts[0]].push_back(child0);
        trianglesWithPoint_[pts[1]].push_back(child0);
        trianglesWithPoint_[pts[1]].push_back(child1);
        trianglesWithPoint_[pts[2]].push_back(child1);
        trianglesWithPoint_[newVertex].push_back(child0);
        trianglesWithPoint_[newVertex].push_back(child1);
        if (nbs[0] >= 0) triangles_[nbs[0]].updateNeighbor(index, child1);
        if (nbs[2] >= 0) triangles_[nbs[2]].updateNeighbor(index, child0);

        delaunayFlip(child0, newVertex);
        delaunayFlip(child1, newVertex);
        return;
      }

      // add new sub-triangles
      // i'th child is adjacent to i'th neighbor, opposite from i'th vertex
      const int child0 = triangles_.size();
      const int child1 = child0 + 1;
      const int child2 = child0 + 2;
      triangles_.emplace_back(points_, v1, v2, newVertex, child1, child2, n0);
      triangles_.emplace_back(points_, v2, v0, newVertex, child2, child0, n1);
      triangles_.emplace_back(points_, v0, v1, newVertex, child0, child1, n2);
      trianglesWithPoint_[v0].push_back(child1);
      trianglesWithPoint_[v0].push_back(child2);
      trianglesWithPoint_[v1].push_back(child0);
      trianglesWithPoint_[v1].push_back(child2);
      trianglesWithPoint_[v2].push_back(child0);
      trianglesWithPoint_[v2].push_back(child1);
      trianglesWithPoint_[newVertex].push_back(child0);
      trianglesWithPoint_[newVertex].push_back(child1);
      trianglesWithPoint_[newVertex].push_back(child2);

      // fix pre-existing neighbor triangles to point to new triangles
      if (n0 >= 0) {
        triangles_[n0].updateNeighbor(index, child0);
      }
      if (n1 >= 0) {
        triangles_[n1].updateNeighbor(index, child1);
      }
      if (n2 >= 0) {
        triangles_[n2].updateNeighbor(index, child2);
      }

      // point triangle to children
      triangles_[index].setChildren(child0, child1, child2);

      // check delaunay
      delaunayFlip(child0, newVertex);
      delaunayFlip(child1, newVertex);
      delaunayFlip(child2, newVertex);
    }


    void delaunayFlip(const int tri, const int keyPoint)
    {
      const int oppTri = triangles_[tri].neighborAcrossGlobalPoint(keyPoint);
      if (oppTri == -1) return;

      // find point of oppTri that is across from keyPoint
      int oppPoint = -1;
      for (int pt=0; pt<3; ++pt) {
        if (triangles_[oppTri].neighbor(pt) == tri) {
          oppPoint = triangles_[oppTri].vertex(pt);
          break;
        }
      }

      const double keyAngle = triangles_[tri].angleAtPoint(keyPoint);
      const double oppAngle = triangles_[oppTri].angleAtPoint(oppPoint);

      // check delaunay condition
      if (keyAngle + oppAngle > M_PI) {
        int tri1, tri2;
        std::tie(tri1, tri2) = swapTwoTriangles(tri, keyPoint, oppTri, oppPoint);
        delaunayFlip(tri1, keyPoint);
        delaunayFlip(tri2, keyPoint);
      }
    }


    std::tuple<int, int> swapTwoTriangles(const int tri1, const int pt1,
        const int tri2, const int pt2)
    {
      const int child0 = triangles_.size();
      const int child1 = child0 + 1;

      // TODO a more elegant way for the entire thing
      int localIndex1 = -1;
      for (int i=0; i<3; ++i) {
        if (triangles_[tri1].vertex(i) == pt1) {
          localIndex1 = i;
          break;
        }
      }
      const int shared1 = triangles_[tri1].vertex((localIndex1+1)%3);
      const int shared2 = triangles_[tri1].vertex((localIndex1+2)%3);

      const int n1_1 = triangles_[tri1].neighborAcrossGlobalPoint(shared1);
      const int n1_2 = triangles_[tri1].neighborAcrossGlobalPoint(shared2);
      const int n2_1 = triangles_[tri2].neighborAcrossGlobalPoint(shared1);
      const int n2_2 = triangles_[tri2].neighborAcrossGlobalPoint(shared2);

      triangles_.emplace_back(points_, pt1, pt2, shared2, n2_1, n1_1, child1);
      triangles_.emplace_back(points_, pt1, shared1, pt2, n2_2, child0, n1_2);
      trianglesWithPoint_[pt1].push_back(child0);
      trianglesWithPoint_[pt1].push_back(child1);
      trianglesWithPoint_[pt2].push_back(child0);
      trianglesWithPoint_[pt2].push_back(child1);
      trianglesWithPoint_[shared1].push_back(child1);
      trianglesWithPoint_[shared2].push_back(child0);

      if (n1_1 >= 0) triangles_[n1_1].updateNeighbor(tri1, child0);
      if (n1_2 >= 0) triangles_[n1_2].updateNeighbor(tri1, child1);
      if (n2_1 >= 0) triangles_[n2_1].updateNeighbor(tri2, child0);
      if (n2_2 >= 0) triangles_[n2_2].updateNeighbor(tri2, child1);

      triangles_[tri1].setChildren(child0, child1);
      triangles_[tri2].setChildren(child0, child1);
      return std::make_tuple(child0, child1);
    }


  private:
    std::vector<Point> points_;
    std::vector<Triangle> triangles_;
    std::vector<std::vector<int>> trianglesWithPoint_;
    const double xmin_, xmax_;
    const double ymin_, ymax_;
};


#endif // DELAUNAY_DELAUNAYTRI_H
