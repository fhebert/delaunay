
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
      assert(index >= 0 and "can't index an invalid triangle");

      // most likely scenario: newVertex is contained in triangle "index"
      // then, edge has its default value of -1
      if (edge == -1) {
        // i'th child is adjacent to i'th neighbor, opposite from i'th vertex
        const int v0 = triangles_[index].vertex(0);
        const int v1 = triangles_[index].vertex(1);
        const int v2 = triangles_[index].vertex(2);
        const int n0 = triangles_[index].neighbor(0);
        const int n1 = triangles_[index].neighbor(1);
        const int n2 = triangles_[index].neighbor(2);
        const int child0 = triangles_.size();
        const int child1 = child0 + 1;
        const int child2 = child0 + 2;

        constructTriangle(child0, v1, v2, newVertex, child1, child2, n0);
        constructTriangle(child1, v2, v0, newVertex, child2, child0, n1);
        constructTriangle(child2, v0, v1, newVertex, child0, child1, n2);

        triangles_[index].setChildren(child0, child1, child2);
        if (n0 >= 0) triangles_[n0].updateNeighbor(index, child0);
        if (n1 >= 0) triangles_[n1].updateNeighbor(index, child1);
        if (n2 >= 0) triangles_[n2].updateNeighbor(index, child2);

        // check delaunay
        delaunayFlip(child0, newVertex);
        delaunayFlip(child1, newVertex);
        delaunayFlip(child2, newVertex);

      }
      // special case for external edges, make two triangles in one triangle
      else if (edge >= 0 and triangles_[index].neighbor(edge) == -1) {
        // going around the triangle in the positive (CCW) direction, starting
        // with the edge 'edge' (which contains newVertex), permutation holds
        // the correct ordering of the triangle's local vertex indices {0,1,2}
        const std::array<int,3> permutation = {{(edge+2)%3, edge, (edge+1)%3}};

        // from the permutation, get the global indices of points and nghbrs
        const int v0 = triangles_[index].vertex(permutation[0]);
        const int v1 = triangles_[index].vertex(permutation[1]);
        const int v2 = triangles_[index].vertex(permutation[2]);
        const int n0 = triangles_[index].neighbor(permutation[0]);
        const int n1 = -1; // by construction; this is the extern. edge
        const int n2 = triangles_[index].neighbor(permutation[2]);
        const int child0 = triangles_.size();
        const int child1 = child0 + 1;

        constructTriangle(child0, v0, v1, newVertex, child1, n1, n2);
        constructTriangle(child1, v1, v2, newVertex, n1, child0, n0);

        triangles_[index].setChildren(child0, child1);
        if (n0 >= 0) triangles_[n0].updateNeighbor(index, child1);
        if (n2 >= 0) triangles_[n2].updateNeighbor(index, child0);

        delaunayFlip(child0, newVertex);
        delaunayFlip(child1, newVertex);
      }
      // special case for internal edges, make four triangles in two triangles
      else {
        // index of "other" triangle to split
        const int neighb = triangles_[index].neighbor(edge);

        // permutation of "this" triangle
        const std::array<int,3> permutation = {{(edge+2)%3, edge, (edge+1)%3}};

        // from the permutation, get the global indices of points and nghbrs
        const int v0 = triangles_[index].vertex(permutation[0]);
        const int v1 = triangles_[index].vertex(permutation[1]);
        const int v2 = triangles_[index].vertex(permutation[2]);
        const int v3 = triangles_[neighb].pointOppositeFromNeighbor(index);
        const int np0 = triangles_[index].neighbor(permutation[0]);
        const int np2 = triangles_[index].neighbor(permutation[2]);
        const int nn0 = triangles_[neighb].neighbor(permutation[0]);
        const int nn2 = triangles_[neighb].neighbor(permutation[2]);
        const int child0 = triangles_.size();
        const int child1 = child0 + 1;
        const int child2 = child0 + 2;
        const int child3 = child0 + 3;

        constructTriangle(child0, v0, v1, newVertex, child1, child3, np2);
        constructTriangle(child1, v1, v2, newVertex, child2, child0, np0);
        constructTriangle(child2, v2, v3, newVertex, child3, child1, nn0);
        constructTriangle(child3, v3, v0, newVertex, child0, child2, nn2);

        triangles_[index].setChildren(child0, child1);
        if (np0 >= 0) triangles_[np0].updateNeighbor(index, child1);
        if (np2 >= 0) triangles_[np2].updateNeighbor(index, child0);

        triangles_[neighb].setChildren(child2, child3);
        if (nn0 >= 0) triangles_[nn0].updateNeighbor(index, child2);
        if (nn2 >= 0) triangles_[nn2].updateNeighbor(index, child3);

        delaunayFlip(child0, newVertex);
        delaunayFlip(child1, newVertex);
        delaunayFlip(child2, newVertex);
        delaunayFlip(child3, newVertex);
      }

      return;
    }


    void constructTriangle(const int index,
        const int p0, const int p1, const int p2,
        const int n0, const int n1, const int n2)
    {
      triangles_.emplace_back(points_, p0, p1, p2, n0, n1, n2);
      trianglesWithPoint_[p0].push_back(index);
      trianglesWithPoint_[p1].push_back(index);
      trianglesWithPoint_[p2].push_back(index);
    }


    void delaunayFlip(const int tri, const int keyPoint)
    {
      const int oppTri = triangles_[tri].neighborAcrossGlobalPoint(keyPoint);
      if (oppTri == -1) return; // it's an external edge

      // find point of oppTri that is across from tri
      const int oppPoint = triangles_[oppTri].pointOppositeFromNeighbor(tri);

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
