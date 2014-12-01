
#ifndef DELAUNAY_TRIANGLE_H
#define DELAUNAY_TRIANGLE_H


#include "Point.h"
#include "VectorOps.h"
#include "Utils.h"

#include <array>
#include <cassert>
#include <string>
#include <tuple>
#include <vector>


class Triangle {

  public:
    Triangle(const std::vector<Point>& points,
        const int a, const int b, const int c, const int na, const int nb, const int nc)
      : isLeaf_(true),
      v0_(points[a]), v1_(points[b]), v2_(points[c]),
      vertices_({{a,b,c}}), neighbors_({{na,nb,nc}}), children_({{-1,-1,-1}}) {}


    int vertex(const int i) const {
      assert(i==0 or i==1 or i==2);
      return vertices_[i];
    }

    int neighbor(const int i) const {
      assert(i==0 or i==1 or i==2);
      return neighbors_[i];
    }

    void setChildren(const int ca, const int cb, const int cc=-1) {
      children_[0] = ca;
      children_[1] = cb;
      children_[2] = cc;
    }

    std::array<int, 3> children() const { return children_; }

    bool isLeaf() const { return children_[0] == -1; }

    void updateNeighbor(const int oldnbr, const int newnbr) {
      assert(contains(neighbors_, oldnbr));
      if (oldnbr == neighbors_[0]) {
        neighbors_[0] = newnbr;
      } else if (oldnbr == neighbors_[1]) {
        neighbors_[1] = newnbr;
      } else {
        neighbors_[2] = newnbr;
      }
    }

    std::tuple<bool,int> isPointInside(const Point& p) const {
      const double a01p = orientedArea(v0_, v1_, p);
      const double a12p = orientedArea(v1_, v2_, p);
      const double a20p = orientedArea(v2_, v0_, p);
      // triangles are set up to all have the same orientation (+), so we need
      // only verify that the sign is appropriate w.r.t. all sides
      const bool inside = (a01p >= 0.0 and a12p >= 0.0 and a20p >= 0.0);

      // if any area is too small, then treat the point as lying on edge
      int edge = -1;
      if (inside) {
        // threshold of (1e-12 * area) is meant to be small, but hte prefactor is arbitarily chosen
        const double threshold = 1e-12 * fabs(orientedArea(v0_, v1_, v2_));
        if (a01p < threshold)
          edge = 2;
        else if (a12p < threshold)
          edge = 0;
        else if (a20p < threshold)
          edge = 1;
      }

      return std::make_tuple(inside, edge);
    }

    std::string toString() const {
      using std::to_string;
      return to_string(v0_[0]) + "  " + to_string(v0_[1])
        + "\n" + to_string(v1_[0]) + "  " + to_string(v1_[1])
        + "\n" + to_string(v2_[0]) + "  " + to_string(v2_[1])
        + "\n" + to_string(v0_[0]) + "  " + to_string(v0_[1])
        + "\n";
    }


  private:
    // triangle defined by 3 vertices, has 3 neighbors
    // i'th neighbor is opposite from i'th vertex
    bool isLeaf_;
    const Point& v0_;
    const Point& v1_;
    const Point& v2_;
    const std::array<int, 3> vertices_;
    std::array<int, 3> neighbors_;
    std::array<int, 3> children_;
};


#endif // DELAUNAY_TRIANGLE_H
