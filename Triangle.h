
#ifndef INTERP_TRIANGLE_H
#define INTERP_TRIANGLE_H


#include "Point.h"
#include "VectorOps.h"
#include "Utils.h"

#include <cassert>
#include <string>
#include <vector>


class Triangle {

  public:
    Triangle(const std::vector<Point>& points,
        const int a, const int b, const int c, const int na, const int nb, const int nc)
      : points_(points), vertices_({{a,b,c}}), neighbors_({{na,nb,nc}}) {}


    int neighbor(const int i) const {
      assert(i==0 or i==1 or i==2);
      return neighbors_[i];
    }

    int vertex(const int i) const {
      assert(i==0 or i==1 or i==2);
      return vertices_[i];
    }

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

    bool isPointInside(const Point& p) const {
      const Point& v0 = points_[vertices_[0]];
      const Point& v1 = points_[vertices_[1]];
      const Point& v2 = points_[vertices_[2]];
      return (sameSide(p,v0, v1,v2) and sameSide(p,v1, v0,v2) and sameSide(p,v2, v0,v1));
    }

    std::string toString() const {
      const Point& v0 = points_[vertices_[0]];
      const Point& v1 = points_[vertices_[1]];
      const Point& v2 = points_[vertices_[2]];
      using std::to_string;
      return to_string(v0[0]) + "  " + to_string(v0[1])
        + "\n" + to_string(v1[0]) + "  " + to_string(v1[1])
        + "\n" + to_string(v2[0]) + "  " + to_string(v2[1])
        + "\n" + to_string(v0[0]) + "  " + to_string(v0[1])
        + "\n";
    }

  private:
    // triangle defined by 3 vertices, has 3 neighbors
    // i'th neighbor is opposite from i'th vertex
    const std::vector<Point>& points_;
    const std::array<int, 3> vertices_;
    std::array<int, 3> neighbors_;
};


#endif // INTERP_TRIANGLE_H
