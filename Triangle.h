
#ifndef INTERP_TRIANGLE_H
#define INTERP_TRIANGLE_H


#include "Point.h"
#include "VectorOps.h"

#include <cassert>
#include <string>


class Triangle {

  public:
    Triangle(const int index, const Point& a, const Point& b, const Point& c,
        const int na, const int nb, const int nc)
      : index_(index), vertices_({{a,b,c}}), neighbors_({{na,nb,nc}}) {}


    int neighbor(const int i) const {
      assert(i==0 or i==1 or i==2);
      return neighbors_[i];
    }

    Point vertex(const int i) const {
      assert(i==0 or i==1 or i==2);
      return vertices_[i];
    }

    void updateNeighbor(const int oldnbr, const int newnbr) {
      assert(oldnbr==neighbors_[0] or oldnbr==neighbors_[1] or oldnbr==neighbors_[2]);
      if (oldnbr == neighbors_[0]) {
        neighbors_[0] = newnbr;
      } else if (oldnbr == neighbors_[1]) {
        neighbors_[1] = newnbr;
      } else {
        neighbors_[2] = newnbr;
      }
    }

    bool isPointEnclosed(const Point& p) const {
      return (sameSide(p,vertices_[0], vertices_[1],vertices_[2])
          and sameSide(p,vertices_[1], vertices_[0],vertices_[2])
          and sameSide(p,vertices_[2], vertices_[0],vertices_[1]))
        ? true : false;
    }

    std::string toString() const {
      using std::to_string;
      return to_string(vertices_[0][0]) + "  " + to_string(vertices_[0][1])
        + "\n" + to_string(vertices_[1][0]) + "  " + to_string(vertices_[1][1])
        + "\n" + to_string(vertices_[2][0]) + "  " + to_string(vertices_[2][1])
        + "\n" + to_string(vertices_[0][0]) + "  " + to_string(vertices_[0][1])
        + "\n";
    }


  private:
    const int index_;
    const std::array<Point, 3> vertices_;
    std::array<int, 3> neighbors_; // i'th neighbor is across from i'th vertex
};


#endif // INTERP_TRIANGLE_H
