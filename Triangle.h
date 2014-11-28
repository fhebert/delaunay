
#ifndef DELAUNAY_TRIANGLE_H
#define DELAUNAY_TRIANGLE_H


#include "Connectivity.h"
#include "Point.h"
#include "VectorOps.h"
#include "Utils.h"

#include <string>
#include <vector>


class Triangle {

  public:
    Triangle(const Connectivity& conn, const std::vector<Point>& points)
      : conn_(conn), v0_(points[conn.vertex(0)]), v1_(points[conn.vertex(1)]), v2_(points[conn.vertex(2)]) {}


    int vertex(const int i) const {
      return conn_.vertex(i);
    }

    int neighbor(const int i) const {
      return conn_.neighbor(i);
    }

    bool isPointInside(const Point& p) const {
      return (sameSide(p,v0_, v1_,v2_) and sameSide(p,v1_, v0_,v2_) and sameSide(p,v2_, v0_,v1_));
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
    const Connectivity& conn_;
    const Point& v0_;
    const Point& v1_;
    const Point& v2_;
};


#endif // DELAUNAY_TRIANGLE_H
