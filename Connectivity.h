
#ifndef DELAUNAY_CONNECTIVITY_H
#define DELAUNAY_CONNECTIVITY_H


#include "Utils.h"

#include <array>
#include <cassert>


class Connectivity {

  public:
    Connectivity(const int a, const int b, const int c, const int na, const int nb, const int nc)
      : vertices_({{a,b,c}}), neighbors_({{na,nb,nc}}) {}


    int vertex(const int i) const {
      assert(i==0 or i==1 or i==2);
      return vertices_[i];
    }

    int neighbor(const int i) const {
      assert(i==0 or i==1 or i==2);
      return neighbors_[i];
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


  private:
    // triangle defined by 3 vertices, has 3 neighbors
    // i'th neighbor is opposite from i'th vertex
    std::array<int, 3> vertices_;
    std::array<int, 3> neighbors_;
};


#endif // DELAUNAY_CONNECTIVITY_H
