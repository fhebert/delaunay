
#ifndef DELAUNAY_INTERPOLATOR_H
#define DELAUNAY_INTERPOLATOR_H

#include "DelaunayTri.h"
#include "InterpolatorHelper.h"
#include "Triangle.h"
#include "Point.h"

#include <cassert>


class Interpolator {

  public:
    Interpolator(const DelaunayTri& tri,
        const std::vector<std::vector<double>>& data)
      : fill_value_(1e234), tri_(tri), data_(data)
    {
      for (const auto dataset : data) {
        assert(dataset.size()==tri.getNumPoints() and
            "datasets incompatible with triangulation");
      }
      gradients_ = estimateGradientsGlobal(tri, data);
    }


    std::vector<double> evaluateAt(const Point& x) const
    {
      int isimplex, iedge;
      std::tie(isimplex,iedge) = tri_.findEnclosingTriangleIndex(x);

      if (isimplex == -1) {
        // outside triangulation
        return std::vector<double>(data_.size(), fill_value_);
      }

      // find coordinates of x
      const Triangle& t = tri_.triangle(isimplex);
      const std::array<double,3> b = t.barycentricCoords(x);

      std::vector<double> result(data_.size(),0);
      for (size_t k=0; k<data_.size(); ++k) {
        // prepare values and gradients
        const std::array<double,3> f = {{
          data_[k][t.vertex(0)],
            data_[k][t.vertex(1)],
            data_[k][t.vertex(2)]}};
        const std::array<std::array<double,2>,3> df = {{
          gradients_[k][t.vertex(0)],
            gradients_[k][t.vertex(1)],
            gradients_[k][t.vertex(2)]}};
        result[k] = scipy_clough_tocher_2d_single(tri_, isimplex, b, f, df);
      }
      return result;
    }


  private:
    std::vector<std::vector<std::array<double,2>>> estimateGradientsGlobal(
        const DelaunayTri& tri,
        const std::vector<std::vector<double>>& data,
        const size_t maxiter=400,
        const double tol=1e-6)
    {
      const size_t ndatasets = data.size();
      std::vector<std::vector<std::array<double,2>>> grad(
          ndatasets, std::vector<std::array<double,2>>(tri.getNumPoints(), {{0,0}}));

      for (size_t k=0; k<ndatasets; ++k) {
        const int ret = scipy_estimate_gradients_2d_global(
            tri,
            data[k],
            maxiter,
            tol,
            grad[k]);

        assert(ret > 0 and "Gradient estimation did not converge, the results may be inaccurate");
      }

      return grad;
    }


  private:
    const double fill_value_;
    const DelaunayTri& tri_;
    std::vector<std::vector<double>> data_;
    std::vector<std::vector<std::array<double,2>>> gradients_;
};


#endif // DELAUNAY_INTERPOLATOR_H
