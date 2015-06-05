
#ifndef DELAUNAY_INTERPOLATOR_HELPER_H
#define DELAUNAY_INTERPOLATOR_HELPER_H

#include "DelaunayTri.h"
#include "Point.h"

#include <cassert>


int scipy_estimate_gradients_2d_global(
    const DelaunayTri& d,
    const std::vector<double> data,
    const size_t maxiter,
    const double tol,
    std::vector<std::array<double,2>>& y)
{
  /*
  Estimate gradients of a function at the vertices of a 2d triangulation.
  Parameters
  ----------
  info : input
    Triangulation in 2D
  data : input
    Function values at the vertices
  maxiter : input
    Maximum number of Gauss-Seidel iterations
  tol : input
    Absolute / relative stop tolerance
  y : output, shape (npoints, 2)
    Derivatives [F_x, F_y] at the vertices
  Returns
  -------
  num_iterations
    Number of iterations if converged, 0 if maxiter reached
    without convergence
  Notes
  -----
  This routine uses a re-implementation of the global approximate
  curvature minimization algorithm described in [Nielson83] and [Renka84].
  References
  ----------
  .. [Nielson83] G. Nielson,
     ''A method for interpolating scattered data based upon a minimum norm
     network''.
     Math. Comp., 40, 253 (1983).
  .. [Renka84] R. J. Renka and A. K. Cline.
     ''A Triangle-based C1 interpolation method.'',
     Rocky Mountain J. Math., 14, 223 (1984).
  */

  std::array<double,4> Q;
  std::array<double,2> r,s;

  for (auto ypoint : y) {
    ypoint = {{0,0}};
  }

  for (size_t iiter=0; iiter<maxiter; ++iiter) {
    double err = 0;
    for (size_t ipoint=0; ipoint<d.getNumPoints(); ++ipoint) {
      Q = {{0, 0, 0, 0}};
      s = {{0, 0}};

      // walk over neighbours of given point
      const auto neighbors = d.getConnectedPoints(ipoint);
      for (const size_t ipoint2 : neighbors) {

        // edge
        const double ex = d.point(ipoint2)[0] - d.point(ipoint)[0];
        const double ey = d.point(ipoint2)[1] - d.point(ipoint)[1];
        const double L = sqrt(ex*ex + ey*ey);
        const double L3 = L*L*L;

        // data at vertices
        const double f1 = data[ipoint];
        const double f2 = data[ipoint2];

        // scaled gradient projections on the edge
        const double df2 = -ex*y[ipoint2][0] - ey*y[ipoint2][1];

        // edge sum
        Q[0] += 4*ex*ex / L3;
        Q[1] += 4*ex*ey / L3;
        Q[3] += 4*ey*ey / L3;

        s[0] += (6*(f1 - f2) - 2*df2) * ex / L3;
        s[1] += (6*(f1 - f2) - 2*df2) * ey / L3;
      }

      Q[2] = Q[1];

      // solve

      const double det = Q[0]*Q[3] - Q[1]*Q[2];
      r[0] = ( Q[3]*s[0] - Q[1]*s[1])/det;
      r[1] = (-Q[2]*s[0] + Q[0]*s[1])/det;

      double change = fmax(fabs(y[ipoint][0] + r[0]), fabs(y[ipoint][1] + r[1]));

      y[ipoint][0] = -r[0];
      y[ipoint][1] = -r[1];

      // relative/absolute error
      change /= fmax(1.0, fmax(fabs(r[0]), fabs(r[1])));
      err = fmax(err, change);
    }

    if (err < tol) {
      return iiter + 1;
    }
  }
  // Didn't converge before maxiter
  return 0;
}


double scipy_clough_tocher_2d_single(
    const DelaunayTri& d,
    const int isimplex,
    const std::array<double,3>& b,
    const std::array<double,3>& f,
    const std::array<std::array<double,2>,3>& df)
{
  /*
  Evaluate Clough-Tocher interpolant on a 2D triangle.
  Parameters
  ----------
  d :
    Delaunay info
  isimplex : int
    Triangle to evaluate on
  b : shape (3,)
    Barycentric coordinates of the point on the triangle
  f : shape (3,)
    Function values at vertices
  df : shape (3, 2)
    Gradient values at vertices
  Returns
  -------
  w :
    Value of the interpolant at the given point
  References
  ----------
  .. [CT] See, for example,
     P. Alfeld,
     ''A trivariate Clough-Tocher scheme for tetrahedral data''.
     Computer Aided Geometric Design, 1, 169 (1984);
     G. Farin,
     ''Triangular Bernstein-Bezier patches''.
     Computer Aided Geometric Design, 3, 83 (1986).
  */

  // XXX: optimize + refactor this!

  const double e12x = (+ d.point(d.triangle(isimplex).vertex(1))[0]
      - d.point(d.triangle(isimplex).vertex(0))[0]);
  const double e12y = (+ d.point(d.triangle(isimplex).vertex(1))[1]
      - d.point(d.triangle(isimplex).vertex(0))[1]);

  const double e23x = (+ d.point(d.triangle(isimplex).vertex(2))[0]
      - d.point(d.triangle(isimplex).vertex(1))[0]);
  const double e23y = (+ d.point(d.triangle(isimplex).vertex(2))[1]
      - d.point(d.triangle(isimplex).vertex(1))[1]);

  const double e31x = (+ d.point(d.triangle(isimplex).vertex(0))[0]
      - d.point(d.triangle(isimplex).vertex(2))[0]);
  const double e31y = (+ d.point(d.triangle(isimplex).vertex(0))[1]
      - d.point(d.triangle(isimplex).vertex(2))[1]);

  const double f1 = f[0];
  const double f2 = f[1];
  const double f3 = f[2];

  const double df12 = +(df[0][0]*e12x + df[0][1]*e12y);
  const double df21 = -(df[1][0]*e12x + df[1][1]*e12y);
  const double df23 = +(df[1][0]*e23x + df[1][1]*e23y);
  const double df32 = -(df[2][0]*e23x + df[2][1]*e23y);
  const double df31 = +(df[2][0]*e31x + df[2][1]*e31y);
  const double df13 = -(df[0][0]*e31x + df[0][1]*e31y);

  const double c3000 = f1;
  const double c2100 = (df12 + 3*c3000)/3;
  const double c2010 = (df13 + 3*c3000)/3;
  const double c0300 = f2;
  const double c1200 = (df21 + 3*c0300)/3;
  const double c0210 = (df23 + 3*c0300)/3;
  const double c0030 = f3;
  const double c1020 = (df31 + 3*c0030)/3;
  const double c0120 = (df32 + 3*c0030)/3;

  const double c2001 = (c2100 + c2010 + c3000)/3;
  const double c0201 = (c1200 + c0300 + c0210)/3;
  const double c0021 = (c1020 + c0120 + c0030)/3;

  //
  // Now, we need to impose the condition that the gradient of the spline
  // to some direction `w` is a linear function along the edge.
  //
  // As long as two neighbouring triangles agree on the choice of the
  // direction `w`, this ensures global C1 differentiability.
  // Otherwise, the choice of the direction is arbitrary (except that
  // it should not point along the edge, of course).
  //
  // In [CT]_, it is suggested to pick `w` as the normal of the edge.
  // This choice is given by the formulas
  //
  //    w_12 = E_24 + g1 * E_23
  //    w_23 = E_34 + g2 * E_31
  //    w_31 = E_14 + g3 * E_12
  //
  //    g1 = -(e24x*e23x + e24y*e23y) / (e23x**2 + e23y**2)
  //    g2 = -(e34x*e31x + e34y*e31y) / (e31x**2 + e31y**2)
  //    g3 = -(e14x*e12x + e14y*e12y) / (e12x**2 + e12y**2)
  //
  // However, this choice gives an interpolant that is *not*
  // invariant under affine transforms. This has some bad
  // consequences: for a very narrow triangle, the spline can
  // develops huge oscillations. For instance, with the input data
  //
  //     [(0, 0), (0, 1), (eps, eps)],   eps = 0.01
  //     F  = [0, 0, 1]
  //     dF = [(0,0), (0,0), (0,0)]
  //
  // one observes that as eps -> 0, the absolute maximum value of the
  // interpolant approaches infinity.
  //
  // So below, we aim to pick affine invariant `g1`, `g2`, `g3`.
  // We choose
  //
  //     w = V_4' - V_4
  //
  // where V_4 is the centroid of the current triangle, and V_4' the
  // centroid of the neighbour. Since this quantity transforms similarly
  // as the gradient under affine transforms, the resulting interpolant
  // is affine-invariant. Moreover, two neighbouring triangles clearly
  // always agree on the choice of `w` (sign is unimportant), and so
  // this choice also makes the interpolant C1.
  //
  // The drawback here is a performance penalty, since we need to
  // peek into neighbouring triangles.
  //

  std::array<double,3> c;
  std::array<double,2> y;
  double g1 = 0;
  double g2 = 0;
  double g3 = 0;
  for (size_t k=0; k<3; ++k) {
    const int itri = d.triangle(isimplex).neighborAcrossGlobalPoint(
        d.triangle(isimplex).vertex(k));

    if (itri == -1) {
      // No neighbour.
      // Compute derivative to the centroid direction (e_12 + e_13)/2.
      if (k == 0)
        g1 = -2./3;
      else if (k == 1)
        g2 = -2./3;
      else if (k == 2)
        g3 = -2./3;
      continue;
    }

    // Centroid of the neighbour, in our local barycentric coordinates

    y[0] = (+ d.point(d.triangle(itri).vertex(0))[0]
        + d.point(d.triangle(itri).vertex(1))[0]
        + d.point(d.triangle(itri).vertex(2))[0]) / 3;

    y[1] = (+ d.point(d.triangle(itri).vertex(0))[1]
        + d.point(d.triangle(itri).vertex(1))[1]
        + d.point(d.triangle(itri).vertex(2))[1]) / 3;

    c = d.triangle(isimplex).barycentricCoords({{y[0],y[1]}});

    // Rewrite V_4'-V_4 = const*[(V_4-V_2) + g_i*(V_3 - V_2)]

    // Now, observe that the results can be written *in terms of
    // barycentric coordinates*. Barycentric coordinates stay
    // invariant under affine transformations, so we can directly
    // conclude that the choice below is affine-invariant.

    if (k == 0)
      g1 = (2*c[2] + c[1] - 1) / (2 - 3*c[2] - 3*c[1]);
    else if (k == 1)
      g2 = (2*c[0] + c[2] - 1) / (2 - 3*c[0] - 3*c[2]);
    else if (k == 2)
      g3 = (2*c[1] + c[0] - 1) / (2 - 3*c[1] - 3*c[0]);
  }

  const double c0111 = (g1*(-c0300 + 3*c0210 - 3*c0120 + c0030)
      + (-c0300 + 2*c0210 - c0120 + c0021 + c0201))/2;
  const double c1011 = (g2*(-c0030 + 3*c1020 - 3*c2010 + c3000)
      + (-c0030 + 2*c1020 - c2010 + c2001 + c0021))/2;
  const double c1101 = (g3*(-c3000 + 3*c2100 - 3*c1200 + c0300)
      + (-c3000 + 2*c2100 - c1200 + c2001 + c0201))/2;

  const double c1002 = (c1101 + c1011 + c2001)/3;
  const double c0102 = (c1101 + c0111 + c0201)/3;
  const double c0012 = (c1011 + c0111 + c0021)/3;

  const double c0003 = (c1002 + c0102 + c0012)/3;

  // extended barycentric coordinates
  double minval = b[0];
  for (size_t k=0; k<3; ++k) {
    if (b[k] < minval)
      minval = b[k];
  }

  const double b1 = b[0] - minval;
  const double b2 = b[1] - minval;
  const double b3 = b[2] - minval;
  const double b4 = 3*minval;

  // evaluate the polynomial -- the stupid and ugly way to do it,
  // one of the 4 coordinates is in fact zero
  const double w = (b1*b1*b1*c3000 + 3*b1*b1*b2*c2100 + 3*b1*b1*b3*c2010 +
      3*b1*b1*b4*c2001 + 3*b1*b2*b2*c1200 +
      6*b1*b2*b4*c1101 + 3*b1*b3*b3*c1020 + 6*b1*b3*b4*c1011 +
      3*b1*b4*b4*c1002 + b2*b2*b2*c0300 + 3*b2*b2*b3*c0210 +
      3*b2*b2*b4*c0201 + 3*b2*b3*b3*c0120 + 6*b2*b3*b4*c0111 +
      3*b2*b4*b4*c0102 + b3*b3*b3*c0030 + 3*b3*b3*b4*c0021 +
      3*b3*b4*b4*c0012 + b4*b4*b4*c0003);

  return w;
}


#endif // DELAUNAY_INTERPOLATOR_HELPER_H
