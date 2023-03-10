#include <algorithm>
#include <cassert>
#include <cmath>

#include "typical-bezier.hh"

using namespace Geometry;

static double angleBetween(const Point2D &u, const Point2D &v) {
  auto angle = std::acos(std::min(1.0, std::max(-1.0, u.normalized() * v.normalized())));
  return angle * (u[0] * v[1] < u[1] * v[0] ? -1 : 1);
}

static Vector2D rotate(const Vector2D &p, double alpha) {
  auto c = std::cos(alpha), s = std::sin(alpha);
  return { c * p[0] - s * p[1], s * p[0] + c * p[1] };
}

static Vector2DVector generateRotations(Vector2D u, double alpha, size_t k) {
  Vector2DVector result;
  for (size_t i = 0; i < k; ++i) {
    result.push_back(u);
    u = rotate(u, alpha);
  }
  return result;
}

static double solve(const DoubleVector &a) {
  assert(a.size() == 3); // only quadratic is supported for now
  auto D = a[1] * a[1] - 4 * a[0] * a[2];
  assert(D >= 0);
  D = std::sqrt(D);
  auto x1 = (-a[1] - D) / (2 * a[2]);
  auto x2 = (-a[1] + D) / (2 * a[2]);
  // If one is negative, return the other
  if (x1 * x2 < 0)
    return x1 < 0 ? x2 : x1;
  // Return the one closer to 0
  if (std::abs(x1) < std::abs(x2))
    return x1;
  return x2;
}

static Point2D otherEnd(const Vector2DVector &coeffs, double s) {
  Point2D p(0, 0);
  double si = 1;
  for (const auto &v : coeffs) {
    p += v * si;
    si *= s;
  }
  return p;
}

static std::pair<double, double> optimize(const Vector2DVector &coeffs, const Vector2D &target) {
  DoubleVector vs;
  Vector2D perp(-target[1], target[0]);
  perp.normalize();
  for (const auto &p : coeffs)
    vs.push_back(p * perp);
  auto s = solve(vs);
  return { s, target.norm() / otherEnd(coeffs, s).norm() };
}

static Point2DVector generateBezierControls(Point2D p, Vector2D v,
                                            double s, double alpha, size_t degree) {
  Point2DVector result;
  for (size_t i = 0; i <= degree; ++i) {
    result.push_back(p);
    p += v;
    v = rotate(v, alpha) * s;
  }
  return result;
}

Point2DVector typicalBezier(const Point2DVector &a, size_t degree) {
  auto angle = angleBetween(a[1] - a[0], a[2] - a[1]);
  auto alpha = angle / (degree - 1);
  auto u = (a[1] - a[0]).normalize();
  auto vs = generateRotations(u, alpha, degree);
  const auto &[s, b0] = optimize(vs, a[2] - a[0]);
  auto v = u * b0;
  return generateBezierControls(a[0], v, s, alpha, degree);
}
