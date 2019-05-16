#include "point_triangle_distance.h"
#include <Eigen/Geometry>
#include <iostream>
#include <algorithm>

void _point_line_projection(
  const Eigen::RowVector3d & x,
  const Eigen::RowVector3d & a,
  const Eigen::RowVector3d & b,
  Eigen::RowVector3d & p
  )
{
  // Project point x onto the line segment between points a and b.
  auto d = (x - a).dot(b - a) / (b - a).dot(b - a) * (b - a);
  int sign = (b - a).dot(d) < 0 ? -1 : 1;
  double norm_d = sign * d.norm();
  auto unit_d = d.normalized();
  norm_d = std::min(std::max(0., norm_d), (b - a).norm());
  p = a + norm_d * unit_d;
}

void point_triangle_distance(
  const Eigen::RowVector3d & x,
  const Eigen::RowVector3d & a,
  const Eigen::RowVector3d & b,
  const Eigen::RowVector3d & c,
  double & d,
  Eigen::RowVector3d & p)
{
  Eigen::Vector3d e1, e2, v, n;
  double alpha, beta, gamma;

  // Adapted from 
  // Christer Ericson. 2004. Real-Time Collision Detection. pp. 47-48. CRC Press, Inc., Boca Raton, FL, USA.
  e1 = b - a;
  e2 = c - a;
  v = x - a;
  double d01 = e1.dot(e2);
  double d00 = e1.dot(e1);
  double d11 = e2.dot(e2);
  double d20 = v.dot(e1);
  double d21 = v.dot(e2);
  double denom = d00 * d11 - d01 * d01;
  beta = (d11 * d20 - d01 * d21) / denom;
  gamma = (d00 * d21 - d01 * d20) / denom;
  alpha = 1. - beta - gamma;

  if (gamma < 0) {
    _point_line_projection(x, a, b, p);
  }
  else if (beta < 0) {
    _point_line_projection(x, a, c, p);
  }
  else if (alpha < 0) {
    _point_line_projection(x, b, c, p);
  }
  else {
    p = alpha * a + beta * b + gamma * c;
  }
  d = (x - p).norm();
}
