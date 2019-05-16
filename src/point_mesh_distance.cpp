#include "point_mesh_distance.h"
#include "point_triangle_distance.h"
#include "igl/per_face_normals.h"
#include <iostream>

void point_mesh_distance(
  const Eigen::MatrixXd & X,
  const Eigen::MatrixXd & VY,
  const Eigen::MatrixXi & FY,
  Eigen::VectorXd & D,
  Eigen::MatrixXd & P,
  Eigen::MatrixXd & N)
{
  Eigen::MatrixXd NY;
  igl::per_face_normals(VY, FY, NY);
  P.resizeLike(X);
  N = Eigen::MatrixXd::Zero(X.rows(),X.cols());
  Eigen::RowVector3d x, p, min_p, v1, v2, v3;
  Eigen::RowVector3i f;
  double d, min_d;
  int min_f;
  for(int i = 0; i < X.rows(); i++) {
    x = X.row(i);
    min_d = 100000;
    for (int j = 0; j < FY.rows(); j++) {
      f = FY.row(j);
      point_triangle_distance(x, VY.row(f(0)), VY.row(f(1)), VY.row(f(2)), d, p);
      if (d < min_d) {
        min_p = p;
        min_d = d;
        min_f = j;
      }
    }
    P.row(i) = min_p;
    N.row(i) = NY.row(min_f);
  }
  D = (X-P).rowwise().norm();
}
