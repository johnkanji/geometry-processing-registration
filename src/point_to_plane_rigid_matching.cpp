#include "point_to_plane_rigid_matching.h"
#include "closest_rotation.h"
#include "igl/cat.h"

void point_to_plane_rigid_matching(
  const Eigen::MatrixXd & X,
  const Eigen::MatrixXd & P,
  const Eigen::MatrixXd & N,
  Eigen::Matrix3d & R,
  Eigen::RowVector3d & t)
{
  Eigen::MatrixXd A, N_diag, N1, N2, N3, NA, tmp;
  Eigen::Matrix3d M;
  Eigen::VectorXd u, zero, ones, X_vec, P_vec, b, Nb;
  
  zero = Eigen::VectorXd::Zero(X.rows());
  ones = Eigen::VectorXd::Ones(X.rows());
  A.resize(3 * X.rows(), 6);
  A << zero, X.col(2), -X.col(1), ones, zero, zero,
       -X.col(2), zero, X.col(0), zero, ones, zero,
       X.col(1), -X.col(0), zero, zero, zero, ones;
  N1 = N.col(0).asDiagonal();
  N2 = N.col(1).asDiagonal();
  N3 = N.col(2).asDiagonal();
  igl::cat(2, N1, N2, tmp);
  igl::cat(2, tmp, N3, N_diag);
  NA = N_diag * A;


  X_vec.resize(X.rows() * 3);
  X_vec << X.col(0), X.col(1), X.col(2);
  P_vec.resize(P.rows() * 3);
  P_vec << P.col(0), P.col(1), P.col(2);
  b = X_vec - P_vec;
  Nb = N_diag * b;

  u = (NA.transpose()*NA).inverse() * (-NA.transpose() * Nb);

  M << 1, u(2), -u(1),
       -u(2), 1, u(0),
       u(1), -u(0), 1;
  closest_rotation(M, R);
  t << u(3), u(4), u(5);
}
