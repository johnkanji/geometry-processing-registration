#include "point_to_point_rigid_matching.h"
#include "closest_rotation.h"
#include <Eigen/Geometry>
#include <iostream>

void point_to_point_rigid_matching(
  const Eigen::MatrixXd & X,
  const Eigen::MatrixXd & P,
  Eigen::Matrix3d & R,
  Eigen::RowVector3d & t)
{
  Eigen::RowVector3d P_mean = P.colwise().mean();
  Eigen::RowVector3d X_mean = X.colwise().mean();
  Eigen::Matrix3d M = ((P.rowwise() - P_mean).transpose() * (X.rowwise() - X_mean)).transpose();
  closest_rotation(M, R);
  t = P_mean.transpose() - R * X_mean.transpose();
}

