#include "closest_rotation.h"
#include <Eigen/Geometry>

void closest_rotation(
  const Eigen::Matrix3d & M,
  Eigen::Matrix3d & R)
{
  Eigen::JacobiSVD<Eigen::Matrix3d> svd(M, Eigen::ComputeThinU | Eigen::ComputeThinV);
  Eigen::Matrix3d U, V, O;
  U = svd.matrixU();
  V = svd.matrixV();
  O << 1, 0, 0,
       0, 1, 0,
       0, 0, (U*V.transpose()).determinant();
  R = U * O * V.transpose();
}
