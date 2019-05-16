#include "random_points_on_mesh.h"
#include "igl/cumsum.h"
#include "igl/doublearea.h"
#include <algorithm>
#include <random>
#include <vector>

void random_points_on_mesh(
  const int n,
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & X)
{
  // REPLACE WITH YOUR CODE:
  X.resize(n,3);
  Eigen::VectorXd A, sum_vec;
  igl::doublearea(V, F, A);
  igl::cumsum(A, 1, sum_vec);
  std::vector<double> sum(sum_vec.data(), sum_vec.data() + sum_vec.size());

  std::uniform_real_distribution<double> uniform_sum(0, sum.back());
  std::uniform_real_distribution<double> uniform_unit(0, 1);
  std::default_random_engine re;

  double r, a, b;
  int fi;
  Eigen::Vector3i f;
  Eigen::Vector3d v1, v2, v3, x;
  for(int i = 0; i < X.rows(); i++) {
    r = uniform_sum(re);
    fi = (std::upper_bound(sum.begin(), sum.end(), r) - sum.begin()) - 1;
    f = F.row(fi);
    v1 = V.row(f(0));
    v2 = V.row(f(1));
    v3 = V.row(f(2));

    a = uniform_unit(re);
    b = uniform_unit(re);
    if (a + b > 1) {
      a = 1 - a;
      b = 1 - b;
    }

    x = v1 + a * (v2 - v1) + b * (v3 - v1);
    X.row(i) = x;

    
  }
}
