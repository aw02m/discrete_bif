#include <eigen3/Eigen/Dense>

Eigen::MatrixXd bialt_prod(const Eigen::MatrixXd &A, const Eigen::MatrixXd &B,
                           int bialt_dim) {
  Eigen::MatrixXd ret(bialt_dim, bialt_dim);

  int row = 0;
  int col = 0;
  Eigen::Matrix2d temp;
  double aug_det;
  for (int p = 1; p < bialt_dim; p++) {
    for (int q = 0; q < p; q++) {
      for (int r = 1; r < bialt_dim; r++) {
        for (int s = 0; s < r; s++) {
          temp(0, 0) = A(p, r);
          temp(0, 1) = A(p, s);
          temp(1, 0) = B(q, r);
          temp(1, 1) = B(q, s);
          aug_det = temp.determinant();
          temp(0, 0) = B(p, r);
          temp(0, 1) = B(p, s);
          temp(1, 0) = A(q, r);
          temp(1, 1) = A(q, s);
          aug_det += temp.determinant();
          aug_det /= 2;
          ret(row, col++) = aug_det;
        }
      }
      col = 0;
      row++;
    }
  }

  return ret;
}

Eigen::MatrixXd bialt_prod_square(const Eigen::MatrixXd &A, int bialt_dim) {
  Eigen::MatrixXd ret(bialt_dim, bialt_dim);

  int row = 0;
  int col = 0;
  Eigen::Matrix2d temp;
  for (int p = 1; p < bialt_dim; p++) {
    for (int q = 0; q < p; q++) {
      for (int r = 1; r < bialt_dim; r++) {
        for (int s = 0; s < r; s++) {
          temp(0, 0) = A(p, r);
          temp(0, 1) = A(p, s);
          temp(1, 0) = A(q, r);
          temp(1, 1) = A(q, s);
          ret(row, col++) = temp.determinant();
        }
      }
      col = 0;
      row++;
    }
  }

  return ret;
}

double det_derivative(const Eigen::MatrixXd &A, const Eigen::MatrixXd &dA,
                      int dim) {
  Eigen::MatrixXd temp(dim, dim);
  double ret = 0;

  for (int i = 0; i < dim; i++) {
    temp = A;
    temp.col(i) = dA.col(i);
    ret += temp.determinant();
  }

  return ret;
}