#include "dynamical_system.hpp"
#include "essential.hpp"

dynamical_system::dynamical_system(nlohmann::json json) {
  xdim = json["x0"].size();
  bialt_dim = 0;
  for (int i = 0; i < xdim; i++)
    bialt_dim += i;
  size_dphidx = xdim * xdim;
  size_dphidlambda = xdim;
  size_dphidxdx = xdim * xdim * xdim;
  size_dphidxdlambda = xdim * xdim;

  period = json["period"];
  inc_param = json["inc_param"];
  var_param = json["var_param"];
  delta_inc = json["delta_inc"];
  inc_iter = json["inc_iter"];
  max_iter = json["max_iter"];
  numerical_diff = json["numerical_diff"];
  dif_strip = json["dif_strip"];
  eps = json["eps"];
  explode = json["explode"];
  mode = json["mode"];

  out_path = json["out_path"];
  json_out_path = json["json_out_path"];

  /* These json array should be casted to the STL container type*/
  std::vector<double> fixed_arr = json["x0"];
  Eigen::Map<Eigen::VectorXd> x0(fixed_arr.data(), fixed_arr.size());
  this->x0 = x0;

  std::vector<double> params_arr = json["params"];
  Eigen::Map<Eigen::VectorXd> params(params_arr.data(), params_arr.size());
  this->p = params;

  f = Eigen::VectorXd::Zero(xdim);
  dfdx = Eigen::MatrixXd::Zero(xdim, xdim);
  dfdlambda = Eigen::VectorXd::Zero(xdim);
  dfdxdlambda = Eigen::MatrixXd::Zero(xdim, xdim);
  dfdxdx =
      std::vector<Eigen::MatrixXd>(xdim, Eigen::MatrixXd::Zero(xdim, xdim));

  dphidx = Eigen::MatrixXd::Zero(xdim, xdim);
  dphidlambda = Eigen::VectorXd::Zero(xdim);
  dphidxdlambda = Eigen::MatrixXd::Zero(xdim, xdim);
  dphidxdx =
      std::vector<Eigen::MatrixXd>(xdim, Eigen::MatrixXd::Zero(xdim, xdim));

  xk = std::vector<Eigen::VectorXd>(period + 1, Eigen::VectorXd::Zero(xdim));
}

void dynamical_system::store_states(const Eigen::VectorXd &v) {
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(xdim, xdim);
  int ve_dim = xdim + size_dphidx + size_dphidlambda + size_dphidxdx +
               size_dphidxdlambda;
  Eigen::VectorXd O = Eigen::VectorXd::Zero(size_dphidlambda + size_dphidxdx +
                                            size_dphidxdlambda);
  Eigen::VectorXd state(ve_dim);

  xk[0] = v(Eigen::seqN(0, xdim));
  p(var_param) = v(xdim);

  I.resize(I.cols() * I.rows(), 1);
  state << xk[0], I, O;

  for (int i = 0; i < period; i++) {
    function(i, state);
    xk[i + 1] = state(Eigen::seqN(0, xdim));
  }

  unsigned int counter = xdim;
  Eigen::MatrixXd dphidx_buf;
  dphidx_buf = state(Eigen::seqN(counter, size_dphidx));
  dphidx_buf.resize(xdim, xdim);
  dphidx = dphidx_buf;
  counter += size_dphidx;

  dphidlambda = state(Eigen::seqN(counter, size_dphidlambda));
  counter += size_dphidlambda;

  Eigen::MatrixXd dphidxdx_buf;
  for (int i = 0; i < xdim; i++) {
    dphidxdx_buf = state(Eigen::seqN(counter + size_dphidx * i, size_dphidx));
    dphidxdx_buf.resize(xdim, xdim);
    dphidxdx[i] = dphidxdx_buf;
    dphidxdx_buf.resize(size_dphidx, 1);
  }
  counter += size_dphidxdx;

  Eigen::MatrixXd dphidxdlambda_buf;
  dphidxdlambda_buf = state(Eigen::seqN(counter, size_dphidxdlambda));
  dphidxdlambda_buf.resize(xdim, xdim);
  dphidxdlambda = dphidxdlambda_buf;
  counter += size_dphidxdlambda;

  dTldx = dphidx;
  dTldlambda = dphidlambda;
  dTldxdx = dphidxdx;
  dTldxdlambda = dphidxdlambda;

  eigvals = Eigen::EigenSolver<Eigen::MatrixXd>(dTldx).eigenvalues();

  switch (mode) {
  case 1:
    chara_poly = dTldx - Eigen::MatrixXd::Identity(xdim, xdim);
    break;
  case 2:
    chara_poly = dTldx + Eigen::MatrixXd::Identity(xdim, xdim);
    break;
  case 3:
    chara_poly = bialt_prod_square(dTldx, xdim, bialt_dim) -
                 Eigen::MatrixXd::Identity(bialt_dim, bialt_dim);
    break;
  }
}

void dynamical_system::store_states_numeric(const Eigen::VectorXd &v) {
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(xdim, xdim);
  int ve_dim = xdim + size_dphidx + size_dphidlambda;
  Eigen::VectorXd O = Eigen::VectorXd::Zero(size_dphidlambda);
  Eigen::VectorXd state(ve_dim);
  Eigen::VectorXd state_init(ve_dim);
  Eigen::VectorXd h = Eigen::VectorXd::Zero(xdim);
  Eigen::MatrixXd Ah;

  xk[0] = v(Eigen::seqN(0, xdim));
  p(var_param) = v(xdim);

  I.resize(I.cols() * I.rows(), 1);
  state << xk[0], I, O;
  state_init = state;

  for (int i = 0; i < period; i++) {
    function(i, state);
    xk[i + 1] = state(Eigen::seqN(0, xdim));
  }

  unsigned int counter = xdim;
  Eigen::MatrixXd dphidx_buf;
  dphidx_buf = state(Eigen::seqN(counter, size_dphidx));
  dphidx_buf.resize(xdim, xdim);
  dphidx = dphidx_buf;
  counter += size_dphidx;

  dphidlambda = state(Eigen::seqN(counter, size_dphidlambda));
  counter += size_dphidlambda;

  // dTldxdx (numerical differenciation)
  for (int i = 0; i < xdim; i++) {
    state = state_init;
    state(i) += dif_strip;
    // make Ah dphidx(x_k + h, lambda)
    for (int j = 0; j < period; j++) {
      function(j, state);
    }
    Ah = state(Eigen::seqN(xdim, size_dphidx));
    Ah.resize(xdim, xdim);
    dphidxdx[i] = (Ah - dphidx) / dif_strip;
    Ah.resize(size_dphidx, 1);
  }

  // dTldxdlambda (numerical differenciation)
  static double lambda_temp;
  lambda_temp = p(var_param);
  p(var_param) += dif_strip;
  // make Ah dphidx(x_k, lambda + h)
  state = state_init;
  for (int i = 0; i < period; i++) {
    function(i, state);
  }
  Ah = state(Eigen::seqN(xdim, size_dphidx));
  Ah.resize(xdim, xdim);
  dphidxdlambda = (Ah - dphidx) / dif_strip;
  p(var_param) = lambda_temp;

  dTldx = dphidx;
  dTldlambda = dphidlambda;
  dTldxdx = dphidxdx;
  dTldxdlambda = dphidxdlambda;

  eigvals = Eigen::EigenSolver<Eigen::MatrixXd>(dTldx).eigenvalues();

  switch (mode) {
  case 1:
    chara_poly = dTldx - Eigen::MatrixXd::Identity(xdim, xdim);
    break;
  case 2:
    chara_poly = dTldx + Eigen::MatrixXd::Identity(xdim, xdim);
    break;
  case 3:
    chara_poly = bialt_prod_square(dTldx, xdim, bialt_dim) -
                 Eigen::MatrixXd::Identity(bialt_dim, bialt_dim);
    break;
  }
}

Eigen::VectorXd dynamical_system::newton_F_G() {
  Eigen::VectorXd F(xdim + 1);

  F(Eigen::seqN(0, xdim)) = xk[period] - xk[0];
  F(xdim) = chara_poly.determinant();

  return F;
}

Eigen::MatrixXd dynamical_system::newton_J_G() {
  Eigen::MatrixXd J(xdim + 1, xdim + 1);
  Eigen::MatrixXd dchidx(1, xdim);
  double dchidlambda;
  for (int i = 0; i < xdim; i++) {
    dchidx(0, i) = det_derivative(chara_poly, dTldxdx[i], xdim);
  }
  dchidlambda = det_derivative(chara_poly, dTldxdlambda, xdim);

  J(Eigen::seqN(0, xdim), Eigen::seqN(0, xdim)) =
      dTldx - Eigen::MatrixXd::Identity(xdim, xdim);
  J(Eigen::seqN(0, xdim), xdim) = dTldlambda;
  J(xdim, Eigen::seqN(0, xdim)) = dchidx;
  J(xdim, xdim) = dchidlambda;

  return J;
}

Eigen::VectorXd dynamical_system::newton_F_PD() {
  Eigen::VectorXd F(xdim + 1);

  F(Eigen::seqN(0, xdim)) = xk[period] - xk[0];
  F(xdim) = chara_poly.determinant();

  return F;
}

Eigen::MatrixXd dynamical_system::newton_J_PD() {
  Eigen::MatrixXd J(xdim + 1, xdim + 1);
  Eigen::MatrixXd dchidx(1, xdim);
  double dchidlambda;
  for (int i = 0; i < xdim; i++) {
    dchidx(0, i) = det_derivative(chara_poly, dTldxdx[i], xdim);
  }
  dchidlambda = det_derivative(chara_poly, dTldxdlambda, xdim);

  J(Eigen::seqN(0, xdim), Eigen::seqN(0, xdim)) =
      dTldx - Eigen::MatrixXd::Identity(xdim, xdim);
  J(Eigen::seqN(0, xdim), xdim) = dTldlambda;
  J(xdim, Eigen::seqN(0, xdim)) = dchidx;
  J(xdim, xdim) = dchidlambda;

  return J;
}

Eigen::VectorXd dynamical_system::newton_F_NS() {
  Eigen::VectorXd F(xdim + 1);

  F(Eigen::seqN(0, xdim)) = xk[period] - xk[0];
  F(xdim) = chara_poly.determinant();

  return F;
}

Eigen::MatrixXd dynamical_system::newton_J_NS() {
  Eigen::MatrixXd J(xdim + 1, xdim + 1);
  Eigen::MatrixXd dchidx(1, xdim);
  double dchidlambda;
  Eigen::MatrixXd temp(bialt_dim, bialt_dim);
  for (int i = 0; i < xdim; i++) {
    temp = bialt_prod_square_derivative(dTldx, dTldxdx[i], xdim, bialt_dim);
    dchidx(0, i) = det_derivative(chara_poly, temp, bialt_dim);
  }
  temp = bialt_prod_square_derivative(dTldx, dTldxdlambda, xdim, bialt_dim);
  dchidlambda = det_derivative(chara_poly, temp, bialt_dim);

  J(Eigen::seqN(0, xdim), Eigen::seqN(0, xdim)) =
      dTldx - Eigen::MatrixXd::Identity(xdim, xdim);
  J(Eigen::seqN(0, xdim), xdim) = dTldlambda;
  J(xdim, Eigen::seqN(0, xdim)) = dchidx;
  J(xdim, xdim) = dchidlambda;

  return J;
}

void dynamical_system::store_states_fix(const Eigen::VectorXd &v) {
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(xdim, xdim);
  Eigen::VectorXd state(xdim + size_dphidx);

  xk[0] = v;

  I.resize(I.cols() * I.rows(), 1);
  state << xk[0], I;

  for (int i = 0; i < period; i++) {
    function(i, state);
    xk[i + 1] = state(Eigen::seqN(0, xdim));
  }
  Eigen::MatrixXd dphidx_buf;
  dphidx_buf = state(Eigen::seqN(xdim, size_dphidx));
  dphidx_buf.resize(xdim, xdim);

  dphidx = dphidx_buf;

  dTldx = dphidx;

  eigvals = Eigen::EigenSolver<Eigen::MatrixXd>(dTldx).eigenvalues();
}

Eigen::VectorXd dynamical_system::newton_F_fix() { return xk[period] - xk[0]; }

Eigen::MatrixXd dynamical_system::newton_J_fix() {
  return dTldx - Eigen::MatrixXd::Identity(xdim, xdim);
}

std::tuple<Eigen::VectorXd, Eigen::MatrixXd>
dynamical_system::newton_FJ(const Eigen::VectorXd &v) {
  Eigen::VectorXd F;
  Eigen::MatrixXd J;

  if (mode != 0) {
    if (!numerical_diff) {
      store_states(v);
    } else {
      store_states_numeric(v);
    }
    switch (mode) {
    case 1:
      F = newton_F_G();
      J = newton_J_G();
      break;
    case 2:
      F = newton_F_PD();
      J = newton_J_PD();
      break;
    case 3:
      F = newton_F_NS();
      J = newton_J_NS();
      break;
    }
  } else {
    store_states_fix(v);
    F = newton_F_fix();
    J = newton_J_fix();
  }

  return std::make_tuple(F, J);
}