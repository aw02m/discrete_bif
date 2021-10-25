#include "dynamical_system.hpp"

dynamical_system::dynamical_system(nlohmann::json json) {
  xdim = json["x0"].size();
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
  fix_mode = json["fix_mode"];

  out_path = json["out_path"];
  json_out_path = json["json_out_path"];

  /* These json array should be casted to the STL container type*/
  std::vector<double> fixed_arr = json["x0"];
  Eigen::Map<Eigen::VectorXd> x0(fixed_arr.data(), fixed_arr.size());
  this->x0 = x0;

  std::vector<double> params_arr = json["params"];
  Eigen::Map<Eigen::VectorXd> params(params_arr.data(), params_arr.size());
  this->p = params;

  theta = json["theta"];
  mu = Eigen::dcomplex(std::cos(theta), std::sin(theta));
  // mu = Eigen::dcomplex(json["sigma"], json["omega"]);

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
  theta = v(xdim + 1);

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
  mu = Eigen::dcomplex(std::cos(theta), std::sin(theta));

  I = Eigen::MatrixXd::Identity(xdim, xdim);
  chara_poly = dTldx - mu * I;
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
  theta = v(xdim + 1);

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
  mu = Eigen::dcomplex(std::cos(theta), std::sin(theta));

  I = Eigen::MatrixXd::Identity(xdim, xdim);
  chara_poly = dTldx - mu * I;
}

Eigen::VectorXd dynamical_system::newton_F() {
  Eigen::VectorXd F(xdim + 2);
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(xdim, xdim);
  Eigen::dcomplex chi = chara_poly.determinant();

  F(Eigen::seqN(0, xdim)) = xk[period] - xk[0];
  F(xdim) = chi.real();
  F(xdim + 1) = chi.imag();

  return F;
}

Eigen::MatrixXd dynamical_system::newton_J() {
  Eigen::MatrixXd J(xdim + 2, xdim + 2);
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(xdim, xdim);
  Eigen::MatrixXcd dchidx(1, xdim);
  Eigen::dcomplex dchidlambda(0, 0);
  Eigen::dcomplex dchidtheta(0, 0);
  for (int i = 0; i < xdim; i++) {
    dchidx(0, i) = det_derivative(chara_poly, dTldxdx[i]);
  }
  dchidlambda = det_derivative(chara_poly, dTldxdlambda);
  Eigen::MatrixXcd dpolydtheta =
      Eigen::dcomplex(std::sin(theta), -std::cos(theta)) * I;
  dchidtheta = det_derivative(chara_poly, dpolydtheta);

  J(Eigen::seqN(0, xdim), Eigen::seqN(0, xdim)) = dTldx - I;
  J(Eigen::seqN(0, xdim), xdim) = dTldlambda;
  J(Eigen::seqN(0, xdim), xdim + 1) = Eigen::VectorXd::Zero(xdim);
  J(xdim, Eigen::seqN(0, xdim)) = dchidx.real();
  J(xdim, xdim) = dchidlambda.real();
  J(xdim, xdim + 1) = dchidtheta.real();
  J(xdim + 1, Eigen::seqN(0, xdim)) = dchidx.imag();
  J(xdim + 1, xdim) = dchidlambda.imag();
  J(xdim + 1, xdim + 1) = dchidtheta.imag();

  return J;
}

void dynamical_system::store_states_fix(const Eigen::VectorXd &v) {
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(xdim, xdim);
  Eigen::VectorXd state(xdim + size_dphidx);

  xk[0] = v(Eigen::seqN(0, xdim));

  I.resize(I.cols() * I.rows(), 1);
  state << xk[0], I;

  for (int i = 0; i < period; i++) {
    function(i, state);
    xk[i + 1] = state(Eigen::seqN(0, xdim));
  }

  unsigned int counter = xdim;
  Eigen::MatrixXd dphidx_buf;
  dphidx_buf = state(Eigen::seqN(counter, size_dphidx));
  dphidx_buf.resize(xdim, xdim);
  dphidx = dphidx_buf;

  dTldx = dphidx;

  // Find the argument <theta> of the characteristic constant
  // whose absolute value is closest to 1.
  eigvals = Eigen::EigenSolver<Eigen::MatrixXd>(dTldx).eigenvalues();
  unsigned int target_index = 0;
  double delta = std::abs(eigvals(0)) - 1.0;
  double delta_buf = 0;
  for (int i = 1; i < xdim; i++) {
    delta_buf = std::abs(eigvals(i)) - 1.0;
    if (delta_buf < delta) {
      target_index = i;
    }
  }
  theta = std::arg(eigvals(target_index));
}

Eigen::VectorXd dynamical_system::newton_fix_F() {
  Eigen::VectorXd F(xdim);
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(xdim, xdim);

  F(Eigen::seqN(0, xdim)) = xk[period] - xk[0];

  return F;
}

Eigen::MatrixXd dynamical_system::newton_fix_J() {
  Eigen::MatrixXd J(xdim, xdim);
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(xdim, xdim);

  J(Eigen::seqN(0, xdim), Eigen::seqN(0, xdim)) = dTldx - I;

  return J;
}

std::tuple<Eigen::VectorXd, Eigen::MatrixXd>
dynamical_system::newton_FJ(const Eigen::VectorXd &v) {
  if (!fix_mode) {
    if (!numerical_diff) {
      store_states(v);
    } else {
      store_states_numeric(v);
    }
    return std::make_tuple(newton_F(), newton_J());
  } else {
    store_states_fix(v);
    return std::make_tuple(newton_fix_F(), newton_fix_J());
  }
}

Eigen::dcomplex dynamical_system::det_derivative(const Eigen::MatrixXcd &A,
                                                 const Eigen::MatrixXcd &dA) {
  Eigen::MatrixXcd temp(xdim, xdim);
  Eigen::dcomplex ret(0, 0);

  for (int i = 0; i < xdim; i++) {
    temp = A;
    temp.col(i) = dA.col(i).cast<Eigen::dcomplex>();
    ret += temp.determinant();
  }

  return ret;
}
