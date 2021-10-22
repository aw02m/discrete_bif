#include "dynamical_system.hpp"

dynamical_system::dynamical_system(nlohmann::json json) {
  xdim = json["x0"].size();
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

  T = Eigen::VectorXd::Zero(xdim);
  dTdx = Eigen::MatrixXd::Zero(xdim, xdim);
  dTdlambda = Eigen::VectorXd::Zero(xdim);
  dTdxdlambda = Eigen::MatrixXd::Zero(xdim, xdim);
  dTdxdx =
      std::vector<Eigen::MatrixXd>(xdim, Eigen::MatrixXd::Zero(xdim, xdim));

  xk = std::vector<Eigen::VectorXd>(period + 1, Eigen::VectorXd::Zero(xdim));
  dTdx_arr =
      std::vector<Eigen::MatrixXd>(period, Eigen::MatrixXd::Zero(xdim, xdim));
  dTdlambda_arr =
      std::vector<Eigen::VectorXd>(period, Eigen::VectorXd::Zero(xdim));
  dTdxdx_arr = std::vector<std::vector<Eigen::MatrixXd>>(
      period,
      std::vector<Eigen::MatrixXd>(xdim, Eigen::MatrixXd::Zero(xdim, xdim)));
  dTdxdlambda_arr =
      std::vector<Eigen::MatrixXd>(period, Eigen::MatrixXd::Zero(xdim, xdim));
  dTldxdx =
      std::vector<Eigen::MatrixXd>(xdim, Eigen::MatrixXd::Zero(xdim, xdim));

  frwd_prod_arr = std::vector<Eigen::MatrixXd>(
      period, Eigen::MatrixXd::Identity(xdim, xdim));
  bkwd_prod_arr = std::vector<Eigen::MatrixXd>(
      period, Eigen::MatrixXd::Identity(xdim, xdim));
  dTkdlambda_arr =
      std::vector<Eigen::VectorXd>(period, Eigen::VectorXd::Zero(xdim));
}

void dynamical_system::store_states(const Eigen::VectorXd &v) {
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(xdim, xdim);

  xk[0] = v(Eigen::seqN(0, xdim));
  p(var_param) = v(xdim);
  theta = v(xdim + 1);

  for (int i = 0; i < period; i++) {
    function(xk[i]);
    xk[i + 1] = T;
    dTdx_arr[i] = dTdx;
    dTdlambda_arr[i] = dTdlambda;
    dTdxdx_arr[i] = dTdxdx;
    dTdxdlambda_arr[i] = dTdxdlambda;
  }

  // dTldx
  dTldx = Eigen::MatrixXd::Identity(xdim, xdim);
  for (int i = period - 1; i >= 0; i--) {
    dTldx *= dTdx_arr[i];
  }
  // dTldlambda
  dTldlambda = dTdlambda_arr[0];
  for (int i = 1; i < period; i++) {
    dTldlambda = dTdx_arr[i] * dTldlambda + dTdlambda_arr[i];
    dTkdlambda_arr[i] = dTldlambda; // buffer for second derivatives
  }
  // dTldxdx
  dTldxdx =
      std::vector<Eigen::MatrixXd>(xdim, Eigen::MatrixXd::Zero(xdim, xdim));
  frwd_prod_arr = std::vector<Eigen::MatrixXd>(
      period, Eigen::MatrixXd::Identity(xdim, xdim));
  bkwd_prod_arr = std::vector<Eigen::MatrixXd>(
      period, Eigen::MatrixXd::Identity(xdim, xdim));
  for (int i = 0; i < period; i++) {
    for (int j = period - 1; j >= 0; j--) {
      if (j > i) {
        frwd_prod_arr[i] *= dTdx_arr[j];
      } else if (j < i) {
        bkwd_prod_arr[i] *= dTdx_arr[j];
      }
    }
  }
  for (int i = 0; i < xdim; i++) {
    for (int j = 0; j < period; j++) {
      dTldxdx[i] += frwd_prod_arr[j] * dTdxdx_arr[j][i] * bkwd_prod_arr[j] *
                    bkwd_prod_arr[j];
    }
  }
  // dTldxdlambda
  dTldxdlambda = dTdxdlambda_arr[0];
  Eigen::MatrixXd temp = Eigen::MatrixXd::Zero(xdim, xdim);
  for (int i = 1; i < period; i++) {
    for (int j = 0; j < xdim; j++) {
      temp.col(j) = dTdxdx_arr[i][j] * bkwd_prod_arr[i] * dTkdlambda_arr[i];
    }
    dTldxdlambda = temp + dTdx_arr[i] * dTldxdlambda +
                   dTdxdlambda_arr[i] * bkwd_prod_arr[i];
  }

  eigvals = Eigen::EigenSolver<Eigen::MatrixXd>(dTldx).eigenvalues();
  mu = Eigen::dcomplex(std::cos(theta), std::sin(theta));

  chara_poly = dTldx - mu * I;
}

void dynamical_system::store_states_numeric(const Eigen::VectorXd &v) {
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(xdim, xdim);
  Eigen::VectorXd h = Eigen::VectorXd::Zero(xdim);
  Eigen::MatrixXd Ah(xdim, xdim);

  xk[0] = v(Eigen::seqN(0, xdim));
  p(var_param) = v(xdim);
  theta = v(xdim + 1);

  for (int i = 0; i < period; i++) {
    function(xk[i]);
    xk[i + 1] = T;
    dTdx_arr[i] = dTdx;
    dTdlambda_arr[i] = dTdlambda;
  }

  // dTldx
  dTldx = Eigen::MatrixXd::Identity(xdim, xdim);
  for (int i = period - 1; i >= 0; i--) {
    dTldx *= dTdx_arr[i];
  }
  // dTldlambda
  dTldlambda = dTdlambda_arr[0];
  for (int i = 1; i < period; i++) {
    dTldlambda = dTdx_arr[i] * dTldlambda + dTdlambda_arr[i];
    dTkdlambda_arr[i] = dTldlambda; // buffer for second derivatives
  }
  // dTldxdx (numerical differenciation)
  for (int i = 0; i < xdim; i++) {
    Ah = Eigen::MatrixXd::Identity(xdim, xdim);
    h(i) = dif_strip;
    // make Ah dTldx(x_k + h)
    for (int j = period - 1; j >= 0; j--) {
      function(xk[j] + h);
      Ah *= dTdx;
    }
    dTldxdx[i] = (Ah - dTldx) / dif_strip;
    h(i) = 0;
  }
  // dTldxdlambda (numerical differenciation)
  static double lambda_temp;
  lambda_temp = p(var_param);
  p(var_param) += dif_strip;
  Ah = Eigen::MatrixXd::Identity(xdim, xdim);
  // make Ah dTldx(x_k, lambda + h)
  for (int i = period - 1; i >= 0; i--) {
    function(xk[i]);
    Ah *= dTdx;
  }
  dTldxdlambda = (Ah - dTldx) / dif_strip;
  p(var_param) = lambda_temp;

  eigvals = Eigen::EigenSolver<Eigen::MatrixXd>(dTldx).eigenvalues();
  mu = Eigen::dcomplex(std::cos(theta), std::sin(theta));

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

  xk[0] = v(Eigen::seqN(0, xdim));

  for (int i = 0; i < period; i++) {
    function(xk[i]);
    xk[i + 1] = T;
    dTdx_arr[i] = dTdx;
  }

  // dTldx
  dTldx = Eigen::MatrixXd::Identity(xdim, xdim);
  for (int i = period - 1; i >= 0; i--) {
    dTldx *= dTdx_arr[i];
  }

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
