#include "ds_derivatives.hpp"

Eigen::VectorXd T(const Eigen::VectorXd &x, const dynamical_system &ds)
{
    Eigen::VectorXd ret(ds.xdim);
    double a, b, c, d;
    unsigned int m;
    a = ds.params(0);
    b = ds.params(1);
    c = ds.params(2);
    d = ds.params(3);
    m = ds.params(4);

    ret(0) = d * x(1);
    ret(1) = a * std::pow(x(0), m) * (x(0)*x(0) - b*b) + c*x(0);
    
    return ret;
}

Eigen::MatrixXd dTdx(const Eigen::VectorXd &x, const dynamical_system &ds)
{
    Eigen::MatrixXd ret(ds.xdim, ds.xdim);
    double a, b, c, d;
    unsigned int m;
    a = ds.params(0);
    b = ds.params(1);
    c = ds.params(2);
    d = ds.params(3);
    m = ds.params(4);

    ret(0,0) = 0;
    ret(0,1) = d;
    ret(1,0) = a * m * std::pow(x(0),m-1) * (x(0)*x(0)-b*b) + a * std::pow(x(0),m) * 2*x(0) + c;
    ret(1,1) = 0;

    return ret;
}