library(Rcpp)

sourceCpp(code = '
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double langevin_nll_cpp(NumericVector par, NumericVector x, NumericVector g1, NumericVector g2, double dt, IntegerVector track_start) {
    // 1. Unpack parameters
    double b00 = par[0], b01 = par[1];
    double b10 = par[2], b11 = par[3];
    double sig0 = exp(par[4]), sig1 = exp(par[5]);

    double q12 = exp(par[6]);
    double q21 = exp(par[7]);
    double lambda = q12 + q21;
    double exp_lt = exp(-lambda * dt);

    double gamma[2][2];
    if (lambda > 0.0) {
        gamma[0][0] = (q21 + q12 * exp_lt) / lambda;
        gamma[0][1] = (q12 - q12 * exp_lt) / lambda;
        gamma[1][0] = (q21 - q21 * exp_lt) / lambda;
        gamma[1][1] = (q12 + q21 * exp_lt) / lambda;
    } else {
        gamma[0][0] = 1.0; gamma[0][1] = 0.0;
        gamma[1][0] = 0.0; gamma[1][1] = 1.0;
    }

    double v[2] = {(sig0 * sig0 / 2.0) * dt, (sig1 * sig1 / 2.0) * dt};
    double m[2] = {sig0 * sqrt(dt), sig1 * sqrt(dt)};

    int n = x.length();
    double foo[2] = {0.0, 0.0};
    double lscale = 0.0;

    for (int i = 0; i < n; i++) {
        if (track_start[i] == 1) {
            foo[0] = 0.5;
            foo[1] = 0.5;
        } else {

            double drift0 = v[0] * (b00 * g1[i] + b01 * g2[i]);
            double exp_x0 = x[i-1] + drift0;
            double p0 = std::max(R::dnorm(x[i], exp_x0, m[0], 0), 1e-300);

            double drift1 = v[1] * (b10 * g1[i] + b11 * g2[i]);
            double exp_x1 = x[i-1] + drift1;
            double p1 = std::max(R::dnorm(x[i], exp_x1, m[1], 0), 1e-300);

            double v_f[2] = {0.0, 0.0};
            for(int j = 0; j < 2; j++) {
                for(int k = 0; k < 2; k++) {
                    v_f[j] += foo[k] * gamma[k][j];
                }
            }
            v_f[0] *= p0;
            v_f[1] *= p1;

            double sumfoo = v_f[0] + v_f[1];
            lscale += log(sumfoo);
            foo[0] = v_f[0] / sumfoo;
            foo[1] = v_f[1] / sumfoo;
        }
    }
    return -lscale;
}

// [[Rcpp::export]]
IntegerVector viterbi_cpp(NumericVector par, NumericVector x, NumericVector g1, NumericVector g2, double dt, IntegerVector track_start) {
    double b00 = par[0], b01 = par[1];
    double b10 = par[2], b11 = par[3];
    double sig0 = exp(par[4]), sig1 = exp(par[5]);

    double q12 = exp(par[6]);
    double q21 = exp(par[7]);
    double lambda = q12 + q21;
    double exp_lt = exp(-lambda * dt);

    double gamma[2][2];
    if (lambda > 0.0) {
        gamma[0][0] = (q21 + q12 * exp_lt) / lambda;
        gamma[0][1] = (q12 - q12 * exp_lt) / lambda;
        gamma[1][0] = (q21 - q21 * exp_lt) / lambda;
        gamma[1][1] = (q12 + q21 * exp_lt) / lambda;
    } else {
        gamma[0][0] = 1.0; gamma[0][1] = 0.0;
        gamma[1][0] = 0.0; gamma[1][1] = 1.0;
    }

    double v[2] = {(sig0 * sig0 / 2.0) * dt, (sig1 * sig1 / 2.0) * dt};
    double m[2] = {sig0 * sqrt(dt), sig1 * sqrt(dt)};

    int n = x.length();
    NumericMatrix xi(n, 2);

    for(int i = 0; i < n; i++) {
        if (track_start[i] == 1) {
            xi(i, 0) = 0.5;
            xi(i, 1) = 0.5;
        } else {
            double drift0 = v[0] * (b00 * g1[i] + b01 * g2[i]);
            double p0 = std::max(R::dnorm(x[i], x[i-1] + drift0, m[0], 0), 1e-300);

            double drift1 = v[1] * (b10 * g1[i] + b11 * g2[i]);
            double p1 = std::max(R::dnorm(x[i], x[i-1] + drift1, m[1], 0), 1e-300);

            double v0 = std::max(xi(i-1, 0) * gamma[0][0], xi(i-1, 1) * gamma[1][0]) * p0;
            double v1 = std::max(xi(i-1, 0) * gamma[0][1], xi(i-1, 1) * gamma[1][1]) * p1;

            double sumv = v0 + v1;
            xi(i, 0) = v0 / sumv;
            xi(i, 1) = v1 / sumv;
        }
    }

    IntegerVector st(n);
    for(int i = n - 1; i >= 0; i--) {
        if (i == n - 1 || track_start[i+1] == 1) {
            st[i] = (xi(i, 0) > xi(i, 1)) ? 1 : 2;
        } else {
            double v0 = gamma[0][st[i+1]-1] * xi(i, 0);
            double v1 = gamma[1][st[i+1]-1] * xi(i, 1);
            st[i] = (v0 > v1) ? 1 : 2;
        }
    }
    return st;
}
')
