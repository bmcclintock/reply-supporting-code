library(Rcpp)

# Replaced source("helpers/jointLangevin.R") with the numerically stable C++ engine
sourceCpp(code = '
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double fitjoint_cpp(NumericVector par, NumericVector x,
                    NumericVector g1, NumericVector g2,
                    NumericVector h1, NumericVector h2,
                    double dt, IntegerVector track_start) {

    double beta_00 = par[0], beta_01 = par[1];
    double beta00 = par[2], beta10 = par[3];
    double beta01 = par[4], beta11 = par[5];

    double v0 = (exp(par[6]*2) / 2.0) * dt;
    double v1 = (exp(par[7]*2) / 2.0) * dt;
    double m0 = exp(par[6]) * sqrt(dt);
    double m1 = exp(par[7]) * sqrt(dt);

    double Psi = exp(par[8]);
    double delta0 = 1.0 / (1.0 + exp(-par[9]));
    double delta1 = 1.0 - delta0;

    int n = x.length();
    double foo[2] = {0.0, 0.0};
    double lscale = 0.0;

    for (int i = 0; i < n; i++) {
        if (track_start[i] == 1) {
            foo[0] = delta0;
            foo[1] = delta1;
        } else {
            double drift0 = v0 * (beta00 * g1[i] + beta01 * g2[i]);
            double drift1 = v1 * (beta10 * g1[i] + beta11 * g2[i]);

            double dx = x[i] - x[i-1];
            double p0 = std::max(R::dnorm(dx, drift0, m0, 0), 1e-300);
            double p1 = std::max(R::dnorm(dx, drift1, m1, 0), 1e-300);

            double eta0 = beta_00 + beta00 * h1[i] + beta01 * h2[i];
            double eta1 = beta_01 + beta10 * h1[i] + beta11 * h2[i];

            double L12 = Psi * exp(0.5 * (eta1 - eta0));
            double L21 = Psi * exp(0.5 * (eta0 - eta1));

            double lambda = L12 + L21;
            double exp_lt = exp(-lambda * dt);

            double G11 = 1.0, G12 = 0.0, G21 = 0.0, G22 = 1.0;
            if (lambda > 0.0) {
                G11 = (L21 + L12 * exp_lt) / lambda;
                G12 = (L12 - L12 * exp_lt) / lambda;
                G21 = (L21 - L21 * exp_lt) / lambda;
                G22 = (L12 + L21 * exp_lt) / lambda;
            }

            double v_f0 = (foo[0] * G11 + foo[1] * G21) * p0;
            double v_f1 = (foo[0] * G12 + foo[1] * G22) * p1;

            double sumfoo = v_f0 + v_f1;
            lscale += log(sumfoo);
            foo[0] = v_f0 / sumfoo;
            foo[1] = v_f1 / sumfoo;
        }
    }
    return -lscale;
}

// [[Rcpp::export]]
IntegerVector jointviterbi_cpp(NumericVector par, NumericVector x,
                               NumericVector g1, NumericVector g2,
                               NumericVector h1, NumericVector h2,
                               double dt, IntegerVector track_start) {

    double beta_00 = par[0], beta_01 = par[1];
    double beta00 = par[2], beta10 = par[3];
    double beta01 = par[4], beta11 = par[5];

    double v0 = (exp(par[6]*2) / 2.0) * dt;
    double v1 = (exp(par[7]*2) / 2.0) * dt;
    double m0 = exp(par[6]) * sqrt(dt);
    double m1 = exp(par[7]) * sqrt(dt);

    double Psi = exp(par[8]);
    double delta0 = 1.0 / (1.0 + exp(-par[9]));
    double delta1 = 1.0 - delta0;

    int n = x.length();
    NumericMatrix xi(n, 2);

    for(int i = 0; i < n; i++) {
        if (track_start[i] == 1) {
            xi(i, 0) = delta0;
            xi(i, 1) = delta1;
        } else {
            double drift0 = v0 * (beta00 * g1[i] + beta01 * g2[i]);
            double drift1 = v1 * (beta10 * g1[i] + beta11 * g2[i]);

            double dx = x[i] - x[i-1];
            double p0 = std::max(R::dnorm(dx, drift0, m0, 0), 1e-300);
            double p1 = std::max(R::dnorm(dx, drift1, m1, 0), 1e-300);

            double eta0 = beta_00 + beta00 * h1[i] + beta01 * h2[i];
            double eta1 = beta_01 + beta10 * h1[i] + beta11 * h2[i];

            double L12 = Psi * exp(0.5 * (eta1 - eta0));
            double L21 = Psi * exp(0.5 * (eta0 - eta1));

            double lambda = L12 + L21;
            double exp_lt = exp(-lambda * dt);

            double G11 = 1.0, G12 = 0.0, G21 = 0.0, G22 = 1.0;
            if (lambda > 0.0) {
                G11 = (L21 + L12 * exp_lt) / lambda;
                G12 = (L12 - L12 * exp_lt) / lambda;
                G21 = (L21 - L21 * exp_lt) / lambda;
                G22 = (L12 + L21 * exp_lt) / lambda;
            }

            double v0_val = std::max(xi(i-1, 0) * G11, xi(i-1, 1) * G21) * p0;
            double v1_val = std::max(xi(i-1, 0) * G12, xi(i-1, 1) * G22) * p1;

            double sumv = v0_val + v1_val;
            xi(i, 0) = v0_val / sumv;
            xi(i, 1) = v1_val / sumv;
        }
    }

    IntegerVector st(n);
    for(int i = n - 1; i >= 0; i--) {
        if (i == n - 1 || track_start[i+1] == 1) {
            st[i] = (xi(i, 0) > xi(i, 1)) ? 1 : 2;
        } else {
            double eta0 = beta_00 + beta00 * h1[i+1] + beta01 * h2[i+1];
            double eta1 = beta_01 + beta10 * h1[i+1] + beta11 * h2[i+1];

            double L12 = Psi * exp(0.5 * (eta1 - eta0));
            double L21 = Psi * exp(0.5 * (eta0 - eta1));

            double lambda = L12 + L21;
            double exp_lt = exp(-lambda * dt);

            double G11 = 1.0, G12 = 0.0, G21 = 0.0, G22 = 1.0;
            if (lambda > 0.0) {
                G11 = (L21 + L12 * exp_lt) / lambda;
                G12 = (L12 - L12 * exp_lt) / lambda;
                G21 = (L21 - L21 * exp_lt) / lambda;
                G22 = (L12 + L21 * exp_lt) / lambda;
            }

            if (st[i+1] == 1) {
                st[i] = (G11 * xi(i, 0) > G21 * xi(i, 1)) ? 1 : 2;
            } else {
                st[i] = (G12 * xi(i, 0) > G22 * xi(i, 1)) ? 1 : 2;
            }
        }
    }
    return st;
}
')
