#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// Update the two-age edge density and integrate its old-state sources in
// storage order (degree k varies fastest). Do not change the R discretization.
// [[Rcpp::export]]
List pde_ii_step(const NumericVector& II, const NumericVector& SI,
                 const NumericMatrix& PT, const NumericVector& Lambda,
                 const NumericMatrix& SI_in, const NumericMatrix& survival,
                 double common, double beta, double alpha, double gamma,
                 double delta_p, double dt, double da, int active_ages) {
  IntegerVector dim = II.attr("dim");
  if (dim.size() != 4 || dim[0] != dim[1] || dim[2] != dim[3])
    stop("II must have dimensions K x K x A x A");
  const int K = dim[0], A = dim[2];
  const R_xlen_t KK = static_cast<R_xlen_t>(K) * K;
  IntegerVector si_dim = SI.attr("dim");
  if (K < 1 || A < 2 || active_ages < 1 || active_ages > A ||
      si_dim.size() != 3 || si_dim[0] != K || si_dim[1] != K || si_dim[2] != A ||
      PT.nrow() != K || PT.ncol() != A || Lambda.size() != K ||
      SI_in.nrow() != K || SI_in.ncol() != K ||
      survival.nrow() != K || survival.ncol() != A)
    stop("Incompatible dimensions for the II update");
  NumericVector next(II.size()); // Zero outside the occupied age range.
  next.attr("dim") = dim;
  NumericVector W(KK * A), Q(KK * A);
  W.attr("dim") = IntegerVector::create(K, K, A);
  Q.attr("dim") = IntegerVector::create(K, K, A);
  const double common2 = common * common;
  for (int b = 0; b < active_ages; ++b) {
    if (b % 16 == 0) checkUserInterrupt();
    for (int a = 0; a < active_ages; ++a) {
      const R_xlen_t age_offset = KK * (a + static_cast<R_xlen_t>(A) * b);
      for (int j = 0; j < K; ++j) {
        const double trace_j = PT(j, b);
        const double survive_j2 = survival(j, b) * survival(j, b);
        for (int k = 0; k < K; ++k) {
          const R_xlen_t pair = k + static_cast<R_xlen_t>(K) * j;
          const double value = II[pair + age_offset];
          if (!std::isfinite(value) || value < 0) stop("Invalid old II density");
          const R_xlen_t row = pair + KK * a;
          W[row] += value;
          Q[row] += value * trace_j;
          if (a < A-1 && b < A-1) {
            const double updated = value * common2 *
              ((survival(k, a) * survival(k, a)) * survive_j2);
            if (!std::isfinite(updated) || updated < 0) stop("Invalid new II density");
            next[pair + age_offset + KK * (A + 1)] = updated;
          }
        }
      }
    }
  }
  for (R_xlen_t i = 0; i < W.size(); ++i) { W[i] *= da; Q[i] *= da; }

  // The two incoming face strips and their shared cell use precisely the
  // same midpoint and half-cell rules as the R reference implementation.
  const int incoming_ages = std::min(active_ages, A-1);
  for (int j = 0; j < K; ++j) {
    for (int k = 0; k < K; ++k) {
      const R_xlen_t pair = k + static_cast<R_xlen_t>(K) * j;
      const R_xlen_t reverse = j + static_cast<R_xlen_t>(K) * k;
      const double factor_k = beta * (1 + k * Lambda[k]);
      const double factor_j = beta * (1 + j * Lambda[j]);
      for (int a = 0; a < incoming_ages; ++a) {
        const double loss = beta + alpha + gamma + delta_p * j * PT(j, a) +
          beta * k * Lambda[k];
        const double reverse_loss = beta + alpha + gamma + delta_p * k * PT(k, a) +
          beta * j * Lambda[j];
        const double first = factor_k * SI[pair + KK*a] * std::exp(-loss*dt/2) *
          common * survival(k, 0) * survival(j, a);
        const double second = factor_j * SI[reverse + KK*a] * std::exp(-reverse_loss*dt/2) *
          common * survival(k, a) * survival(j, 0);
        if (!std::isfinite(first) || first < 0 || !std::isfinite(second) || second < 0)
          stop("Invalid II boundary density");
        next[pair + KK * A * (a+1)] = first;
        next[pair + KK * (a+1)] = second;
      }
      const double corner = 0.5 * (factor_k * SI_in(k,j) + factor_j * SI_in(j,k)) *
        common * survival(k,0) * survival(j,0);
      if (!std::isfinite(corner) || corner < 0) stop("Invalid II corner density");
      next[pair] = corner;
    }
  }
  return List::create(_["II"] = next, _["W"] = W, _["Q"] = Q);
}

// Check symmetry without materializing a second four-dimensional array.
// [[Rcpp::export]]
double pde_ii_symmetry(const NumericVector& II) {
  IntegerVector dim = II.attr("dim");
  if (dim.size() != 4 || dim[0] != dim[1] || dim[2] != dim[3])
    stop("II must have dimensions K x K x A x A");
  const int K = dim[0], A = dim[2];
  const R_xlen_t KK = static_cast<R_xlen_t>(K) * K;
  double maximum = 1, difference = 0;
  for (int b = 0; b < A; ++b) {
    if (b % 16 == 0) checkUserInterrupt();
    for (int a = 0; a < A; ++a) {
      for (int j = 0; j < K; ++j) {
        for (int k = 0; k < K; ++k) {
          double value = II[k + K*j + KK*(a + static_cast<R_xlen_t>(A)*b)];
          double reverse = II[j + K*k + KK*(b + static_cast<R_xlen_t>(A)*a)];
          maximum = std::max(maximum, value);
          difference = std::max(difference, std::abs(value-reverse));
        }
      }
    }
  }
  return difference/maximum;
}
