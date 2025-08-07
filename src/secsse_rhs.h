//  Copyright (c) 2021 - 2023, Thijs Janzen
//  Copyright (c) 2023, Hanno Hildenbrandt
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)

#pragma once
#include <Rcpp.h>
#include <RcppParallel.h>
#include <type_traits>
#include <vector>

#include "sq_matrix.h"

namespace loglik {

template <typename T>
using rvector = RcppParallel::RVector<T>;

template <typename T>
using rmatrix = RcppParallel::RMatrix<T>;

class interval1 {

  const rvector<const double> lc_; // cladogenesis rates
  const rvector<const double> m_; //  extinction rates

  const rvector<const double> la_; // anagenesis rates
  const sq_matrix q_; // transition rates
  const double p_;

  const size_t n_; // number of unique states

  const std::vector<double> t_vec;


  const sq_matrix g_mat_; //  colonization rates
  const rvector<const double> g_;

public:

  // constructor
  interval1(const Rcpp::NumericVector& lc,
            const Rcpp::NumericVector& la,
            const Rcpp::NumericVector& m,
            const Rcpp::NumericVector& g,
            const Rcpp::NumericMatrix& q,
            const double p,
            const size_t n)
    : lc_(lc),
      m_(m),
      la_(la),
      q_(q),
      p_(p),
      n_(n),
      t_vec(q_.row_sums()),
      g_mat_(g, n_),
      g_(g) {
  }

  size_t size() const noexcept {
    // (DE + DM3 + E) * n + DA3
    return 3 * n_ + 1;
  }

  void mergebranch(const std::vector<double>& N,
                   const std::vector<double>& M,
                   std::vector<double>& out) const {
    /*
     out[DE_0]  = lc_[0] * N[DE_0] * M[DE_0];
     out[DE_1]  = lc_[1] * N[DE_1] * M[DE_1];

     out[DM3_0] = N[DM3_0];
     out[DM3_1] = N[DM3_1];

     out[E_0]   = N[E_0];
     out[E_1]   = N[E_1];
     out[DA_3]  = N[DA_3];
     */

    out = N;

    for (size_t i = 0; i < n_; ++i) {
      // out[DE_0]  = lc_[0] * N[DE_0] * M[DE_0]
      out[i] = lc_[i] * N[i] * M[i];
    }
  }

  // this is the dx/dt calculation // true rhs that gets integrated
  // along the branches
  void operator()(const std::vector<double>& x,
                  std::vector<double>& dxdt,
                  const double /* t */) const
  {
    auto DA3 = x.back();
    auto gamma_nonself = g_mat_.non_self();

    auto DE  = vector_view_t<const double>(x.data() , n_);
    auto DM3 = vector_view_t<const double>(x.data() + n_, n_);
    auto E   = vector_view_t<const double>(x.data() + n_ + n_, n_);

    auto q_mult_E   = q_ * E;
    auto q_mult_DE  = q_ * DE;
    auto q_mult_DM3 = q_ * DM3;

    auto g_row_sums = g_mat_.row_sums();

    double s_g_DM3 = 0.0;

    for (size_t i = 0; i < n_; ++i) {
      auto lambda_c_mu_t_vec_sum = lc_[i] + m_[i] + t_vec[i];

      // DE
      dxdt[i] = -(lambda_c_mu_t_vec_sum) * DE[i] +
                  2 * lc_[i] * DE[i] * E[i] +
                  q_mult_DE[i];
      // DM3
      dxdt[i + n_] = -(lambda_c_mu_t_vec_sum + gamma_nonself[i] + la_[i]) * DM3[i] +
                  (m_[i] + la_[i] * E[i] + lc_[i] * E[i] * E[i] + p_ * q_mult_E[i]) * DA3 +
                  (1 - p_) * q_mult_DM3[i] +
                  gamma_nonself[i] * DM3[i];
      // E
      dxdt[i + n_ + n_] = m_[i] - (lambda_c_mu_t_vec_sum) * E[i] +
                       lc_[i] * E[i] * E[i] +
                       q_mult_E[i];

      s_g_DM3 += g_[i] * DM3[i];
    }

    // DA3
    auto a = std::accumulate(g_.begin(), g_.end(), 0.0);
    dxdt.back() = -a * DA3 + s_g_DM3;
  }
};

class interval2 {
  const rvector<const double> lc_; // cladogenesis rates
  const rvector<const double> m_; //  extinction rates

  const rvector<const double> la_; // anagenesis rates
  const sq_matrix q_; // transition rates
  const double p_;

  const size_t n_; // number of unique states

  const std::vector<double> t_vec;


  const sq_matrix g_mat_; //  colonization rates
  const rvector<const double> g_;

public:

  // constructor
  interval2(const Rcpp::NumericVector& lc,
            const Rcpp::NumericVector& la,
            const Rcpp::NumericVector& m,
            const Rcpp::NumericVector& g,
            const Rcpp::NumericMatrix& q,
            const double p,
            const size_t n)
    : lc_(lc),
      m_(m),
      la_(la),
      q_(q),
      p_(p),
      n_(n),
      t_vec(q_.row_sums()),
      g_mat_(g, n_),
      g_(g) {
  }

  size_t size() const noexcept {
    // (DE + DM3 + E) * n + DA3
    return 4 * n_ + 1;
  }

  // this is the dx/dt calculation // true rhs that gets integrated
  // along the branches
  void operator()(const std::vector<double>& x,
                std::vector<double>& dxdt,
                const double /* t */) const
  {
    // substitute with your own code below:

    // vector x is:
    // [
    // DE_0, DE_1, ... , DE_n
    // DM3_0, DM3_1, ..., DM3_n
    // E_0, E_1, ..., E_n
    // DA3
    // ]

    auto DA3 = x.back();
    auto gamma_nonself = g_mat_.non_self();

    auto DE  = vector_view_t<const double>(x.data() + 0 * n_, n_);
    auto DM2 = vector_view_t<const double>(x.data() + 1 * n_, n_);
    auto DM3 = vector_view_t<const double>(x.data() + 2 * n_, n_);
    auto E   = vector_view_t<const double>(x.data() + 3 * n_, n_);

    auto q_mult_E   = q_ * E;
    auto q_mult_DE  = q_ * DE;
    auto q_mult_DM2 = q_ * DM2;
    auto q_mult_DM3 = q_ * DM3;

    auto g_row_sums = g_mat_.row_sums();

    double s_g_DM3 = 0.0;

    for (size_t i = 0; i < n_; ++i) {
      auto lambda_c_mu_t_vec_sum = lc_[i] + m_[i] + t_vec[i];

      // DE
      dxdt[i] = -(lambda_c_mu_t_vec_sum) * DE[i] +
                  2 * lc_[i] * DE[i] * E[i] +
                  q_mult_DE[i];

      // DM2
      dxdt[i + n_] = -(lambda_c_mu_t_vec_sum + g_[i] + la_[i]) * DM2[i] +
                      (la_[i] * DE[i] + 2 * lc_[i] * DE[i] * E[i] + p_ * q_mult_DE[i]) * DA3 +
                      (1 - p_) * q_mult_DM2[i] + gamma_nonself[i] * DM2[i];

      // DM3
      dxdt[i + n_ + n_] = -(lambda_c_mu_t_vec_sum + gamma_nonself[i] + la_[i]) * DM3[i] +
        (m_[i] + la_[i] * E[i] + lc_[i] * E[i] * E[i] + p_ * q_mult_E[i]) * DA3 +
        (1 - p_) * q_mult_DM3[i] +
        gamma_nonself[i] * DM3[i];
      // E
      dxdt[i + n_ + n_ + n_] = m_[i] - (lambda_c_mu_t_vec_sum) * E[i] +
        lc_[i] * E[i] * E[i] +
        q_mult_E[i];

      s_g_DM3 += g_[i] * DM3[i];
    }

    // DA3
    auto a = std::accumulate(g_.begin(), g_.end(), 0.0);
    dxdt.back() = -a * DA3 + s_g_DM3;
  }
};

class interval3 {
  const rvector<const double> lc_; // cladogenesis rates
  const rvector<const double> m_; //  extinction rates

  const rvector<const double> la_; // anagenesis rates
  const sq_matrix q_; // transition rates
  const double p_;

  const size_t n_; // number of unique states

  const std::vector<double> t_vec;


  const sq_matrix g_mat_; //  colonization rates
  const rvector<const double> g_;

public:

  // constructor
  interval3(const Rcpp::NumericVector& lc,
            const Rcpp::NumericVector& la,
            const Rcpp::NumericVector& m,
            const Rcpp::NumericVector& g,
            const Rcpp::NumericMatrix& q,
            const double p,
            const size_t n)
    : lc_(lc),
      m_(m),
      la_(la),
      q_(q),
      p_(p),
      n_(n),
      t_vec(q_.row_sums()),
      g_mat_(g, n_),
      g_(g) {
  }

  size_t size() const noexcept {
    // (DE + DM3 + E) * n + DA3
    return 5 * n_ + 2;
  }
  void operator()(const std::vector<double>& x,
                std::vector<double>& dxdt,
                const double /* t */) const
  {

    auto DA3 = x.back();
    auto DA2 = x[x.size() - 2];
    auto gamma_nonself = g_mat_.non_self();

    auto DE  = vector_view_t<const double>(x.data() , n_);
    auto DM1 = vector_view_t<const double>(x.data() + 1 * n_, n_);
    auto DM2 = vector_view_t<const double>(x.data() + 2 * n_, n_);
    auto DM3 = vector_view_t<const double>(x.data() + 3 * n_, n_);
    auto E   = vector_view_t<const double>(x.data() + 4 * n_, n_);

    auto q_mult_E   = q_ * E;
    auto q_mult_DE  = q_ * DE;
    auto q_mult_DM1 = q_ * DM1;
    auto q_mult_DM2 = q_ * DM2;
    auto q_mult_DM3 = q_ * DM3;

    auto g_row_sums = g_mat_.row_sums();

    double s_g_DM3 = 0.0;
    double s_g_DM2 = 0.0;

    for (size_t i = 0; i < n_; ++i) {
      auto lambda_c_mu_t_vec_sum = lc_[i] + m_[i] + t_vec[i];

      // DE
      dxdt[i] = -(lambda_c_mu_t_vec_sum) * DE[i] +
        2 * lc_[i] * DE[i] * E[i] +
        q_mult_DE[i];

      // DM1
      dxdt[i + 1 * n_] = -(lambda_c_mu_t_vec_sum + g_[i] + la_[i]) * DM1[i] +
                         (m_[i] + la_[i] * E[i] + lc_[i] * E[i] * E[i] + p_ * q_mult_E[i]) * DA2 +
                         (1 - p_) * q_mult_DM1[i] + g_[i] * DM2[i];
      // DM2
      dxdt[i + 2 * n_] = -(lambda_c_mu_t_vec_sum + gamma_nonself[i] + la_[i]) * DM2[i] +
        (m_[i] + la_[i] * E[i] + lc_[i] * E[i] * E[i] + p_ * q_mult_E[i]) * DA2 +
        (la_[i] * DE[i] + 2 * lc_[i] * DE[i]  + p_ * q_mult_DE[i]) * DA3 +
        (1 - p_) * q_mult_DM2[i] + gamma_nonself[i] * DM2[i];

      // DM3
      dxdt[i + 3 * n_] = -(lambda_c_mu_t_vec_sum + gamma_nonself[i] + la_[i]) * DM3[i] +
        (m_[i] + la_[i] * E[i] + lc_[i] * E[i] * E[i] + p_ * q_mult_E[i]) * DA3 +
        (1 - p_) * q_mult_DM3[i] +
        gamma_nonself[i] * DM3[i];

      // E
      dxdt[i + 4 * n_] = m_[i] - (lambda_c_mu_t_vec_sum) * E[i] +
        lc_[i] * E[i] * E[i] +
        q_mult_E[i];

      s_g_DM3 += g_[i] * DM3[i];
      s_g_DM2 += g_[i] * DM2[i];
    }

    // DA3
    auto a = std::accumulate(g_.begin(), g_.end(), 0.0);
    dxdt.back() = -a * DA3 + s_g_DM3;

    // DA2
    dxdt[dxdt.size() - 2] = -a * DA2 + s_g_DM2;
  }
};

class interval4 {
  const rvector<const double> lc_; // cladogenesis rates
  const rvector<const double> m_; //  extinction rates

  const rvector<const double> la_; // anagenesis rates
  const sq_matrix q_; // transition rates
  const double p_;

  const size_t n_; // number of unique states

  const std::vector<double> t_vec;


  const sq_matrix g_mat_; //  colonization rates
  const rvector<const double> g_;

public:

  // constructor
  interval4(const Rcpp::NumericVector& lc,
            const Rcpp::NumericVector& la,
            const Rcpp::NumericVector& m,
            const Rcpp::NumericVector& g,
            const Rcpp::NumericMatrix& q,
            const double p,
            const size_t n)
    : lc_(lc),
      m_(m),
      la_(la),
      q_(q),
      p_(p),
      n_(n),
      t_vec(q_.row_sums()),
      g_mat_(g, n_),
      g_(g) {
  }

  size_t size() const noexcept {
    return 2 * n_ + 1;
  }

  void operator()(const std::vector<double>& x,
                std::vector<double>& dxdt,
                const double /* t */) const
  {

    auto DA1 = x.back();
    auto gamma_nonself = g_mat_.non_self();

    auto DM1  = vector_view_t<const double>(x.data() , n_);
    auto E   = vector_view_t<const double>(x.data() + n_, n_);

    auto q_mult_E   = q_ * E;
    auto q_mult_DM1  = q_ * DM1;

    auto g_row_sums = g_mat_.row_sums();

    double s_g_DM1 = 0.0;

    for (size_t i = 0; i < n_; ++i) {
      auto lambda_c_mu_t_vec_sum = lc_[i] + m_[i] + t_vec[i];

      // DM1
      dxdt[i] = -(lambda_c_mu_t_vec_sum + gamma_nonself[i] + la_[i]) * DM1[i] +
                (m_[i] + la_[i] * E[i] + lc_[i] * E[i] * E[i] + p_ * q_mult_E[i]) * DA1 +
                (1 - p_) * q_mult_DM1[i] + gamma_nonself[i] * DM1[i];

      // E
      dxdt[i + n_] = m_[i] - (lambda_c_mu_t_vec_sum) * E[i] +
        lc_[i] * E[i] * E[i] +
        q_mult_E[i];

      s_g_DM1 += g_[i] * DM1[i];
    }

    // DA3
    auto a = std::accumulate(g_.begin(), g_.end(), 0.0);
    dxdt.back() = -a * DA1 + s_g_DM1;
  }
};

} // namespace secsse
