// Copyright 2023 Hanno Hildenbrandt
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)


#include <cstdlib>    // std::getenv, std::atoi
#include <vector>
#include <chrono>
#include "config.h"    // NOLINT [build/include_subdir]
#include <Rcpp.h>
#include <RcppParallel.h>
#include "secsse_loglik.h"    // NOLINT [build/include_subdir]


// probably the cleanest way to retrieve RcppParallel's concurrency setting
// set by RcppParallel::setThreadOptions(numThreads)
size_t get_rcpp_num_threads() {
  auto* nt_env = std::getenv("RCPP_PARALLEL_NUM_THREADS");
  return (nullptr == nt_env)
    ? tbb::task_arena::automatic  // -1
  : static_cast<size_t>(std::atoi(nt_env));
}

template <typename ODE>
Rcpp::List calc_ll(std::unique_ptr<ODE> od,
                   const Rcpp::IntegerVector& ances,
                   const Rcpp::NumericMatrix& states,
                   const Rcpp::NumericMatrix& forTime,
                   const std::string& method,
                   double atol,
                   double rtol,
                   bool see_states,
                   bool use_normalization) {
  auto num_threads = get_rcpp_num_threads();
  auto global_control = tbb::global_control(tbb::global_control::max_allowed_parallelism, num_threads);

  auto T0 = std::chrono::high_resolution_clock::now();
  std::vector<std::vector<double>> tstates{};
  for (int i = 0; i < states.nrow(); ++i) {
    tstates.emplace_back(states.row(i).begin(), states.row(i).end());
  }
  const auto phy_edge = make_phy_edge_vector(loglik::rmatrix<const double>(forTime));
  auto inodes = find_inte_nodes(phy_edge, loglik::rvector<const int>(ances), tstates);

  calc_ll_res ll_res;
  if (use_normalization) {
    ll_res  = calc_ll(Integrator<ODE, odeintcpp::normalize>(      std::move(od), method, atol, rtol), inodes, tstates);
  } else {
    ll_res = calc_ll(Integrator<ODE, odeintcpp::no_normalization>(std::move(od), method, atol, rtol), inodes, tstates);
  }


  auto T1 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> DT = (T1 - T0);
  Rcpp::NumericMatrix states_out;
  if (see_states) {
    // R side expect full states back.
    states_out = Rcpp::NumericMatrix(states.nrow(), states.ncol());
    for (int i = 0; i < states.nrow(); ++i) {
      std::copy(std::begin(tstates[i]), std::end(tstates[i]),
                states_out.row(i).begin());
    }
  }
  return Rcpp::List::create(Rcpp::Named("loglik") = ll_res.loglik,
                            Rcpp::Named("node_M") = ll_res.node_M,
                            Rcpp::Named("merge_branch") = ll_res.merge_branch,
                            Rcpp::Named("states") = states_out,
                            Rcpp::Named("duration") = DT.count());
}

// [[Rcpp::export]]
Rcpp::List calc_ll_cpp(const Rcpp::IntegerVector& ances,
                       const Rcpp::NumericMatrix& states,
                       const Rcpp::NumericMatrix& forTime,
                       const Rcpp::NumericVector& lambda_cs,
                       const Rcpp::NumericVector& lambda_as,
                       const Rcpp::NumericVector& mus,
                       const Rcpp::NumericVector& gammas,
                       const Rcpp::NumericMatrix& qs,
                       const double& p,
                       const std::string& method,
                       double atol,
                       double rtol,
                       bool see_states,
                       bool use_normalization) {
  try {
    size_t num_unique_states = (states.ncol() - 1) / 3;

    return calc_ll(std::make_unique<loglik::interval1>(lambda_cs,
                                                       lambda_as,
                                                       mus,
                                                       gammas,
                                                       qs,
                                                       p,
                                                       num_unique_states),
                                                       ances, states, forTime, method, atol, rtol, see_states, use_normalization);
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (const char* msg) {
    Rcpp::Rcout << msg << std::endl;
  } catch(...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
  return NA_REAL;
}
