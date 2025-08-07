#include <cstdlib>    // std::getenv, std::atoi
#include <vector>
#include <chrono>
#include "config.h"    // NOLINT [build/include_subdir]
#include <Rcpp.h>
#include <RcppParallel.h>
#include "secsse_loglik.h"    // NOLINT [build/include_subdir]

template <typename ODE>
Rcpp::List calc_ll_single_branch(std::unique_ptr<ODE> od,
                                 const Rcpp::NumericVector& states,
                                 const Rcpp::NumericVector& forTime,
                                 const std::string& method,
                                 double atol,
                                 double rtol) {
  try {

  auto t0 = std::min(forTime[0], forTime[1]);
  auto t1 = std::max(forTime[0], forTime[1]);

  auto T0 = std::chrono::high_resolution_clock::now();

  auto states_out = std::vector<double>(states.begin(), states.end());



  auto workhorse = Integrator<ODE, odeintcpp::no_normalization>(std::move(od), method, atol, rtol);

  workhorse(states_out, t0, t1);


  auto T1 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> DT = (T1 - T0);


  return Rcpp::List::create(Rcpp::Named("states") = states_out,
                            Rcpp::Named("duration") = DT.count());
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (const char* msg) {
    Rcpp::Rcout << msg << std::endl;
  } catch(...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
  return NA_REAL;
}



// [[Rcpp::export]]
Rcpp::List cpp_solve(const Rcpp::NumericVector& lambda_cs,
                                     const Rcpp::NumericVector& lambda_as,
                                     const Rcpp::NumericVector& mus,
                                     const Rcpp::NumericVector& gammas,
                                     const Rcpp::NumericMatrix& qs,
                                     const double& p,
                                     const std::string& chosen_interval,
                                     const std::string& inte_method,
                                     const Rcpp::NumericVector& init_states,
                                     const Rcpp::NumericVector& time,
                                     double atol,
                                     double rtol) {

  auto num_unique_states = lambda_cs.size();

  if (chosen_interval == "interval2") {
    return calc_ll_single_branch(std::make_unique<loglik::interval2>(lambda_cs,
                                                                     lambda_as,
                                                                     mus,
                                                                     gammas,
                                                                     qs,
                                                                     p,
                                                                     num_unique_states),
                                 init_states,
                                 time,
                                 inte_method,
                                 atol,
                                 rtol);
  }

  if (chosen_interval == "interval3") {
    return calc_ll_single_branch(std::make_unique<loglik::interval3>(lambda_cs,
                                                                     lambda_as,
                                                                     mus,
                                                                     gammas,
                                                                     qs,
                                                                     p,
                                                                     num_unique_states),
                                                                     init_states,
                                                                     time,
                                                                     inte_method,
                                                                     atol,
                                                                     rtol);
  }

  if (chosen_interval == "interval4") {
    return calc_ll_single_branch(std::make_unique<loglik::interval4>(lambda_cs,
                                                                     lambda_as,
                                                                     mus,
                                                                     gammas,
                                                                     qs,
                                                                     p,
                                                                     num_unique_states),
                                                                     init_states,
                                                                     time,
                                                                     inte_method,
                                                                     atol,
                                                                     rtol);
  }

  return NA_REAL;
}
