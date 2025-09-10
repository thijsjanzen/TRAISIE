#' branch solving
#' @description
#' solve along branch
#' @inheritParams default_params_doc
#' @param interval_func chosen function for interval, can also be string if
#' using Rcpp
#' @param initial_conditions vector of initial conditions
#' @param time vector with two time points
#' @export
solve_branch <- function(interval_func,
                         initial_conditions,
                         time,
                         parameter,
                         trait_mainland_ancestor,
                         methode = "ode45",
                         rcpp_methode = "odeint::runge_kutta_cash_karp54",
                         atol,
                         rtol,
                         use_Rcpp = 0) {
  solution <- c()
  if (use_Rcpp <= 1) {
    parameter[[7]] <- trait_mainland_ancestor
    solution <- deSolve::ode(
      y = initial_conditions,
      times = time,
      func = interval_func,
      parms = parameter,
      method = methode,
      atol = atol,
      rtol = rtol
    )
    solution <- matrix(solution[, -1], nrow = 2)
  } else {
    interval_name <- as.character(substitute(interval_func))
    solution <- solve_branch_cpp(interval_name,
                                 initial_conditions,
                                 time,
                                 parameter,
                                 trait_mainland_ancestor,
                                 rcpp_methode,
                                 atol,
                                 rtol)
  }
  return(solution)
}
