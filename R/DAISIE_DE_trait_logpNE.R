
#' This function calculates the likelihood of observing a non-endemic lineage with specified species trait states,
#'
#'
#' @inheritParams default_params_doc
#'
#' @export
#'
#' @examples
#' library(DAISIE)
#' library(secsse)
#' data("Galapagos_datalist")
#' datalist <- Galapagos_datalist
#'
#' i <- 3
#' brts <- datalist[[i]]$branching_times
#' trait <- 0
#'
#'
#' parameter <- list(
#'   c(2.546591, 1.2, 1, 0.2),
#'   c(2.678781, 2, 1.9, 3),
#'   c(0.009326754, 0.003, 0.002, 0.2),
#'   c(1.008583, 1, 2, 1.5),
#'   matrix(c(
#'     0,    1,    0.5,  0,
#'     0,    0,    0.002,0.005,
#'     rep(0, 8)
#'   ), nrow = 4),
#'   0
#' )
#'
#' parameter <- list(
#'   c(2.546591, 2.546591),
#'   c(2.678781, 2.678781),
#'   c(0.009326754, 0.009326754),
#'   c(1.008583, 1.008583),
#'   matrix(c(
#'     rep(0, 4)
#'   ), nrow = 2),
#'   0
#' )
#' brts <- datalist[[2]]$branching_times
#' DAISIE_DE_trait_logpNE(
#'   brts                    = brts,
#'   trait                   = trait,
#'   status                  = 1,
#'   parameter               = parameter,
#'   num_observed_states     = 2,
#'   num_hidden_states       = 1,
#'   atol                    = 1e-15,
#'   rtol                    = 1e-15,
#'   methode                 = "ode45"
#' )
DAISIE_DE_trait_logpNE <- function(brts,
                                   status,
                                   trait,
                                   trait_mainland_ancestor = FALSE,
                                   num_observed_states,
                                   num_hidden_states,
                                   parameter,
                                   atol  = 1e-15,
                                   rtol  = 1e-15,
                                   methode                 = "ode45",
                                   rcpp_methode = "odeint::bulirsch_stoer",
                                   use_Rcpp = 0) {

  check_arguments(brts = brts,
                  parameter = parameter,
                  phy = 1,
                  traits = trait,
                  num_observed_states = num_observed_states,
                  num_hidden_states = num_hidden_states,
                  status = status,
                  sampling_fraction = 1)

  # Unpack times from brts
  t0   <- brts[1]
  tmax <- brts[2]
  t1   <- brts[2]

  tp   <- 0


  # Time intervals

  time2 <- c(tp, t1)
  time3 <- c(tp, tmax)
  time4 <- c(tmax, t0)

  # Solve for interval [tp, t2] (stem phase)

  # Run appropriate sequence of intervals
  if (status == 4) {
    initial_conditions2 <- get_initial_conditions2(status = status,
                                                   trait = trait,
                                                   num_observed_states = num_observed_states,
                                                   num_hidden_states = num_hidden_states,
                                                   brts = brts,
                                                   trait_mainland_ancestor = trait_mainland_ancestor)
    solution2 <- solve_branch(interval_func = interval2,
                              initial_conditions = initial_conditions2,
                              time = time2,
                              parameter = parameter,
                              methode = methode,
                              rcpp_methode = rcpp_methode,
                              atol = atol,
                              rtol = rtol,
                              use_Rcpp = use_Rcpp)

    initial_conditions4 <- get_initial_conditions4(status = status,
                                                   solution = solution2,
                                                   parameter = parameter,
                                                   trait_mainland_ancestor = trait_mainland_ancestor,
                                                   num_observed_states = num_observed_states,
                                                   num_hidden_states = num_hidden_states)
    solution4 <- solve_branch(interval_func = interval4,
                              initial_conditions = initial_conditions4,
                              time = time4,
                              parameter = parameter,
                              methode = methode,
                              rcpp_methode = rcpp_methode,
                              atol = atol,
                              rtol = rtol,
                              use_Rcpp = use_Rcpp)
  }

  if (status == 1) {
    initial_conditions3 <- get_initial_conditions3(status = status,
                                                   num_observed_states = num_observed_states,
                                                   num_hidden_states = num_hidden_states,
                                                   trait = trait)
    solution3 <- solve_branch(interval_func = interval3,
                              initial_conditions = initial_conditions3,
                              time = time3,
                              parameter = parameter,
                              methode = methode,
                              rcpp_methode = rcpp_methode,
                              atol = atol,
                              rtol = rtol,
                              use_Rcpp = use_Rcpp)

    initial_conditions4 <- get_initial_conditions4(status = status,
                                                   solution = solution3,
                                                   parameter = parameter,
                                                   trait_mainland_ancestor = trait_mainland_ancestor,
                                                   num_observed_states = num_observed_states,
                                                   num_hidden_states = num_hidden_states)
    solution4 <- solve_branch(interval_func = interval4,
                              initial_conditions = initial_conditions4,
                              time = time4,
                              parameter = parameter,
                              methode = methode,
                              rcpp_methode = rcpp_methode,
                              atol = atol,
                              rtol = rtol,
                              use_Rcpp = use_Rcpp)
  }

  # Extract log-likelihood from final solution
  Lk <- solution4[2, length(solution4[2, ])]
  logLkb <- log(Lk)
  return(logLkb)
}
