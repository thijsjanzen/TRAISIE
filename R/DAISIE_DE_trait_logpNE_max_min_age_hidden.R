#' testing fuction, for comparison with DAISIE
#' @description
#' this function calculates the likelihood of observing a singleton species on an island
#' with the trait state `i`, either non-endemic or rendered endemic by a trait change, and
#' for which only the estimated maximum age of colonization is known.
#' @export
#' @inheritParams default_params_doc
#' @examples
#' # load DAISIE package and data
#' library(DAISIE)
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
#' status <- 8
#' parameter <- list(2.546591, 2.678781, 0.009326754, 1.008583, matrix(c(0), nrow = 1), 0 )
#'
#'
#' DAISIE_DE_trait_logpNE_max_min_age_hidden(
#'   brts                  = c(23.0000,  8.3540,  8.3308),
#'   trait                 = 0,
#'   status                = 8,
#'   parameter             = parameter,
#'   num_observed_states   = 1,
#'   num_hidden_states     = 1,
#'   atol                  = 1e-15,
#'   rtol                  = 1e-15,
#'   methode               = "ode45"
#' )

DAISIE_DE_trait_logpNE_max_min_age_hidden <- function(brts,
                                                      trait,
                                                      status,
                                                      parameter,
                                                      num_observed_states,
                                                      num_hidden_states,
                                                      atol = 1e-15,
                                                      rtol = 1e-15,
                                                      methode = "ode45",
                                                      rcpp_methode = "odeint::bulirsch_stoer",
                                                      use_Rcpp = 0) {
  t0   <- brts[1]
  tmax <- brts[2]
  tmin <- brts[3]
  tp   <- 0

  # number of unique state
  n <- num_observed_states * num_hidden_states

  #########interval2 [t_p, tmin]

  m = length(parameter[[1]])

  initial_conditions2 <- get_initial_conditions2(status = status,
                                                 num_observed_states = num_observed_states,
                                                 num_hidden_states = num_hidden_states,
                                                 trait = trait,
                                                 brts = brts)



  # Time sequence for interval [tp, tmin]
  time2 <- c(tp, tmin)

  # Solve the system for interval [tp, tmin]
  solution2 <- solve_branch(interval_func = interval2,
                            initial_conditions = initial_conditions2,
                            time = time2,
                            parameter = parameter,
                            methode = methode,
                            rcpp_methode = rcpp_methode,
                            atol = atol,
                            rtol = rtol,
                            use_Rcpp = use_Rcpp)

  #########interval3 [tmin, tmax]

  # Initial conditions

  # only use second row, because the first row of solution2 is the initial state
  initial_conditions3_max_min <- c(solution2[2,][1:n],
                                   rep(0, n),       ### DE: select DE in solution2
                                   solution2[2,][(n + 1):(n + n)],         ### DM2: select DM2 in solution2
                                   solution2[2,][(n + n + 1):(n + n + n)],         ### DM3: select DM3 in solution2
                                   solution2[2,][(n + n + n + 1):(n + n + n + n)],         ### E: select E in solution2
                                   0,
                                   solution2[2,][length(solution2[2,])])                       ### DA3: select DA3 in solution2

  initial_conditions3_max_min <- matrix(initial_conditions3_max_min, nrow = 1)

  # Time sequence for interval [tmin, tmax]
  time3 <- c(tmin, tmax)

  # Solve the system for interval [tp, tmax]
  solution3 <- solve_branch(interval_func = interval3,
                            initial_conditions = initial_conditions3_max_min,
                            time = time3,
                            parameter = parameter,
                            methode = methode,
                            rcpp_methode = rcpp_methode,
                            atol = atol,
                            rtol = rtol,
                            use_Rcpp = use_Rcpp)

  #########interval4 [tmax, t0]

  # Initial conditions

  # only use second row, because the first row of solution3 is the initial state
  initial_conditions4_max_min <- c(solution3[2,][(n + 1):(n + n)],                                 ### DM1: select DM2 in solution3
                                   solution3[2,][(n + n + n + n + 1):(n + n + n + n + n)],         ### E: select E in solution3
                                   solution3[2,][length(solution3[2,]) - 1])                       ### DA1: select DA2 in solution3

  initial_conditions4_max_min <- matrix(initial_conditions4_max_min, nrow = 1)

  # Time sequence for interval [tmax, t0]
  time4 <- c(tmax, t0)

  # Solve the system for interval [tmax, t0]
  solution4 <- solve_branch(interval_func = interval4,
                            initial_conditions = initial_conditions4_max_min,
                            time = time4,
                            parameter = parameter,
                            methode = methode,
                            rcpp_methode = rcpp_methode,
                            atol = atol,
                            rtol = rtol,
                            use_Rcpp = use_Rcpp)

  # Extract log-likelihood
  Lk <- solution4[2,][length(solution4[2,])]
  logLkb <- log(Lk)
  return(logLkb)
}
