#' testing fuction, for comparison with DAISIE
#' @description
#' This function compute the likelihood that all species that colonize the island
#' have gone extinct prior to the present.
#' @export
#' @inheritParams default_params_doc
#' @examples
#' #Load DAISIE package and data
#' library(DAISIE)
#' data("Galapagos_datalist")
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
#'   c(2.546591, 2.546591, 2.546591, 2.546591),
#'   c(2.678781, 2.678781, 2.678781, 2.678781),
#'   c(0.009326754, 0.009326754, 0.009326754, 0.009326754),
#'   c(1.008583, 1.008583, 1.008583, 1.008583),
#'   matrix(c(
#'     rep(0, 16)
#'   ), nrow = 4),
#'   0
#' )
#' # Compute log‐likelihood under the DE‐trait model
#' DAISIE_DE_trait_logp0(
#'   datalist            = datalist,
#'   parameter           = parameter,
#'   num_observed_states     = 1,
#'   num_hidden_states       = 1,
#'   atol                = 1e-5,
#'   rtol                = 1e-5,
#'   methode             = "ode45"
#' )


DAISIE_DE_trait_logp0 <- function(datalist,
                                  parameter,
                                  atol = 1e-15,
                                  rtol = 1e-15,
                                  num_observed_states,
                                  num_hidden_states,
                                  methode = "ode45",
                                  rcpp_methode = "odeint::bulirsch_stoer",
                                  use_Rcpp = 0) {

  n <- num_observed_states * num_hidden_states
  t0 <- datalist[[1]]$island_age
  tp <- 0
  #########interval4 [t_p, t_0]

  initial_conditions40 <- c(rep(0, n), ### DM1
                           rep(0, n), ### E
                           1)                         ### DA1

  # Time sequence for interval [tp, t0]
  time4 <- c(tp, t0)

  # Solve the system for interval [tp, t1]
  solution4 <- solve_branch(interval_func = interval4,
                            initial_conditions = initial_conditions40,
                            time = time4,
                            parameter = parameter,
                            methode = methode,
                            rcpp_methode = rcpp_methode,
                            atol = atol,
                            rtol = rtol,
                            use_Rcpp = use_Rcpp)

  # Extract log-likelihood
  Lk <- solution4[2, ][length(solution4[2, ])]
  logLkb <- log(Lk)
  return(logLkb)
}
