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
#' datalist <- Galapagos_datalist
#' datalist[[1]]$M0 <- 500
#' datalist[[1]]$M1 <- 400
#'
#' parameter <- list(
#'   c(2.546591, 2.546591, 2.546591, 2.546591),
#'   c(2.678781, 2.678781, 2.678781, 2.678781),
#'   c(0.009326754, 0.009326754, 0.009326754, 0.009326754),
#'   c(1.008583, 1.008583, 1.008583, 1.008583),
#'   matrix(c(
#'     0,    0,    0,  0,
#'     0,    0,    0.00,0.00,
#'     rep(0, 8)
#'   ), nrow = 4),
#'   1
#' )
#'
#' parameter <- list(
#'   c(2.546591, 1.2, 1, 0.2),
#'   c(2.678781, 2, 1.9, 3),
#'   c(0.009326754, 0.003, 0.002, 0.2),
#'   c(1.008583, 1, 2, 1.5),
#'   matrix(c(
#'     0,    .001,    0.005,  0,
#'     .001,    0,    0.000,0.005,
#'     0.005,    000,    0,  0.005,
#'     0,   0.005,  0.005,0.00
#'   ), nrow = 4),
#'   1
#' )
#'
#' DAISIE_DE_trait_logp0(
#'   datalist,
#'   parameter               = parameter,
#'   num_observed_states     = 2,
#'   num_hidden_states       = 2,
#'   atol                    = 1e-15,
#' trait_mainland_ancestor   =  NA,
#'   rtol                    = 1e-15,
#'   methode                 = "ode45",
#'   rcpp_methode ="odeint::runge_kutta_cash_karp54",
#'   use_Rcpp                = 2)
#'


DAISIE_DE_trait_logp0 <- function(
    datalist,
    parameter,
    atol = 1e-15,
    rtol = 1e-15,
    num_observed_states,
    num_hidden_states,
    trait_mainland_ancestor = NA,
    methode = "ode45",
    rcpp_methode ="odeint::runge_kutta_cash_karp54",
    use_Rcpp = 2) {

  loglik_func <- function(i) {
    trait_mainland_ancestor_extended <- rep(0,num_observed_states * num_hidden_states)
    trait_mainland_ancestor_extended[i] <- 1 #set only the trait of interest to 1

    Lk <- DAISIE_DE_trait_logp0_core(datalist,
                                     parameter,
                                     atol = 1e-15,
                                     rtol = 1e-15,
                                     num_observed_states,
                                     num_hidden_states,
                                     trait_mainland_ancestor = trait_mainland_ancestor_extended,
                                     methode = "ode45",
                                     rcpp_methode = rcpp_methode,
                                     use_Rcpp = use_Rcpp)
    return(Lk)
  }
  indices <-  seq_len(num_observed_states * num_hidden_states)
  Lk_vec <- sapply(indices, loglik_func)

  ## added !all(is.na(trait_mainland_ancestor)) because when trait_mainland_ancestor = NA,  length(trait_mainland_ancestor) = length(trait_mainland_ancestor_extended) = 1
  if(!any(is.na(trait_mainland_ancestor)) && length(trait_mainland_ancestor) == num_observed_states * num_hidden_states) { #this is the case where a full probability distribution is specified across all observed and hidden states
    weights <- trait_mainland_ancestor/sum(trait_mainland_ancestor)
  } else {

    # Determine probabilities for observed states
    if (all(is.numeric(trait_mainland_ancestor))) {
      # User provided trait_mainland_ancestor, e.g. c(1, 0)
      probs <- trait_mainland_ancestor
      # Replicate each probability across the hidden states
      s <- unlist(lapply(probs, function(p) rep(p, num_hidden_states)))



    } else {

      # Nothing provided -> compute probabilities from the transition matrix
      Q =  parameter[[5]]
      diag(Q ) <- 0
      diag(Q ) <- -rowSums(Q)

      pi <- pracma::null(t(Q))
      if (pi[which.max(abs(pi))] < 0) pi <- -pi

      pi <-  pmax (rep (0, 4) , pi)
      s <- pi
    }
    # Normalize weights
    if (sum(s) == 0) {
      weights <- s
    } else {
      weights <- s / sum(s)
    }

  }

  log_Lk <- log(sum(Lk_vec * weights))
  return( list (loglik = log_Lk, lik_states = Lk_vec, weights = weights))
}

DAISIE_DE_trait_logp0_core <- function(datalist,
                                       parameter,
                                       atol = 1e-15,
                                       rtol = 1e-15,
                                       num_observed_states,
                                       num_hidden_states,
                                       trait_mainland_ancestor= NA,
                                       methode = "ode45",
                                       rcpp_methode =
                                         "odeint::runge_kutta_cash_karp54",
                                       use_Rcpp = 0) {

  n <- num_observed_states * num_hidden_states
  t0 <- datalist[[1]]$island_age
  tp <- 0

  #########interval4 [t_p, t_0]

  initial_conditions40 <- c(rep(0, n),  ### DM1
                            rep(0, n),  ### E
                            1)          ### DA1

  # Time sequence for interval [tp, t0]
  time4 <- c(tp, t0)

  # Solve the system for interval [tp, t1]
  solution4 <- solve_branch(interval_func = interval4,
                            initial_conditions = initial_conditions40,
                            time = time4,
                            parameter = parameter,
                            trait_mainland_ancestor = trait_mainland_ancestor,
                            methode = methode,
                            rcpp_methode = rcpp_methode,
                            atol = atol,
                            rtol = rtol,
                            use_Rcpp = use_Rcpp)

  # Extract log-likelihood
  Lk <- solution4[2, ][length(solution4[2, ])]

  return(Lk)
}
