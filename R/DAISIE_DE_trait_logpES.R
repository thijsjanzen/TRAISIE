#' Testing function for comparison with DAISIE
#'
#' @description
#' This function calculates the likelihood of observing a clade with specified species trait states,
#' given known colonization time. It is designed for comparison with DAISIE-based models.
#'
#' @inheritParams default_params_doc
#'
#' @export
#'
#' @examples
#' library(DAISIE)
#' data("Galapagos_datalist")
#' datalist <- Galapagos_datalist
#'
#' i <- 9
#' brts <- datalist[[i]]$branching_times
#' trait <- 0
#' sampling_fraction <- c(1,1)
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
#' status = 2
#' DAISIE_DE_trait_logpES(
#'   brts                    = brts,
#'   trait                   = trait,
#'   status                  = status,
#'   sampling_fraction       = sampling_fraction,
#'   parameter               = parameter,
#'   trait_mainland_ancestor = NA,
#'   num_observed_states     = 2,
#'   num_hidden_states       = 2,
#'   atol                    = 1e-15,
#'   rtol                    = 1e-15,
#'   methode                 = "ode45",
#'   use_Rcpp                = 2)



DAISIE_DE_trait_logpES <- function(
    datalist,
    brts,
    parameter,
    trait,
    num_observed_states,
    num_hidden_states,
    trait_mainland_ancestor = NA, #this should contain either a full probability distribution across all states, only the observed states, or NA
    status,
    sampling_fraction,
    atol = 1e-15,
    rtol = 1e-15,
    methode = "ode45",
    rcpp_methode = "odeint::runge_kutta_cash_karp54",
    use_Rcpp = 2
) {

  lik_func <- function(i) {
    trait_mainland_ancestor_extended <- rep(0,num_observed_states * num_hidden_states)
    trait_mainland_ancestor_extended[i] <- 1 #set only the trait of interest to 1

    Lk_log <- DAISIE_DE_trait_logpES_core (brts                    = brts,
                                           parameter               = parameter,
                                           trait                   = trait,
                                           num_observed_states     = num_observed_states,
                                           num_hidden_states       = num_hidden_states,
                                           trait_mainland_ancestor = trait_mainland_ancestor_extended,
                                           status                  = status,
                                           sampling_fraction       = sampling_fraction,
                                           atol                    = atol,
                                           rtol                    = rtol,
                                           methode                 = "ode45",
                                           rcpp_methode            = rcpp_methode,
                                           use_Rcpp                = use_Rcpp)
    return(Lk_log)
  }

  indices <- seq_len(num_observed_states * num_hidden_states)
  Lk_vec <- sapply(indices, lik_func)

  ## added !all(is.na(trait_mainland_ancestor)) because when trait_mainland_ancestor = NA,  length(trait_mainland_ancestor) = length(trait_mainland_ancestor_extended) = 1
  if(!all(is.na(trait_mainland_ancestor)) && length(trait_mainland_ancestor) == num_observed_states * num_hidden_states) { #this is the case where a full probability distribution is specified across all observed and hidden states
    weights <- trait_mainland_ancestor / sum(trait_mainland_ancestor)
  } else {
    if(all(is.numeric(trait_mainland_ancestor))) { # this is the case when only a probability distribution is specified for the observed states; this could be c(M0/M, M1/M)

      ### the following calculates the terms before the + sign
      s <- numeric(num_observed_states * num_hidden_states)
      weights1 <- c()
      for(j in 1:length(trait_mainland_ancestor)) {
        s[((j - 1) * num_hidden_states + 1):(j * num_hidden_states)] <- rep(trait_mainland_ancestor[j], num_hidden_states)

        weights1 <- s
        weights1 <- s/sum(s)

      }

      ### the following calculates the terms after the + sign
      weights2 <- 1 - sum(trait_mainland_ancestor)
      weights <- weights1 + weights2

      if (all(weights == 0)) {
        weights <- weights
      } else {
        weights <- weights / sum(weights)
      }

    } else { # this is the case where nothing is provided, i.e. NA
      M0 <- datalist[[1]]$M0
      M1 <- datalist[[1]]$M1
      M <- datalist[[1]]$not_present
      trait_mainland_ancestor <- c(M0/M, M1/M)
      s <- numeric(num_observed_states * num_hidden_states)

      weights1 <- c()
      for(j in 1:length(trait_mainland_ancestor)) {
        s[((j - 1) * num_hidden_states + 1):(j * num_hidden_states)] <- rep(trait_mainland_ancestor[j], num_hidden_states)

        weights1 <- s
        weights1 <- s/sum(s)
      }

      ### the following calculates the terms after the + sign

      weights2 <- 1 - sum(trait_mainland_ancestor)
      weights <- weights1 + weights2

      if (all(weights == 0)) {
        weights <- weights
      } else {
        weights <- weights / sum(weights)
      }
    }
  }
  log_Lk <- log(sum(Lk_vec * weights))
  return( list (loglik = log_Lk, lik_states = Lk_vec, weights = weights))
}





DAISIE_DE_trait_logpES_core <- function(brts,
                                        status,
                                        trait,
                                        sampling_fraction,
                                        num_observed_states,
                                        num_hidden_states,
                                        trait_mainland_ancestor = NA,
                                        parameter,
                                        atol  = 1e-15,
                                        rtol  = 1e-15,
                                        methode                 = "ode45",
                                        rcpp_methode =
                                          "odeint::runge_kutta_cash_karp54",
                                        use_Rcpp = 0) {


  check_arguments(brts, parameter,
                  phy = 0,
                  trait,
                  num_observed_states,
                  num_hidden_states,
                  status,
                  sampling_fraction = sampling_fraction)





  # Unpack times from brts
  t0   <- brts[1]
  tmax <- brts[2]
  t1   <- brts[2]
  tp   <- 0

  # Time intervals

  time2 <- c(tp, t1)
  time3 <- c(tp, tmax)
  time4 <- c(tmax, t0)

  # Number of states in the system
  #n <- num_observed_states * num_hidden_states

  # Solve for interval [tp, t2] (stem phase)


  # Run appropriate sequence of intervals
  if ((status == 2 || status == 3) && length(brts) == 2) {
    initial_conditions2 <- get_initial_conditions2(status = status,
                                                   num_observed_states = num_observed_states,
                                                   num_hidden_states = num_hidden_states,
                                                   trait = trait,
                                                   brts = brts,
                                                   sampling_fraction = sampling_fraction,
                                                   trait_mainland_ancestor = trait_mainland_ancestor)

    solution2 <- solve_branch(interval_func = interval2,
                              initial_conditions = initial_conditions2,
                              time = time2,
                              parameter = parameter,
                              trait_mainland_ancestor = trait_mainland_ancestor,
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
                              trait_mainland_ancestor = trait_mainland_ancestor,
                              methode = methode,
                              rcpp_methode = rcpp_methode,
                              atol = atol,
                              rtol = rtol,
                              use_Rcpp = use_Rcpp)
  }

  if (status == 5) {
    initial_conditions3 <- get_initial_conditions3(status = status,
                                                   num_observed_states = num_observed_states,
                                                   num_hidden_states = num_hidden_states,
                                                   trait = trait,
                                                   sampling_fraction = sampling_fraction)
    solution3 <- solve_branch(interval_func = interval3,
                              initial_conditions = initial_conditions3,
                              time = time3,
                              parameter = parameter,
                              trait_mainland_ancestor = trait_mainland_ancestor,
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
                              trait_mainland_ancestor = trait_mainland_ancestor,
                              methode = methode,
                              rcpp_methode = rcpp_methode,
                              atol = atol,
                              rtol = rtol,
                              use_Rcpp = use_Rcpp)
  }

  # Extract log-likelihood from final solution
  Lk <- solution4[2, length(solution4[2, ])]

  return(Lk)
}
