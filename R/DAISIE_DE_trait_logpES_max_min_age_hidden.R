#' testing function, for comparison with DAISIE
#' @description
#' this function calculates the likelihood of observing a singleton endemic species on an island
#' with the trait state `i`, and for which only the estimated maximum and minimum ages of colonization are known.
#' @export
#' @inheritParams default_params_doc
#' @examples
#' library(DAISIE)
#' data("Biwa_datalist")
#' datalist <- Biwa_datalist
#'
#'
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
#'     0,    0,    0.002,0.005,
#'     0,    .1000,    0,  0,
#'     0,    0,    0.100,0.00
#'   ), nrow = 4),
#'   1
#' )
#'
#'
#'
#'
#' DAISIE_DE_trait_logpES_max_min_age_hidden(
#'   brts                  = data_list1[[2]]$branching_times,
#'   trait                 = 0,
#'   status                = 9,
#'   parameter             = parameter,
#'   num_observed_states   = 2,
#'   num_hidden_states     = 2,
#'   atol                  = 1e-15,
#'   rtol                  = 1e-15,
#'   methode               = "ode45",
#'   trait_mainland_ancestor = c(0, 1),
#'   sampling_fraction     = c(1,1),
#'   use_Rcpp = 2
#' )



DAISIE_DE_trait_logpES_max_min_age_hidden <- function(
    brts,
    parameter,
    trait,
    num_observed_states,
    num_hidden_states,
    trait_mainland_ancestor, #this should contain either a full probability distribution across all states, only the observed states, or NA
    status,
    sampling_fraction,
    Mainland_pool_size_vec = NULL,
    atol = 1e-15,
    rtol = 1e-15,
    methode = "ode45",
    rcpp_methode = "odeint::runge_kutta_cash_karp54",
    use_Rcpp = 2
) {

  Lk_vec <- numeric(num_observed_states * num_hidden_states)

  for (i in seq_len(num_observed_states * num_hidden_states)) { #loop over all possible states, observed and hidden, one by one
    trait_mainland_ancestor_extended <- rep(0,num_observed_states * num_hidden_states)
    trait_mainland_ancestor_extended[i] <- 1 #set only the trait of interest to 1

    Lk_log <- DAISIE_DE_trait_logpES_max_min_age_hidden_core (brts,
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
    Lk_vec[i] <- Lk_log # ideally this should not be needed if the function above does not do logtransformation
  }

  ## added !all(is.na(trait_mainland_ancestor)) because when trait_mainland_ancestor = NA,  length(trait_mainland_ancestor) = length(trait_mainland_ancestor_extended) = 1
  if(!all(is.na(trait_mainland_ancestor)) && length(trait_mainland_ancestor) == length(trait_mainland_ancestor_extended)) { #this is the case where a full probability distribution is specified across all observed and hidden states
    weights <- trait_mainland_ancestor/sum(trait_mainland_ancestor)
  } else {
    if(all(is.numeric(trait_mainland_ancestor))) { # this is the case when only a probability distribution is specified for the observed states; this could be c(M0/M, M1/M)

      s <- numeric(num_observed_states * num_hidden_states)

      weights1 <- c()
      for(j in 1:length(trait_mainland_ancestor)) {
        s[((j - 1) * num_hidden_states + 1):(j * num_hidden_states)] <- rep(trait_mainland_ancestor[j], num_hidden_states)

        weights_j <- Lk_vec[((j - 1) * num_hidden_states + 1):(j * num_hidden_states)]

        if (sum(weights_j) == 0)
        {
          weights_j <- weights_j/1
        }else{
          weights_j <- weights_j/sum(weights_j)
        }
        weights1 <- c(weights1, weights_j)
      }
      weights1 <- weights1 * s/sum(weights1)


      weights2 <- Lk_vec * (1 - sum(trait_mainland_ancestor))/sum(Lk_vec)

      weights <- weights1 + weights2

      if (all(weights == 0)) {
        weights <- weights
      } else {
        weights <- weights / sum(weights)
      }

    } else { # this is the case where nothing is provided, i.e. NA
      weights <- Lk_vec/sum(Lk_vec)
    }
  }
  log_Lk <- log(sum(Lk_vec * weights))
  return( list (loglik = log_Lk, lik_states = Lk_vec, weights = weights))
}








DAISIE_DE_trait_logpES_max_min_age_hidden_core <- function(brts,
                                                      trait,
                                                      status,
                                                      sampling_fraction = 1,
                                                      parameter,
                                                      trait_mainland_ancestor = NA,
                                                      num_observed_states,
                                                      num_hidden_states,
                                                      atol = 1e-15,
                                                      rtol = 1e-15,
                                                      methode = "ode45",
                                                      rcpp_methode = "odeint::runge_kutta_cash_karp54",
                                                      use_Rcpp = 0) {
  t0   <- brts[1]
  tmax <- brts[2]
  tmin <- brts[3]
  tp   <- 0

  # number of unique state
  n <- num_observed_states * num_hidden_states

  #########interval2 [t_p, tmin]

  m = length(parameter[[1]])


  ## SOLVED: can't we call 'get_initial_conditions' here? //NO, because brts > 2
  initial_conditions2 <- get_initial_conditions2(status = status,
                                                 num_observed_states = num_observed_states,
                                                 num_hidden_states = num_hidden_states,
                                                 trait = trait,
                                                 brts = brts,
                                                 sampling_fraction = sampling_fraction,
                                                 trait_mainland_ancestor = trait_mainland_ancestor)

  # Time sequence for interval [tp, tmin]
  time2 <- c(tp, tmin)

  solution2 <- solve_branch(interval_func = interval2,
                            initial_conditions = initial_conditions2,
                            time = time2,
                            parameter = parameter,
                            methode = methode,
                            rcpp_methode = rcpp_methode,
                            trait_mainland_ancestor = trait_mainland_ancestor,
                            atol = atol,
                            rtol = rtol,
                            use_Rcpp = use_Rcpp)

  #########interval3 [tmin, tmax]

  # Initial conditions

  # only use second row, because the first row of solution3 is the initial state
  initial_conditions3_max_min <- c(solution2[2,][1:n],                                             ### DE: select DE in solution2
                                   rep(0, n),                                                      ### DM1: select DE in solution2
                                   solution2[2,][(n + 1):(n + n)],                         ### DM2: select DM2 in solution2
                                   solution2[2,][(n + n + 1):(n + n + n)],                 ### DM3: select DM3 in solution2
                                   solution2[2,][(n + n + n + 1):(n + n + n + n)],         ### E: select E in solution2
                                   0,                                                              ### DA2
                                   solution2[2,][length(solution2[2,])])                           ### DA3: select DA3 in solution2

  initial_conditions3_max_min <- matrix(initial_conditions3_max_min, nrow = 1)

  # Time sequence for interval [tmin, tmax]
  time3 <- c(tmin, tmax)

  solution3 <- solve_branch(interval_func = interval3,
                            initial_conditions = initial_conditions3_max_min,
                            time = time3,
                            parameter = parameter,
                            trait_mainland_ancestor = trait_mainland_ancestor,
                            methode = methode,
                            rcpp_methode = rcpp_methode,
                            atol = atol,
                            rtol = rtol,
                            use_Rcpp = use_Rcpp)

  #########interval4 [tmax, t0]

  # Initial conditions

  # only use second row, because the first row of solution3 is the initial state
  initial_conditions4_max_min <- c(solution3[2,][(n + n + 1):(n + n + n)],                         ### DM1: select DM2 in solution3
                                   solution3[2,][(n + n + n + n + 1):(n + n + n + n + n)],         ### E: select E in solution3
                                   solution3[2,][length(solution3[2,]) - 1])

  initial_conditions4_max_min <- matrix(initial_conditions4_max_min, nrow = 1)

  # Time sequence for interval [tmax, t0]
  time4 <- c(tmax, t0)

  # Solve the system for interval [tmax, t0]
  solution4 <- solve_branch(interval_func = interval4,
                            initial_conditions = initial_conditions4_max_min,
                            time = time4,
                            parameter = parameter,
                            trait_mainland_ancestor = trait_mainland_ancestor,
                            methode = methode,
                            rcpp_methode = rcpp_methode,
                            atol = atol,
                            rtol = rtol,
                            use_Rcpp = use_Rcpp)

  # Extract log-likelihood
  Lk <- solution4[2,][length(solution4[2,])]

  return(Lk)
}
