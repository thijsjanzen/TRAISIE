#' Calculate Next Event Time in Trait-Dependent DAISIE Model (with Hidden States)
#'
#' This function computes the waiting time to the next evolutionary event (immigration,
#' extinction, anagenesis, cladogenesis, or trait transition) in a trait-dependent DAISIE
#' simulation that includes both observed and hidden states.
#'
#' @param max_rates A named list of rates for all possible events (same as used in
#'        \code{DAISIE_sample_event_trait_dep}). It must include at least the first
#'        \code{4 * (num_observed_states * num_hidden_states) +
#'        (num_observed_states * num_hidden_states)^2} entries.
#' @param timeval Current simulation time.
#' @param total_time Total duration of the simulation (end time).
#' @param num_observed_states Integer. Number of observed trait states.
#' @param num_hidden_states Integer. Number of hidden trait states per observed state.
#'
#' @return A named list with:
#' \describe{
#'   \item{timeval}{Updated simulation time after sampling a waiting time}
#'   \item{dt}{Sampled waiting time (time increment)}
#' }
#'
#' @details
#' The total event rate is calculated as the sum of the first
#' \code{4 * (num_observed_states * num_hidden_states) +
#' (num_observed_states * num_hidden_states)^2} entries of the \code{max_rates} list,
#' which corresponds to immigration, extinction, anagenesis, cladogenesis, and all possible
#' trait transitions across combined states.
#'
#' The next waiting time is drawn from an exponential distribution with rate equal to
#' the total event rate. If the total rate is 0 (i.e., no possible events), the simulation
#' time jumps directly to \code{total_time}.
#'
#' @examples
#' max_rates <- as.list(runif(4 * 4 + 16, min = 0.01, max = 0.1))  # 4 combined states
#' t0 <- 0
#' total_time <- 10
#' res <- calc_next_timeval(max_rates, t0, total_time, num_observed_states = 2, num_hidden_states = 2)
#' print(res)
#'


##########

calc_next_timeval <- function(max_rates, timeval, total_time, num_observed_states, num_hidden_states) {


  # In the vector ot rates, exclude all keys that start with the Ks and p
  # or sum only the possible events

  number_of_possible_events <- 4*(num_observed_states * num_hidden_states) + ((num_observed_states * num_hidden_states)^2)

  totalrate <- sum(unlist(max_rates[1:number_of_possible_events]))


  # Calculate the next time value
  if (totalrate != 0) {
    dt <- stats::rexp(1, totalrate)
    timeval <- timeval + dt
  } else {
    timeval <- total_time
  }

  return(list(timeval = timeval, dt = dt))
}

