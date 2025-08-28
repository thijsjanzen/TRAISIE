#' @keywords internal
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


