#' Simulate trait-dependent island replicates
#'
#' @description
#' Simulates island communities under a **trait-dependent** diversification
#' process with optional hidden states. This function is a convenience wrapper
#' that runs multiple replicates by calling
#' \code{DAISIE_sim_core_mult_trait_dep()} once per replicate. Time is held
#' constant over each run (no explicit time-varying rates), and baseline
#' geography/ontogeny settings are fixed inside the wrapper (see **Details**).
#'
#' If species share a single macro evolutionary process, supply the corresponding
#' entries in \code{trait_pars} for one observed state and set
#' \code{num_observed_states = 1}. If multiple observed states (and optionally
#' hidden states) are modeled, supply rate/transition parameters per state in
#' \code{trait_pars} and set \code{num_observed_states} and
#' \code{num_hidden_states} accordingly.
#'
#' @usage
#' DAISIE_sim_trait_dep(
#'   time,
#'   mainland,
#'   trait_pars,
#'   replicates,
#'   num_observed_states,
#'   num_hidden_states
#' )
#'
#' @param time Numeric scalar. Simulation horizon (e.g., island age in Myr).
#'   \code{time = 4} simulates the full 4 Myr history; \code{time = 2} stops at mid-life.
#' @param mainland List or numeric vector describing the mainland source pool
#'   (e.g., counts per trait state or clade). The core simulator expects a
#'   strictly positive total pool size.
#' @param trait_pars List of trait-dependent parameters consumed by
#'   \code{DAISIE_sim_core_mult_trait_dep()}, including per-state
#'   cladogenesis, extinction, colonization, anagenesis, and state-transition
#'   matrices (for observed and, if used, hidden states), and the parameter \code{p}
#' (\code{0} if a trait transition is not accompanied by anagenesis, \code{1} if it is).
#' @param replicates Integer. Number of independent island replicates to simulate.
#' @param num_observed_states Integer (>= 1). Number of **observed** trait states.
#' @param num_hidden_states Integer (>= 1). Number of **hidden** trait states; set
#'   to \code{1} if no hidden state is used.
#'
#' @details
#' The wrapper fixes several inputs to the core function to represent a
#' time-constant, no-ontogeny baseline:
#' \itemize{
#'   \item \code{island_ontogeny = 0} (no ontogeny),
#'   \item \code{sea_level = "const"} (constant sea level),
#'   \item \code{hyper_pars = create_hyper_pars(d = 0.027, x = 0.15)},
#'   \item \code{extcutoff = 1000},
#'   \item \code{area_pars = create_area_pars(max_area = 100, current_area = 90,
#'         proportional_peak_t = 0.5, total_island_age = 4,
#'         sea_level_amplitude = 5, sea_level_frequency = 10,
#'         island_gradient_angle = 0)}.
#' }
#'
#' @value
#' A list of length \code{replicates}. Each element is the return value from
#' \code{DAISIE_sim_core_mult_trait_dep()} for that replicate (or \code{NULL} if
#' the replicate failed). When at least one failure occurs, the returned list
#' carries an attribute \code{"failed"} with:
#' \itemize{
#'   \item \code{$indices}: integer vector of failed replicate indices,
#'   \item \code{$messages}: character vector of corresponding error messages.
#' }
#' A completion summary is emitted via \code{message()}.
#'
#' @section Failure handling:
#' Replicates are executed inside \code{tryCatch()}. On error, the replicate index
#' and message are appended to the \code{"failed"} attribute while other replicates
#' continue. Use \code{attr(x, "failed")} to inspect failures.
#'
#' @seealso
#' \code{\link{DAISIE_sim_core_mult_trait_dep}},
#' \code{\link{DAISIE_sim_cr}},
#' \code{\link{create_hyper_pars}},
#' \code{\link{create_area_pars}}
#'
#'
#' @examples
#' \dontrun{
#'trait_pars = list(immig_rate1 = 0.09,
#'                  ext_rate1 = 0.95,
#'                  ana_rate1 = 1.4,
#'                  clado_rate1 = 0.64,
#'                  immig_rate2 = 0.09,
#'                 ext_rate2 = 0.35,
#'                 ana_rate2 = 0.4,
#'                 clado_rate2 = 0.32,
#'                 trans_rate1 = 0.0,
#'                  trans_rate2 = 1.6,
#'                 trans_rate3 = 2.1,
#'                 trans_rate4 = 0.,
#'                  K1 = Inf,
#'                  K2 = Inf,
#'                  p = 0)
#'DAISIE_sim_trait_dep (  time = 20,
#'                         mainland = list(M1 = 100, M2 = 150),
#'                         trait_pars = trait_pars,
#'                         replicates = 1,
#'                         num_observed_states = 2,
#'                         num_hidden_states = 1)
#' }
#' @export


DAISIE_sim_trait_dep <- function(time,
                                 mainland,
                                 trait_pars,
                                 replicates,
                                 num_observed_states,
                                 num_hidden_states) {

  island_replicates <- vector("list", replicates)
  fail_idx  <- integer(0)
  fail_msg  <- character(0)

  for (rep in seq_len(replicates)) {
    cat(sprintf("replicate %d/%d ... ", rep, replicates)); flush.console()

    res <- tryCatch(
      {
        DAISIE_sim_core_mult_trait_dep(
          time = time,
          mainland = mainland,
          island_ontogeny = 0,
          sea_level = "const",
          hyper_pars = create_hyper_pars(d = 0.027, x = 0.15),
          extcutoff = 1000,
          area_pars = create_area_pars(
            max_area = 100,
            current_area = 90,
            proportional_peak_t = 0.5,
            total_island_age = 4,
            sea_level_amplitude = 5,
            sea_level_frequency = 10,
            island_gradient_angle = 0
          ),
          trait_pars = trait_pars,
          num_observed_states = num_observed_states,
          num_hidden_states = num_hidden_states
        )
      },
      error = function(e) {
        fail_idx <<- c(fail_idx, rep)
        fail_msg <<- c(fail_msg, conditionMessage(e))
        cat("failed\n")
        NULL
      }
    )

    if (!is.null(res)) cat("done\n")
    island_replicates[[rep]] <- res
  }

  if (length(fail_idx) > 0) {
    attr(island_replicates, "failed") <- list(indices = fail_idx, messages = fail_msg)
    message(sprintf("Completed with %d/%d failures: %s",
                    length(fail_idx), replicates, paste(fail_idx, collapse = ", ")))
  } else {
    message("All replicates completed successfully.")
  }

  island_replicates
}


?DAISIE_sim
