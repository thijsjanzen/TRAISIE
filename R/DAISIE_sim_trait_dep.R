#' Simulates island replicates with an clade-specific (CS) diversity-dependent
#' constant-rate process
#'
#' @inheritParams default_params_doc
#'
#' @return A list. The highest level of the least corresponds to each individual
#' replicate. See return for `DAISIE_sim_cr()` for details.
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
