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

  island_replicates <- list()

  for (rep in 1:replicates) {
    island_replicates[[rep]] <- DAISIE_sim_core_mult_trait_dep(
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
  }

  return(island_replicates)
}

