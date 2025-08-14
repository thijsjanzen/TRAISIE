

#########with hidden
update_rates_mult_trait <- function(timeval,
                                    total_time,
                                    hyper_pars = hyper_pars,
                                    area_pars = NULL,
                                    peak = NULL,
                                    island_ontogeny = NULL,
                                    sea_level = NULL,
                                    extcutoff,
                                    num_spec,
                                    mainland,
                                    trait_pars,
                                    island_spec = NULL,
                                    num_observed_states,
                                    num_hidden_states) {

  # Function to calculate rates at time = timeval. Returns a list with each rate.
  island_ontogeny <- if (!is.null(area_pars$island_ontogeny)) area_pars$island_ontogeny else 0
  sea_level <- if (!is.null(area_pars$sea_level)) area_pars$sea_level else 0

  # Calculate the island area at the given time
  A <- DAISIE:::island_area(
    timeval = timeval,
    total_time = total_time,
    area_pars = area_pars,
    peak = peak,
    island_ontogeny = island_ontogeny,
    sea_level = sea_level
  )

  # Get immigration rate
  immig_rate <- get_immig_rate(
    A = A,
    num_spec = num_spec,
    mainland = mainland,
    trait_pars = trait_pars,
    island_spec = island_spec,
    num_observed_states = num_observed_states,
    num_hidden_states = num_hidden_states
  )

  # Get extinction rate
  ext_rate <- get_ext_rate(
    hyper_pars = hyper_pars,
    extcutoff = extcutoff,
    num_spec = num_spec,
    A = A,
    trait_pars = trait_pars,
    island_spec = island_spec,
    num_observed_states = num_observed_states,
    num_hidden_states = num_hidden_states
  )

  # Get anagenesis rate
  ana_rate <- get_ana_rate(
    trait_pars = trait_pars,
    island_spec = island_spec,
    num_observed_states = num_observed_states,
    num_hidden_states = num_hidden_states
  )

  # Get cladogenesis rate
  clado_rate <- get_clado_rate(
    hyper_pars = hyper_pars,
    num_spec = num_spec,
    A = A,
    trait_pars = trait_pars,
    island_spec = island_spec,
    num_observed_states = num_observed_states,
    num_hidden_states = num_hidden_states
  )

  # Get transition rate
  trans_rate <- get_trans_rate(
    trait_pars = trait_pars,
    island_spec = island_spec,
    num_observed_states = num_observed_states,
    num_hidden_states = num_hidden_states
  )

  # Initialize rates list
  rates <- list()

  # Dynamically generate the rates for each state
  for (i in 1:(num_observed_states*num_hidden_states)) {
    rates[[paste("immig_rate", i, sep = "")]] <- immig_rate[[paste("immig_rate", i, sep = "")]]
    rates[[paste("ext_rate", i, sep = "")]] <- ext_rate[[paste("ext_rate", i, sep = "")]]
    rates[[paste("ana_rate", i, sep = "")]] <- ana_rate[[paste("ana_rate", i, sep = "")]]
    rates[[paste("clado_rate", i, sep = "")]] <- clado_rate[[paste("clado_rate", i, sep = "")]]
  }

  # Dynamically add transition rates
  for (i in 1:((num_observed_states*num_hidden_states)^2)) {

    rates[[paste("trans_rate", i, sep = "")]] <- trans_rate[[paste("trans_rate", i, sep = "")]]

  }

  return(rates)
}

