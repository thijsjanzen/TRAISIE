#' Core simulation for multi-trait, diversity-dependent DAISIE model
#'
#' Simulates island biodiversity dynamics under a clade-specific,
#' multi-trait, diversity-dependent constant-rate process.
#' This function is an extension of the DAISIE core simulation framework,
#' supporting multiple observed and hidden trait states.
#'
#' @param time Numeric scalar. Total simulation time (island age).
#' @param mainland Named list of mainland abundances (e.g., \code{list(M1 = 100, M2 = 150)}).
#' @param island_ontogeny Island ontogeny type (numeric or string). Defaults to \code{0}.
#' @param sea_level Sea level scenario (string). Defaults to \code{"const"}.
#' @param hyper_pars List of hyperparameters for the simulation model.
#' @param area_pars List of parameters controlling island area change (optional).
#' @param extcutoff Numeric. Extinction cutoff threshold. Defaults to \code{1000}.
#' @param trait_pars List of trait-dependent rates and parameters.
#' @param num_observed_states Integer. Number of observed trait states.
#' @param num_hidden_states Integer. Number of hidden trait states.
#'
#' @return A list representing the simulated island, including:
#'   \item{island_age}{Numeric island age.}
#'   \item{stt_table}{State-through-time table.}
#'   \item{island_spec}{Island species data.}
#'   and other DAISIE-like components.
#' @details
#' The simulation proceeds as a continuous-time Markov process,
#' with event rates (immigration, extinction, anagenesis, cladogenesis,
#' trait transitions) recalculated after each event. The function
#' supports multiple observed and hidden states, with final output
#' collapsing hidden states into observed states.
#'
#' @seealso \code{\link{DAISIE_create_island_trait}},
#'   \code{\link{update_rates_mult_trait}},
#'   \code{\link{DAISIE_sim_mult_trait_update_state_cr}}
#' @export

DAISIE_sim_core_mult_trait_dep <- function(
    time,
    mainland,
    island_ontogeny = 0,
    sea_level = "const",
    hyper_pars,
    area_pars,
    extcutoff = 1000,
    trait_pars,
    num_observed_states,
    num_hidden_states
) {

  #### Initialization ####
  timeval <- 0
  total_time <- time
  island_ontogeny <- DAISIE:::translate_island_ontogeny(island_ontogeny)
  sea_level <- DAISIE:::translate_sea_level(sea_level)


  testit::assert(length(trait_pars) > 5)



  mainland_total <- sum(unlist(mainland))

  testit::assert(mainland_total > 0)
  if(mainland[[1]] != 0){
    mainland_spec <- seq(1, mainland[[1]], 1)
  }else{
    mainland_spec <- c()
  }
  maxspecID <- mainland_total

  island_spec <- c()


  # if  (num_observed_states == 2)
  #  {
  #  stt_table <- matrix(ncol = 7)
  #  colnames(stt_table) <- c("Time","nI","nA","nC","nI2","nA2","nC2")
  #  stt_table[1,] <- c(total_time,0,0,0,0,0,0)

  #  } else if  (num_observed_states == 3) {
  #   stt_table <- matrix(ncol = 10)
  #   colnames(stt_table) <- c("Time","nI","nA","nC","nI2","nA2","nC2","nI3","nA3","nC3")
  #   stt_table[1,] <- c(total_time,0,0,0,0,0,0,0,0,0)
  #}


  stt_table <- matrix(ncol = 3 *(num_observed_states*num_hidden_states) + 1)  # 3 for each state (nI, nA, nC) and 1 for Time
  colnames(stt_table) <- c("Time",
                           paste0("nI", 1:(num_observed_states*num_hidden_states)),
                           paste0("nA", 1:(num_observed_states*num_hidden_states)),
                           paste0("nC", 1:(num_observed_states*num_hidden_states)))

  # Initialize the first row
  stt_table[1,] <- c(total_time, rep(0, 3 * (num_observed_states*num_hidden_states)))


  num_spec <- length(island_spec[, 1])


  #### Start Monte Carlo iterations ####
  while (timeval < total_time) {
    rates <- update_rates_mult_trait(timeval= timeval,
                                     total_time = total_time,
                                     hyper_pars = hyper_pars,
                                     area_pars = NULL,
                                     peak = NULL,
                                     island_ontogeny = NULL,
                                     sea_level = NULL,
                                     extcutoff = extcutoff,
                                     num_spec = num_spec,
                                     mainland = mainland,
                                     trait_pars = trait_pars,
                                     island_spec = island_spec,
                                     num_observed_states = num_observed_states,
                                     num_hidden_states = num_hidden_states)




    timeval_and_dt <- calc_next_timeval(
      max_rates = rates,
      timeval = timeval,
      total_time = total_time,
      num_observed_states = num_observed_states,
      num_hidden_states = num_hidden_states)
    timeval <- timeval_and_dt$timeval

    if (timeval < total_time) {

      possible_event <- DAISIE_sample_event_trait_dep(
        rates = rates,
        num_observed_states = num_observed_states,
        num_hidden_states = num_hidden_states
      )
      #print(possible_event)
      updated_state <- DAISIE_sim_mult_trait_update_state_cr(
        timeval = timeval,
        total_time = total_time,
        possible_event = possible_event,
        maxspecID = maxspecID,
        island_spec = island_spec,
        stt_table = stt_table,
        trait_pars = trait_pars,
        num_observed_states = num_observed_states,
        num_hidden_states = num_hidden_states,
        mainland = mainland
      )



      island_spec <- updated_state$island_spec
      maxspecID   <- updated_state$maxspecID
      stt_table   <- updated_state$stt_table
      num_spec    <- length(island_spec[, 1])

    }
  }

  ### change the true traits to the observed traits because the hidden states are unknown
  for (i in 1:length(island_spec[,1])) {

    state <- as.numeric(island_spec[i,][8])

    if (state >= 1 && state <= num_hidden_states) {
      island_spec[i,][8] = "0"

    } else if (state >= (num_hidden_states + 1) && state <= (2 * num_hidden_states)) {
      island_spec[i,][8] = "1"
      # Colonist species in state
    } else if (state >= (2 * num_hidden_states + 1) && state <= 3 * num_hidden_states) {
      island_spec[i,][8] = "2"
    } else if (state >= (3 * num_hidden_states + 1) && state <= 4 * num_hidden_states) {
      island_spec[i,][8] = "3"
    }
  }


  #### Finalize STT ####
  #stt_table <- rbind(
   # stt_table,
   # c(
    #  0,
    #  stt_table[nrow(stt_table), 2],
     # stt_table[nrow(stt_table), 3],
     # stt_table[nrow(stt_table), 4]
   # )
 # )
  island <- DAISIE_create_island_trait(
    stt_table = stt_table,
    total_time = total_time,
    island_spec = island_spec,
    mainland = mainland,
    trait_pars = trait_pars,
    num_observed_states = num_observed_states,
    num_hidden_states = num_hidden_states)
  # ordered_stt_times <- sort(island$stt_table[, 1], decreasing = TRUE)
  # testit::assert(all(ordered_stt_times == island$stt_table[, 1]))
  return(island)
}
