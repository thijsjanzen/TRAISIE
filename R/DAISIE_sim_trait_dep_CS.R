
DAISIE_sim_trait_dep_CS <- function (total_time, mainland, trait_pars, replicates, nonoceanic_pars, prop_type2_pool,
                              replicates_apply_type2, sample_freq, hyper_pars, area_pars,
                              cond, verbose, files_to_write = 0, num_observed_states, num_hidden_states)
{
  island_replicates <- list()

  for (rep in 1:replicates) {
    island_replicates[[rep]] <- list()
    full_list <- list()
    if (cond == 0) {
      number_present <- -1
    }
    else {
      number_present <- 0
    }

    mainland_total <- sum(unlist(mainland))
    while (number_present < cond) {
      for (m_spec in 1: mainland_total) {
        full_list[[m_spec]] <- DAISIE_sim_core_mult_trait_dep(
          time = total_time,
          mainland = 1 ,
          sea_level = "const",
          extcutoff = 300,
          area_pars = area_pars,
          hyper_pars = hyper_pars,
          trait_pars = trait_pars,
          num_observed_states = num_observed_states,
          num_hidden_states = num_hidden_states
        )
      }


      stac_vec <- unlist(full_list)[which(names(unlist(full_list)) ==
                                            "stac")]
      present <- which(stac_vec != 0)
      number_present <- length(present)
    }
    island_replicates[[rep]] <- full_list
    if (verbose == TRUE) {
      message("Island replicate ", rep)
    }
  }


  if (files_to_write > 0) {
    for (filenum in 1:files_to_write) {
      chunks <- ceiling(seq_along(1:replicates)/files_to_write)
      start <- min(which(chunks == filenum))
      end <- max(which(chunks == filenum))
      island_reps <- island_replicates[start:end]
      save(start, end, island_reps, file = paste0("DAISIE_sims",
                                                  start, "-", end, ".Rdata"))
    }
  }



  if (files_to_write == 0) {
    island_replicates <- DAISIE:::DAISIE_format_CS(island_replicates = island_replicates,
                                                   time = total_time, M = mainland[[1]], sample_freq = sample_freq,
                                                   verbose = verbose)
  }


  if (files_to_write > 0) {
    rm(island_replicates)
    for (filenum in 1:files_to_write) {
      chunks <- ceiling(seq_along(1:replicates)/files_to_write)
      start <- min(which(chunks == filenum))
      end <- max(which(chunks == filenum))
      load(paste0("DAISIE_sims", start, "-", end, ".Rdata"))
      island_replicates <- DAISIE:::DAISIE_format_CS(island_replicates = island_reps,
                                                     time = total_time, M = mainland[[1]], sample_freq = sample_freq,
                                                     verbose = verbose)
      save(start, end, island_replicates, file = paste0("DAISIE_sims_formatted",
                                                        start, "-", end, ".Rdata"))
    }
  }
  return(island_replicates)
}


