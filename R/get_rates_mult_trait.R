#' @keywords internal

############ hidden

get_immig_rate <- function(A = 1,
                           num_spec,
                           mainland,
                           trait_pars,
                           island_spec = NULL,
                           num_observed_states,
                           num_hidden_states) {
  immig_list <- list()

  for (i in 1:(num_observed_states*num_hidden_states)) {


    if (i >= 1 && i <= num_hidden_states) {
      mainland_n <- mainland[[1]]
    } else if (i > num_hidden_states && i <= (2 * num_hidden_states)) {
      mainland_n <- mainland[[2]]
    } else if (i > (2 * num_hidden_states) && i <= (3 * num_hidden_states)) {
      mainland_n <- mainland[[3]]
    } else if (i > (3 * num_hidden_states) && i <= (4 * num_hidden_states)) {
      mainland_n <- mainland[[4]]
    }


    gam <- trait_pars[[paste("immig_rate", i, sep = "")]]
    num_spec_trait <- length(which(island_spec[, 8] == as.character(i)))

    if (paste("K", i, sep = "") %in% names(trait_pars)) {
      immig_rate <- max(c(mainland_n * gam * (1 - (num_spec_trait / trait_pars[[paste("K", i, sep = "")]])),
                          0), na.rm = TRUE)
    } else {

      ### K is the same for all type of lineage
      immig_rate <- max(c(mainland_n * gam * (1 - (num_spec_trait / trait_pars$K1)),
                          0), na.rm = TRUE)
    }

    immig_list[[paste("immig_rate", i, sep = "")]] <- immig_rate
  }

  return(immig_list)
}

##################

get_ext_rate <- function(hyper_pars,
                         extcutoff = 1000,
                         num_spec,
                         A = 1,
                         trait_pars,
                         island_spec = NULL,
                         num_observed_states,
                         num_hidden_states) {
  ext_list <- list()

  for (i in 1:(num_observed_states*num_hidden_states)) {

    num_spec_trait <- length(which(island_spec[, 8] == as.character(i)))
    ext_rate <- trait_pars[[paste("ext_rate", i, sep = "")]] * num_spec_trait
    ext_list[[paste("ext_rate", i, sep = "")]] <- ext_rate

  }

  return(ext_list)
}

################

get_ana_rate <- function(island_spec = NULL,
                         trait_pars,
                         num_observed_states,
                         num_hidden_states) {
  ana_list <- list()

  for (i in 1:(num_observed_states*num_hidden_states)) {

    ana_rate <- trait_pars[[paste("ana_rate", i, sep = "")]] * length(
      intersect(which(island_spec[, 4] == "I"),
                which(island_spec[, 8] == as.character(i)))
    )
    ana_list[[paste("ana_rate", i, sep = "")]] <- ana_rate
  }

  return(ana_list)
}

##############

get_clado_rate <- function(hyper_pars,
                           num_spec,
                           A,
                           trait_pars,
                           island_spec = NULL,
                           num_observed_states,
                           num_hidden_states) {
  clado_list <- list()

  for (i in 1:(num_observed_states*num_hidden_states)) {
    num_spec_trait <- length(which(island_spec[, 8] == as.character(i)))

    if (paste("K", i, sep = "") %in% names(trait_pars)) {
      clado_rate <- max(
        0, trait_pars[[paste("clado_rate", i, sep = "")]] * num_spec_trait * (1 - num_spec_trait / trait_pars[[paste("K", i, sep = "")]]),
        na.rm = TRUE)

    } else {
      clado_rate <- max(
        0, trait_pars[[paste("clado_rate", i, sep = "")]] * num_spec_trait * (1 - num_spec / trait_pars$K1),
        na.rm = TRUE)
    }

    clado_list[[paste("clado_rate", i, sep = "")]] <- clado_rate
  }

  return(clado_list)
}

###############
get_trans_rate <- function(trait_pars,
                           island_spec,
                           num_observed_states,
                           num_hidden_states, p) {
  trans_list <- list()

  for (i in 1:((num_observed_states*num_hidden_states)^2)) {

    n <- (num_observed_states*num_hidden_states)


    if (i %in% 1:n) {
      num_spec_trait <- length(which(island_spec[, 8] == as.character(1)))
    } else if  (i %in% (1 + n):( 2* n)) {
      num_spec_trait <- length(which(island_spec[, 8] == as.character(2)))
    }  else if (i %in% (1 + 2*n):( 3* n)) {
      num_spec_trait <- length(which(island_spec[, 8] == as.character(3)))
    } else if  (i %in% (1 + 3*n):( 4* n)) {
      num_spec_trait <- length(which(island_spec[, 8] == as.character(4)))
    }


    #from_state <- ((i - 1) %/% num_hidden_states) + 1
    # num_spec_trait <- sum(island_spec[,8] == as.character(from_state))

    trans_rate <- trait_pars[[paste("trans_rate", i, sep = "")]] * num_spec_trait
    trans_list[[paste("trans_rate", i, sep = "")]] <- trans_rate
  }

  return(trans_list)
}






