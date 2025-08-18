#' @rawNamespace useDynLib(treeLL, .registration = TRUE)
#' @rawNamespace import(Rcpp)
#' @rawNamespace importFrom(RcppParallel, RcppParallelLibs)
#' @keywords internal
transf_funcdefpar <- function(idparsfuncdefpar,
                              functions_defining_params,
                              idfactorsopt,
                              trparsfix,
                              trparsopt,
                              idparsfix,
                              idparsopt) {
  trparfuncdefpar <- NULL
  ids_all <- c(idparsfix, idparsopt)

  values_all <- c(trparsfix / (1 - trparsfix),
                  trparsopt / (1 - trparsopt))
  a_new_envir <- new.env()
  x <- as.list(values_all)  ## To declare all the ids as variables

  if (is.null(idfactorsopt)) {
    names(x) <- paste0("par_", ids_all)
  } else {
    names(x) <- c(paste0("par_", ids_all), paste0("factor_", idfactorsopt))
  }
  list2env(x, envir = a_new_envir)

  for (jj in seq_along(functions_defining_params)) {
    myfunc <- functions_defining_params[[jj]]
    environment(myfunc) <- a_new_envir
    value_func_defining_parm <- local(myfunc(), envir = a_new_envir)

    ## Now, declare the variable that is just calculated, so it is available
    ## for the next calculation if needed
    y <- as.list(value_func_defining_parm)
    names(y) <- paste0("par_", idparsfuncdefpar[jj])
    list2env(y, envir = a_new_envir)

    if (is.numeric(value_func_defining_parm) == FALSE) {
      stop("Something went wrong with the calculation of
                 parameters in 'functions_param_struct'")
    }
    trparfuncdefpar <- c(trparfuncdefpar, value_func_defining_parm)
  }
  trparfuncdefpar <- trparfuncdefpar / (1 + trparfuncdefpar)
  rm(a_new_envir)
  return(trparfuncdefpar)
}

#' @keywords internal
update_values_transform_cla <- function(trpars,
                                        idparslist,
                                        idpars,
                                        parvals) {
  for (i in seq_along(idpars)) {
    for (j in seq_len(nrow(trpars[[3]]))) {
      id <- which(idparslist[[1]][[j]] == idpars[i])
      trpars[[1]][[j]][id] <- parvals[i]
    }
    for (j in 2:3) {
      id <- which(idparslist[[j]] == idpars[i])
      trpars[[j]][id] <- parvals[i]
    }
  }
  return(trpars)
}

#' @keywords internal
transform_params_cla <- function(idparslist,
                                 idparsfix,
                                 trparsfix,
                                 idparsopt,
                                 trparsopt,
                                 structure_func,
                                 idparsfuncdefpar,
                                 trparfuncdefpar) {
  trpars1 <- idparslist
  for (j in seq_len(nrow(trpars1[[3]]))) {
    trpars1[[1]][[j]][, ] <- NA
  }

  for (j in 2:3) {
    trpars1[[j]][] <- NA
  }

  if (length(idparsfix) != 0) {
    trpars1 <- update_values_transform_cla(trpars1,
                                           idparslist,
                                           idparsfix,
                                           trparsfix)
  }

  trpars1 <- update_values_transform_cla(trpars1,
                                         idparslist,
                                         idparsopt,
                                         trparsopt)
  ## structure_func part
  if (!is.null(structure_func)) {
    trpars1 <- update_values_transform_cla(trpars1,
                                           idparslist,
                                           idparsfuncdefpar,
                                           trparfuncdefpar)
  }

  pre_pars1 <- list()
  pars1 <- list()

  for (j in seq_len(nrow(trpars1[[3]]))) {
    pre_pars1[[j]] <- trpars1[[1]][[j]][, ] / (1 - trpars1[[1]][[j]][, ])
  }

  pars1[[1]] <- pre_pars1
  for (j in 2:3) {
    pars1[[j]] <- trpars1[[j]] / (1 - trpars1[[j]])
  }

  return(pars1)
}

#' @keywords internal
update_values_transform <- function(trpars,
                                    idparslist,
                                    idpars,
                                    parvals) {
  for (i in seq_along(idpars)) {
    for (j in 1:6) {
      id <- which(idparslist[[j]] == idpars[i])
      trpars[[j]][id] <- parvals[i]
    }
  }
  return(trpars)
}

#' @keywords internal
transform_params_normal <- function(idparslist,
                                    idparsfix,
                                    trparsfix,
                                    idparsopt,
                                    trparsopt,
                                    structure_func,
                                    idparsfuncdefpar,
                                    trparfuncdefpar) {
  trpars1 <- idparslist
  for (j in 1:6) {
    trpars1[[j]][] <- NA
  }
  if (length(idparsfix) != 0) {
    trpars1 <- update_values_transform(trpars1,
                                       idparslist,
                                       idparsfix,
                                       trparsfix)
  }

  trpars1 <- update_values_transform(trpars1,
                                     idparslist,
                                     idparsopt,
                                     trparsopt)

  ## if structure_func part
  if (is.null(structure_func) == FALSE) {
    trpars1 <- update_values_transform(trpars1,
                                       idparslist,
                                       idparsfuncdefpar,
                                       trparfuncdefpar)
  }

  pars1 <- list()
  for (j in 1:6) {
    pars1[[j]] <- trpars1[[j]] / (1 - trpars1[[j]])
  }
  return(pars1)
}

#' @keywords internal
transform_parameters <- function(trparsopt,
                                        trparsfix,
                                        idparsopt,
                                        idparsfix,
                                        idparslist,
                                        structure_func = NULL) {
  if (!is.null(structure_func)) {
    idparsfuncdefpar <- structure_func[[1]]
    functions_defining_params <- structure_func[[2]]

    if (length(structure_func[[3]]) > 1) {
      idfactorsopt <- structure_func[[3]]
    } else {
      if (structure_func[[3]] == "noFactor") {
        idfactorsopt <- NULL
      } else {
        idfactorsopt <- structure_func[[3]]
      }
    }

    trparfuncdefpar <- transf_funcdefpar(idparsfuncdefpar =
                                           idparsfuncdefpar,
                                         functions_defining_params =
                                           functions_defining_params,
                                         idfactorsopt = idfactorsopt,
                                         trparsfix = trparsfix,
                                         trparsopt = trparsopt,
                                         idparsfix = idparsfix,
                                         idparsopt = idparsopt)
  }

  if (is.list(idparslist[[1]])) {
    # when the ml function is called from cla_secsse
    pars1 <- transform_params_cla(idparslist,
                                  idparsfix,
                                  trparsfix,
                                  idparsopt,
                                  trparsopt,
                                  structure_func,
                                  idparsfuncdefpar,
                                  trparfuncdefpar)
  } else {
    # when non-cla option is called
    pars1 <- transform_params_normal(idparslist,
                                     idparsfix,
                                     trparsfix,
                                     idparsopt,
                                     trparsopt,
                                     structure_func,
                                     idparsfuncdefpar,
                                     trparfuncdefpar)
  }
  return(pars1)
}


#' Times at which speciation or extinction occurs
#' @title Event times of a (possibly non-ultrametric) phylogenetic tree
#' @param phy phylogenetic tree of class phylo, without polytomies, rooted and
#' with branch lengths. Need not be ultrametric.
#' @return times at which speciation or extinction happens.
#' @note This script has been modified from BAMMtools' internal function
#' NU.branching.times
#' @export
event_times <- function(phy) {
  if (ape::is.ultrametric(phy)) {
    return(ape::branching.times(phy))
  } else {
    if (ape::is.binary(phy) == FALSE) {
      stop("error. Need fully bifurcating (resolved) tree\n")
    }
    phy$begin <- rep(0, nrow(phy$edge))
    phy$end <- rep(0, nrow(phy$edge))
    fx <- function(phy, node) {
      cur_time <- 0
      root <- length(phy$tip.label) + 1
      if (node > root) {
        cur_time <- phy$end[which(phy$edge[, 2] == node)]
      }
      dset <- phy$edge[, 2][phy$edge[, 1] == node]
      i1 <- which(phy$edge[, 2] == dset[1])
      i2 <- which(phy$edge[, 2] == dset[2])
      phy$end[i1] <- cur_time + phy$edge.length[i1]
      phy$end[i2] <- cur_time + phy$edge.length[i2]
      if (dset[1] > length(phy$tip.label)) {
        phy$begin[phy$edge[, 1] == dset[1]] <- phy$end[i1]
        phy <- fx(phy, node = dset[1])
      }
      if (dset[2] > length(phy$tip.label)) {
        phy$begin[phy$edge[, 1] == dset[2]] <- phy$end[i2]
        phy <- fx(phy, node = dset[2])
      }
      return(phy)
    }
    phy <- fx(phy, node = length(phy$tip.label) + 1)
    maxbt <- max(phy$end)
    nodes <- (length(phy$tip.label) + 1):(2 * length(phy$tip.label) - 1)
    bt <- numeric(length(nodes))
    names(bt) <- nodes
    for (i in seq_along(bt)) {
      tt <- phy$begin[phy$edge[, 1] == nodes[i]][1]
      bt[i] <- maxbt - tt
    }
    return(bt)
  }
}

#' Print likelihood for initial parameters
#'
#' @inheritParams default_params_doc
#'
#' @return Invisible `NULL`. Prints a `message()` to the console with the
#'   initial loglikelihood if `verbose >= 1`
#' @noRd
print_init_ll <- function(initloglik) {
  init_ll_msg1 <- "Calculating the likelihood for the initial parameters."
  init_ll_msg2 <-
    paste0("The loglikelihood for the initial parameter values is ",
           initloglik)
  init_ll_msg3 <- c("Optimizing the likelihood - this may take a while.")
  message(paste(init_ll_msg1, init_ll_msg2, init_ll_msg3, sep = "\n"))

  invisible(NULL)
}

#' @keywords internal
check_arguments <- function(brts = NULL,
                            parameter = NULL,
                            phy = NULL,
                            traits = NULL,
                            num_observed_states = NULL,
                            num_hidden_states = NULL,
                            status = NULL,
                            sampling_fraction = NULL) {

  if (is.null(brts)) {
    stop("brts not provided")
  }
  if (is.null(parameter)) {
    stop("parameters not provided at all")
  }
  if (!is.list(parameter)) {
    stop("parameters need to be provided as a list")
  }
  if (is.null(phy)) {
    stop("phy not provided")
  }
  if (is.null(num_observed_states)) {
    stop("number of observed states not provided")
  }
  if (is.null(num_hidden_states)) {
    stop("number of hidden states not provided")
  }
  if (is.null(status)) {
    stop("status not provided")
  }
  if (is.null(sampling_fraction)) {
    stop("sampling_fraction not provided")
  }
  # Only require traits if there is more than one state in total
  total_states <- as.numeric(num_observed_states) * as.numeric(num_hidden_states)
  if (is.na(total_states)) stop("num_observed_states and num_hidden_states must be numeric")
  if (total_states > 1 && is.null(traits)) {
    stop("traits not provided (required when total states > 1)")
  }
}
