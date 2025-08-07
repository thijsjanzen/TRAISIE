#' apply maximum likelihood to a dataset
#' @description ml function
#' @inheritParams default_params_doc
#' @export
calc_ml <- function(datalist,
                      num_observed_states,
                      num_hidden_states,
                      idparslist,
                      idparsopt,
                      initparsopt,
                      idparsfix,
                      parsfix,
                      cond = "proper_cond",
                      tol = c(1e-04, 1e-05, 1e-07),
                      maxiter = 1000 * round((1.25) ^ length(idparsopt)),
                      optimmethod = "simplex",
                      methode = "ode45",
                      rcpp_methode = "odeint::bulirsch_stoer",
                      num_cycles = 1,
                      verbose = FALSE,
                      num_threads = 1,
                      atol = 1e-15,
                      rtol = 1e-15,
                      use_Rcpp = 0
                      ) {
  if (identical(as.numeric(sort(c(idparsopt, idparsfix))),
                as.numeric(sort(unique(unlist(idparslist))))) == FALSE) {
    stop("All elements in idparslist must be included in either
             idparsopt or idparsfix ")
  }

  see_ancestral_states <- FALSE

  trparsopt <- initparsopt / (1 + initparsopt)
  trparsopt[which(initparsopt == Inf)] <- 1
  trparsfix <- parsfix / (1 + parsfix)
  trparsfix[which(parsfix == Inf)] <- 1

  optimpars <- c(tol, maxiter, verbose)

  initloglik <- loglik_choosepar(trparsopt = trparsopt,
                                 trparsfix = trparsfix,
                                 idparsopt = idparsopt,
                                 idparsfix = idparsfix,
                                 idparslist = idparslist,
                                 datalist = datalist,
                                 num_observed_states = num_observed_states,
                                 num_hidden_states = num_hidden_states,
                                 cond = cond,
                                 atol = atol,
                                 rtol = rtol,
                                 methode = methode,
                                 rcpp_methode = rcpp_methode,
                                 verbose = verbose,
                                 use_Rcpp = use_Rcpp,
                                 num_threads = num_threads)
  # Function here
  if (verbose) print_init_ll(initloglik = initloglik)

  if (initloglik == -Inf) {
    stop("The initial parameter values have a likelihood that is
             equal to 0 or below machine precision.
             Try again with different initial values.")
  } else {
    out <- DDD::optimizer(optimmethod = optimmethod,
                          optimpars = optimpars,
                          fun = loglik_choosepar,
                          trparsopt = trparsopt,
                          num_cycles = num_cycles,
                          idparsopt = idparsopt,
                          trparsfix = trparsfix,
                          idparsfix = idparsfix,
                          idparslist = idparslist,
                          datalist = datalist,
                          num_observed_states = num_observed_states,
                          num_hidden_states = num_hidden_states,
                          cond = cond,
                          atol = atol,
                          rtol = rtol,
                          methode = methode,
                          rcpp_methode = rcpp_methode,
                          verbose = verbose,
                          use_Rcpp = use_Rcpp,
                          num_threads = num_threads)
    if (out$conv != 0) {
      stop("Optimization has not converged.
                 Try again with different initial values.")
    } else {
      ml_pars1 <- secsse_transform_parameters(as.numeric(unlist(out$par)),
                                              trparsfix,
                                              idparsopt,
                                              idparsfix,
                                              idparslist)
      out2 <- list(MLpars = ml_pars1,
                   ML = as.numeric(unlist(out$fvalues)),
                   conv = out$conv)
    }
  }
  return(out2)
}

#' loglik choosepar temp
#' @description
#' temporary export for testing
#' @inheritParams default_params_doc
#' @export
loglik_choosepar <- function(trparsopt,
                             trparsfix,
                             idparsopt,
                             idparsfix,
                             idparslist,
                              datalist,
                             num_observed_states,
                             num_hidden_states,
                             cond = cond,
                             atol,
                             rtol,
                             methode,
                             rcpp_methode,
                             verbose,
                             use_Rcpp,
                             num_threads) {
  alltrpars <- c(trparsopt, trparsfix)

  loglik <- NA

  if (max(alltrpars) > 1 || min(alltrpars) < 0) {
    loglik <- -Inf
  } else {
    pars1 <- secsse_transform_parameters(trparsopt, trparsfix,
                                         idparsopt, idparsfix,
                                         idparslist)

    loglik <- DAISIE_DE_trait_loglik_CS(parameter = pars1,
                                        datalist = datalist,
                                        methode = methode,
                                        rcpp_methode = rcpp_methode,
                                        atol = atol,
                                        rtol = rtol,
                                        num_observed_states = num_observed_states,
                                        num_hidden_states = num_hidden_states,
                                        cond = cond,
                                        verbose = verbose,
                                        use_Rcpp = use_Rcpp,
                                        num_threads = num_threads)

    if (is.nan(loglik) || is.na(loglik)) {
      warning("There are parameter values used which cause
                numerical problems.")
      loglik <- -Inf
    }
  }
  #if (verbose) {
  #  out_print <- c(trparsopt / (1 - trparsopt), loglik)
  #  message(paste(out_print, collapse = " "))
  #}
  return(loglik)
}
