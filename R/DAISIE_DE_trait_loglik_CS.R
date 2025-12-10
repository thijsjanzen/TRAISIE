#' apply loglikelihood to a dataset
#' @description umbrella function
#' @inheritParams default_params_doc
#' @export


DAISIE_DE_trait_loglik_CS <- function( parameter,
                                       datalist,
                                       methode = "lsodes",
                                       rcpp_methode =
                                       "odeint::runge_kutta_cash_karp54",
                                       atol = 1e-15,
                                       rtol = 1e-15,
                                       num_observed_states,
                                       num_hidden_states,
                                       cond = 1,
                                       num_threads = 1,
                                       verbose = FALSE,
                                       use_Rcpp = use_Rcpp)

{
  logcond <- 0 # default value gives no effect

  prob_dis_L0 <- c(datalist[[1]]$M0/datalist[[1]]$not_present, datalist[[1]]$M1/datalist[[1]]$not_present)
  if (length(parameter) >= 6) {
    logp0 <- DAISIE_DE_trait_logp0(datalist = datalist,
                                   parameter = parameter,
                                   atol = atol,
                                   rtol = rtol,
                                   num_observed_states = num_observed_states,
                                   num_hidden_states = num_hidden_states,
                                   trait_mainland_ancestor =  prob_dis_L0,
                                   methode = methode,
                                   use_Rcpp = use_Rcpp)
    if (is.null(datalist[[1]]$not_present)) {
      loglik <- (datalist[[1]]$not_present_type1 + datalist[[1]]$not_present_type2) * logp0$loglik
      numimm <- (datalist[[1]]$not_present_type1 + datalist[[1]]$not_present_type2) + length(datalist) - 1
    } else {
      loglik <- datalist[[1]]$not_present * logp0$loglik
      numimm <- datalist[[1]]$not_present + length(datalist) - 1
    }

    ### condition on at least one successful colonization
    logcond <- (cond == 1) * log(1 - exp(numimm * logp0$loglik))
    for (i in 2:length(datalist)) {
      datalist[[i]]$type1or2 <- 1
    }
  }

  loglik <- loglik - logcond

  vec_loglikelihood <- rep(NA, length(datalist) - 1) # first entry is not data
  for (i in 2:length(datalist)) {
    stac <- datalist[[i]]$stac
    brts <- datalist[[i]]$branching_times
    traits <- datalist[[i]]$traits
    trait <- datalist[[i]]$traits

    sampling_fraction <- datalist[[i]]$sampling_fraction

    phy <- datalist[[i]]$phylogeny

    trait_mainland_ancestor <- datalist[[i]]$root_state

    if (stac %in% c(1, 4)) {
      loglikelihood <- DAISIE_DE_trait_logpNE(datalist = datalist,
                                              brts = brts,
                                              status = stac,
                                              trait = trait,
                                              trait_mainland_ancestor = trait_mainland_ancestor,
                                              num_observed_states = num_observed_states,
                                              num_hidden_states = num_hidden_states,
                                              sampling_fraction = sampling_fraction,
                                              parameter = parameter,
                                              atol  = atol,
                                              rtol  = rtol,
                                              methode = methode,
                                              rcpp_methode = rcpp_methode,
                                              use_Rcpp = use_Rcpp)
    } else if (stac %in% c(2, 5)) {
      if (length(brts) == 2) {
        loglikelihood <- DAISIE_DE_trait_logpES(datalist = datalist,
                                                brts = brts,
                                                status = stac,
                                                trait = trait,
                                                sampling_fraction = sampling_fraction,
                                                trait_mainland_ancestor = trait_mainland_ancestor,
                                                num_observed_states = num_observed_states,
                                                num_hidden_states = num_hidden_states,
                                                parameter = parameter,
                                                atol  = atol,
                                                rtol  = rtol,
                                                methode = methode,
                                                rcpp_methode = rcpp_methode,
                                                use_Rcpp = use_Rcpp)
      } else {
        loglikelihood <- DAISIE_DE_trait_logpEC(datalist = datalist,
                                                brts = brts,
                                                parameter = parameter,
                                                phy = phy,
                                                traits = traits,
                                                num_observed_states = num_observed_states,
                                                num_hidden_states = num_hidden_states,
                                                trait_mainland_ancestor = trait_mainland_ancestor,
                                                status = stac,
                                                sampling_fraction = sampling_fraction,
                                                atol  = atol,
                                                rtol  = rtol,
                                                methode = methode,
                                                use_Rcpp = use_Rcpp,
                                                rcpp_methode = rcpp_methode,
                                                num_threads = num_threads)
      }
    } else if (stac == 3) {
      if (length(brts) == 2) {
        loglikelihood <- DAISIE_DE_trait_logpES(datalist = datalist,
                                                brts = brts,
                                                status = stac,
                                                trait = trait,
                                                sampling_fraction = sampling_fraction,
                                                trait_mainland_ancestor = trait_mainland_ancestor,
                                                num_observed_states = num_observed_states,
                                                num_hidden_states = num_hidden_states,
                                                parameter = parameter,
                                                atol  = atol,
                                                rtol  = rtol,
                                                methode = methode,
                                                rcpp_methode = rcpp_methode,
                                                use_Rcpp = use_Rcpp)
      } else {
        loglikelihood <- DAISIE_DE_trait_logpEC( datalist = datalist,
                                                 brts = brts,
                                                 parameter = parameter,
                                                 phy = phy,
                                                 traits = traits,
                                                 num_observed_states = num_observed_states,
                                                 num_hidden_states = num_hidden_states,
                                                 trait_mainland_ancestor = trait_mainland_ancestor,
                                                 status = stac,
                                                 sampling_fraction = sampling_fraction,
                                                 atol  = atol,
                                                 rtol  = rtol,
                                                 methode = methode,
                                                 rcpp_methode = rcpp_methode,
                                                 use_Rcpp = use_Rcpp)
      }
    }
    else if (stac == 6) {
      loglikelihood <- DAISIE_DE_trait_logpEC(datalist = datalist,
                                              brts = brts,
                                              parameter = parameter,
                                              phy = phy,
                                              traits = traits,
                                              num_observed_states = num_observed_states,
                                              num_hidden_states = num_hidden_states,
                                              trait_mainland_ancestor = trait_mainland_ancestor,
                                              status = stac,
                                              sampling_fraction = sampling_fraction,
                                              atol  = atol,
                                              rtol  = rtol,
                                              methode = methode,
                                              use_Rcpp = use_Rcpp,
                                              rcpp_methode = rcpp_methode,
                                              num_threads = num_threads)
    }
    else if (stac == 8) {
      loglikelihood <-
        DAISIE_DE_trait_logpNE_max_min_age_hidden(datalist = datalist,
                                                  brts = brts,
                                                  trait = trait,
                                                  status = stac,
                                                  parameter = parameter,
                                                  trait_mainland_ancestor = trait_mainland_ancestor,
                                                  num_observed_states = num_observed_states,
                                                  num_hidden_states = num_hidden_states,
                                                  sampling_fraction = sampling_fraction,
                                                  atol  = atol,
                                                  rtol  = rtol,
                                                  methode = "ode45",
                                                  rcpp_methode = rcpp_methode,
                                                  use_Rcpp = use_Rcpp)
    } else if (stac == 9) {
      loglikelihood <-
        DAISIE_DE_trait_logpES_max_min_age_hidden(datalist = datalist,
                                                  brts = brts,
                                                  trait = trait,
                                                  sampling_fraction = sampling_fraction,
                                                  status = stac,
                                                  parameter = parameter,
                                                  trait_mainland_ancestor = trait_mainland_ancestor,
                                                  num_observed_states = num_observed_states,
                                                  num_hidden_states = num_hidden_states,
                                                  atol  = atol,
                                                  rtol  = rtol,
                                                  methode = "ode45",
                                                  rcpp_methode = rcpp_methode,
                                                  use_Rcpp = use_Rcpp)
    } else {
      stop("Unknown stac value: ", stac)
    }

    vec_loglikelihood[i - 1] <- loglikelihood$loglik
  }

  loglik <- sum(vec_loglikelihood) + loglik
  return(loglik)
  #return(c(vec_loglikelihood, loglik))
}

