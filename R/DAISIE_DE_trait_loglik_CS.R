#' apply loglikelihood to a dataset
#' @description umbrella function
#' @inheritParams default_params_doc
#' @export
DAISIE_DE_trait_loglik_CS <- function( parameter,
                                       datalist,
                                       methode = "lsodes",
                                       rcpp_methode = "odeint::bulirsch_stoer",
                                       atol = 1e-15,
                                       rtol = 1e-15,
                                       num_observed_states,
                                       num_hidden_states,
                                       cond = 1,
                                       num_threads = 1,
                                       verbose = FALSE,
                                       use_Rcpp = 0)

{
  logcond <- 0 # default value gives no effect

  if (length(parameter) == 6) {
    logp0 <- DAISIE_DE_trait_logp0(datalist = datalist,
                                   parameter = parameter,
                                   atol = atol,
                                   rtol = rtol,
                                   num_observed_states = num_observed_states,
                                   num_hidden_states = num_hidden_states,
                                   methode = methode,
                                   use_Rcpp = use_Rcpp)
    if (is.null(datalist[[1]]$not_present)) {
      loglik <- (datalist[[1]]$not_present_type1 + datalist[[1]]$not_present_type2) * logp0
      numimm <- (datalist[[1]]$not_present_type1 + datalist[[1]]$not_present_type2) + length(datalist) - 1
    } else {
      loglik <- datalist[[1]]$not_present * logp0
      numimm <- datalist[[1]]$not_present + length(datalist) - 1
    }

    ### condition on at least one successful colonization
    logcond <- (cond == 1) * log(1 - exp(numimm * logp0))
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
    trait_mainland_ancestor <- datalist[[i]]$root_state[2]
    sampling_fraction <- datalist[[i]]$sampling_fraction
    phy <- datalist[[i]]$phylogeny
    if (stac %in% c(1, 4)) {
      loglikelihood <- DAISIE_DE_trait_logpNE(brts = brts,
                                              status = stac,
                                              trait = trait,
                                              trait_mainland_ancestor = FALSE,
                                              num_observed_states = num_observed_states,
                                              num_hidden_states = num_hidden_states,
                                              parameter = parameter,
                                              atol  = atol,
                                              rtol  = rtol,
                                              methode = methode,
                                              rcpp_methode = rcpp_methode,
                                              use_Rcpp = use_Rcpp)
    } else if (stac %in% c(2, 5)) {
      if (length(brts) == 2) {
        loglikelihood <- DAISIE_DE_trait_logpES(brts = brts,
                                                status = stac,
                                                trait = trait,
                                                sampling_fraction = sampling_fraction,
                                                trait_mainland_ancestor = FALSE,
                                                num_observed_states = num_observed_states,
                                                num_hidden_states = num_hidden_states,
                                                parameter = parameter,
                                                atol  = atol,
                                                rtol  = rtol,
                                                methode = methode,
                                                rcpp_methode = rcpp_methode,
                                                use_Rcpp = use_Rcpp)
      } else {
        loglikelihood <- DAISIE_DE_trait_logpEC(brts = brts,
                                                parameter = parameter,
                                                phy = phy,
                                                traits = traits,
                                                num_observed_states = num_observed_states,
                                                num_hidden_states = num_hidden_states,
                                                trait_mainland_ancestor = FALSE,
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
        loglikelihood <- DAISIE_DE_trait_logpES(brts = brts,
                                                status = stac,
                                                trait = trait,
                                                sampling_fraction = sampling_fraction,
                                                trait_mainland_ancestor = FALSE,
                                                num_observed_states = num_observed_states,
                                                num_hidden_states = num_hidden_states,
                                                parameter = parameter,
                                                atol  = atol,
                                                rtol  = rtol,
                                                methode = methode,
                                                rcpp_methode = rcpp_methode,
                                                use_Rcpp = use_Rcpp)
      } else {
        loglikelihood <- DAISIE_DE_trait_logpEC( brts = brts,
                                                 parameter = parameter,
                                                 phy = phy,
                                                 traits = traits,
                                                 num_observed_states = num_observed_states,
                                                 num_hidden_states = num_hidden_states,
                                                 trait_mainland_ancestor = FALSE,
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
      loglikelihood <- DAISIE_DE_trait_logpEC(brts = brts,
                                              parameter = parameter,
                                              phy = phy,
                                              traits = traits,
                                              num_observed_states = num_observed_states,
                                              num_hidden_states = num_hidden_states,
                                              trait_mainland_ancestor = FALSE,
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
      loglikelihood <- DAISIE_DE_trait_logpNE_max_min_age_hidden(brts = brts,
                                                                 trait = trait,
                                                                 status = stac,
                                                                 parameter = parameter,
                                                                 num_observed_states = num_observed_states,
                                                                 num_hidden_states = num_hidden_states,
                                                                 atol  = atol,
                                                                 rtol  = rtol,
                                                                 methode = "ode45",
                                                                 rcpp_methode = rcpp_methode,
                                                                 use_Rcpp = use_Rcpp)
    } else if (stac == 9) {
      loglikelihood <- DAISIE_DE_trait_logpES_max_min_age_hidden(brts = brts,
                                                                 trait = trait,
                                                                 sampling_fraction = sampling_fraction,
                                                                 status = stac,
                                                                 parameter = parameter,
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

    vec_loglikelihood[i - 1] <- loglikelihood
  }

  loglik <- sum(vec_loglikelihood) + loglik
 # cat(vec_loglikelihood, loglik,"\n")
  return(loglik)
}
