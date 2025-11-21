
DAISIE_DE_trait_logpEC1(
  brts                    = brts,
  phy                     = phy,
  traits                  = traits,
  status                  = 2,
  sampling_fraction       = sampling_fraction,
  parameter               = parameter,
  trait_mainland_ancestor = c(0,1),
  Mainland_pool_size_vec = c(3000, 3000),   # <--- NEW ARGUMENT
  num_observed_states     = 2,
  num_hidden_states       = 2,
  atol                    = 1e-15,
  rtol                    = 1e-15,
  methode                 = "ode45",
  use_Rcpp                = 2)



 DAISIE_DE_trait_logpEC1 <- function(
    brts,
    parameter,
    phy,
    traits,
    num_observed_states,
    num_hidden_states,
    trait_mainland_ancestor = NA,
    status,
    sampling_fraction,
    Mainland_pool_size_vec = NULL,
    num_threads = 1,
    atol = 1e-15,
    rtol = 1e-15,
    methode = "ode45",
    rcpp_methode = "odeint::runge_kutta_cash_karp54",
    use_Rcpp = 0
 ) {

   # total mainland pool
   M <- sum(Mainland_pool_size_vec)

   # ------------------------------------------------------------------
   # CASE 1 : THE MAINLAND TRAIT IS KNOWN (E.G., c(1,0))
   # ------------------------------------------------------------------

   if (!all(is.na(trait_mainland_ancestor))) {

     obs_state <- which(trait_mainland_ancestor == 1)

     Lk_vec <- numeric(num_hidden_states)

     # Loop over hidden states
     for (h in seq_len(num_hidden_states)) {

       # Create vector of length = observed × hidden
       trait_mainland_extended <- rep(0, num_observed_states * num_hidden_states)

       # index = observed_state block + hidden_state
       index <- (obs_state - 1) * num_hidden_states + h
       trait_mainland_extended[index] <- 1

       # Compute likelihood for this hidden state
       Lk_h_log <- DAISIE_DE_trait_logpEC_core(
         brts                    = brts,
         parameter               = parameter,
         phy                     = phy,
         traits                  = traits,
         num_observed_states     = num_observed_states,
         num_hidden_states       = num_hidden_states,
         trait_mainland_ancestor = trait_mainland_extended,
         status                  = status,
         sampling_fraction       = sampling_fraction,
         num_threads             = num_threads,
         atol                    = atol,
         rtol                    = rtol,
         methode                 = methode,
         rcpp_methode            = rcpp_methode,
         use_Rcpp                = use_Rcpp
       )

       # weighted contribution of (obs_state, hidden_state = h)
       Lk_vec[h] <- exp(Lk_h_log) * Mainland_pool_size_vec[obs_state] / (M*num_hidden_states)
     }

     # sum over hidden states
     return(log(sum(Lk_vec)))
   }



   # ---------------------------------------------------------------------
   # CASE 2 : TRAIT OF MAINLAND ANCESTOR UNKNOWN (trait_mainland_ancestor = NA)
   # ---------------------------------------------------------------------
   if (all(is.na(trait_mainland_ancestor))) {

     Lk_obs_vec <- numeric(num_observed_states)

     # Loop over observed states
     for (i in seq_len(num_observed_states)) {

       Lk_hidden_vec <- numeric(num_hidden_states)

       # Loop over hidden states
       for (h in seq_len(num_hidden_states)) {

         # Create extended vector of length = observed × hidden
         trait_mainland_extended <- rep(0, num_observed_states * num_hidden_states)

         index <- (i - 1) * num_hidden_states + h
         trait_mainland_extended[index] <- 1

         # Compute likelihood for this (observed=i, hidden=h)
         Lk_h_log <- DAISIE_DE_trait_logpEC_core(
           brts                    = brts,
           parameter               = parameter,
           phy                     = phy,
           traits                  = traits,
           num_observed_states     = num_observed_states,
           num_hidden_states       = num_hidden_states,
           trait_mainland_ancestor = trait_mainland_extended,
           status                  = status,
           sampling_fraction       = sampling_fraction,
           num_threads             = num_threads,
           atol                    = atol,
           rtol                    = rtol,
           methode                 = methode,
           rcpp_methode            = rcpp_methode,
           use_Rcpp                = use_Rcpp
         )

         Lk_hidden_vec[h] <- exp(Lk_h_log)
       }

       # Sum over hidden states and apply weight of observed trait i
       Lk_obs_vec[i] <- sum(Lk_hidden_vec) * Mainland_pool_size_vec[i] / (M*num_hidden_states)
     }

     # Sum over observed states
     return(log(sum(Lk_obs_vec)))
   }

 }




 DAISIE_DE_trait_logpEC_core <- function(
    brts,
    parameter,
    phy,
    traits,
    num_observed_states,
    num_hidden_states,
    trait_mainland_ancestor = NA,
    status,
    sampling_fraction,
    num_threads = 1,
    atol = 1e-15,
    rtol = 1e-15,
    methode = "ode45",
    rcpp_methode = "odeint::runge_kutta_cash_karp54",
    use_Rcpp = 0
 ) {



   check_arguments(brts, parameter, phy, traits, num_observed_states,
                   num_hidden_states, status, sampling_fraction)



   if (length(brts) < 3) {
     stop("need at least three branching times")
   }

   # Unpack times from brts
   t0   <- brts[1]
   tmax <- brts[2]
   t1   <- brts[2]
   t2   <- brts[3]
   tp   <- 0

   # Time intervals

   time2 <- c(t2, t1)
   time3 <- c(t2, tmax)
   time4 <- c(tmax, t0)

   # Number of states in the system
   #n <- num_observed_states * num_hidden_states

   # Solve for interval [tp, t2] (stem phase)
   res <- c()

   if (length(phy$tip.label) < 2) {
     stop("Tip too small to calculate tree likelihood")
   }

   if (use_Rcpp == 0) {
     res <- loglik_R_tree(
       parameter = parameter,
       phy = phy,
       traits = traits,
       sampling_fraction = sampling_fraction,
       num_hidden_states = num_hidden_states,
       trait_mainland_ancestor = trait_mainland_ancestor,
       atol = atol,
       rtol = rtol
     )
   } else {
     res <- loglik_cpp_tree(
       parameter = parameter,
       phy = phy,
       traits = traits,
       sampling_fraction = sampling_fraction,
       num_hidden_states = num_hidden_states,
       trait_mainland_ancestor = trait_mainland_ancestor,
       atol = atol,
       rtol = rtol,
       num_threads = num_threads
     )
   }

   # Run appropriate sequence of intervals
   if ((status == 2 || status == 3) && length(brts) > 2) {

     initial_conditions2 <- get_initial_conditions2(status = status,
                                                    res = res,
                                                    trait = traits,
                                                    num_observed_states = num_observed_states,
                                                    num_hidden_states = num_hidden_states,
                                                    brts = brts,
                                                    sampling_fraction = sampling_fraction,
                                                    trait_mainland_ancestor = trait_mainland_ancestor)

     solution2 <- solve_branch(interval_func = interval2,
                               initial_conditions = initial_conditions2,
                               time = time2,
                               parameter = parameter,
                               trait_mainland_ancestor = trait_mainland_ancestor,
                               methode = methode,
                               rcpp_methode = rcpp_methode,
                               atol = atol,
                               rtol =  rtol,
                               use_Rcpp = use_Rcpp)

     initial_conditions4 <- get_initial_conditions4(status = status,
                                                    solution = solution2,
                                                    parameter = parameter,
                                                    trait_mainland_ancestor = trait_mainland_ancestor,
                                                    num_observed_states = num_observed_states,
                                                    num_hidden_states = num_hidden_states)

     solution4 <- solve_branch(interval_func = interval4,
                               initial_conditions = initial_conditions4,
                               time = time4,
                               parameter = parameter,
                               trait_mainland_ancestor = trait_mainland_ancestor,
                               methode = methode,
                               rcpp_methode = rcpp_methode,
                               atol = atol,
                               rtol = rtol,
                               use_Rcpp = use_Rcpp)
   }

   if (status == 6) {
     initial_conditions3 <- get_initial_conditions3(status = status,
                                                    res = res,
                                                    num_observed_states = num_observed_states,
                                                    num_hidden_states = num_hidden_states,
                                                    trait = traits,
                                                    sampling_fraction = sampling_fraction)
     solution3 <- solve_branch(interval_func = interval3,
                               initial_conditions = initial_conditions3,
                               time = time3,
                               parameter = parameter,
                               trait_mainland_ancestor = trait_mainland_ancestor,
                               methode = methode,
                               rcpp_methode = rcpp_methode,
                               atol = atol,
                               rtol = rtol,
                               use_Rcpp = use_Rcpp)


     initial_conditions4 <- get_initial_conditions4(status = status,
                                                    solution = solution3,
                                                    parameter = parameter,
                                                    trait_mainland_ancestor = trait_mainland_ancestor,
                                                    num_observed_states = num_observed_states,
                                                    num_hidden_states = num_hidden_states)
     solution4 <- solve_branch(interval_func = interval4,
                               initial_conditions = initial_conditions4,
                               time = time4,
                               parameter = parameter,
                               trait_mainland_ancestor = trait_mainland_ancestor,
                               methode = methode,
                               rcpp_methode = rcpp_methode,
                               atol = atol,
                               rtol = rtol,
                               use_Rcpp = use_Rcpp)
   }

   # Extract log-likelihood from final solution
   Lk <- solution4[2, length(solution4[2, ])]
   logLkb <- log(Lk)
   return(logLkb)
 }
