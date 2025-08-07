#' @keywords internal
loglik_hidden_rhs <- function(t, state, parameter) {
  with(as.list(c(state, parameter)), {

    lambdac <- parameter[[1]]
    mu      <- parameter[[2]]
    gamma   <- parameter[[3]]
    lambdaa <- parameter[[4]]
    q       <- parameter[[5]]
    p       <- parameter[[6]]


    n <- (length(state) - 1) / 3

    dDE     <- numeric(n)
    dDM3    <- numeric(n)
    dE      <- numeric(n)

    t_vec   <- rowSums(q)

    DE  <- state[1:n]
    DM3 <- state[(n + 1):(n + n)]
    E   <- state[(n + n + 1):(n + n + n)]
    DA3 <- state[length(state)]

    gamma_matrix <- matrix(parameter[[3]], nrow = n, ncol = length(parameter[[3]]), byrow = TRUE)

    if (trait_mainland_ancestor == FALSE) {
      gamma <- gamma  # do nothing
      dist_gamma <- gamma/n
    } else {
      s <- ((trait_mainland_ancestor * num_hidden_states) + 1) :
        ((trait_mainland_ancestor + 1) * num_hidden_states)
      gamma <- gamma[s]
      dist_gamma <- gamma/num_hidden_states
    }

    q_mult_E   <- t(q %*% E)
    q_mult_DE  <- t(q %*% DE)
    q_mult_DM3 <- t(q %*% DM3)


    lambda_c_mu_t_vec_sum <- lambdac + mu + t_vec
    E_sq <- E * E

    dDE <- -(lambda_c_mu_t_vec_sum) * DE +
      2 * lambdac * DE * E +
      q_mult_DE

    dDM3 <-  -(lambda_c_mu_t_vec_sum + dist_gamma + lambdaa) * DM3 +
      (mu + lambdaa * E + lambdac * E_sq + p * q_mult_E) * DA3 +
      (1 - p) * q_mult_DM3 +
      dist_gamma * DM3

    dE <- mu - (lambda_c_mu_t_vec_sum) * E +
      lambdac * E_sq +
      q_mult_E

    dDA3 <- -sum(dist_gamma) * DA3 + sum(dist_gamma * DM3)

    return(list(c(dDE, dDM3, dE, dDA3)))
  })
}


#' @keywords internal
calcThruNodes_hidden <- function(
    ances,
    states,
    loglik,
    forTime,
    parameter,
    methode,
    phy,
    rhs_func,
    reltol,
    abstol
) {

  nb_node <- phy$Nnode
  focal <- ances
  desRows <- which(phy$edge[, 1] == focal)
  desNodes <- phy$edge[desRows, 2]

  nodeM <- numeric()
  nodeN <- numeric()

  for (desIndex in 1:2) {
    y <- states[desNodes[desIndex],]
    #
    timeInte <- forTime[which(forTime[, 2] == desNodes[desIndex]), 3]
    ##  To do the calculation in both lineages

    nodeMN <- deSolve::ode(y = y,
                           func = rhs_func,
                           times = c(0, timeInte),
                           parms = parameter,
                           rtol = reltol,
                           atol = abstol,
                           method = methode)

    if (desIndex == 1) {
      nodeN <- nodeMN
    }
    if (desIndex == 2) {
      nodeM <- nodeMN
    }
  }
  ## At the node
  nodeM <- as.numeric(nodeM[2,-1])
  nodeN <- as.numeric(nodeN[2,-1])


  lambda_c <- parameter[[1]]
  n <- length(lambda_c)

  DE_N <- nodeN[1:n]
  DE_M <- nodeM[1:n]

  combined_state <- nodeN
  combined_state[1:n] <- lambda_c * DE_N * DE_M

  states[focal,] <- combined_state
  return(list(states = states,
              loglik = loglik,
              combined_state = combined_state,
              nodeM = nodeM,
              nodeN = nodeN))
}
  # sf = sampling fraction
#' @keywords internal
calc_init_state_hidden <- function(trait,
                                   sampling_fraction,
                                   num_unique_states,
                                   num_hidden_states,
                                   mainland = FALSE,
                                   trait_mainland_ancestor) {

  DE  <- rep(0, num_unique_states)
  DM3 <- rep(0, num_unique_states)
  E   <- rep(0, num_unique_states)
  DA3 <- 1


  if (trait == "FALSE") {
    DE[c(1, n)] <- sampling_fraction
    E[c(1, n)]  <- 1 - sampling_fraction
  }

  else if (trait == trait) {
    DE[c((num_hidden_states * trait + 1), num_hidden_states + trait * num_hidden_states)] <- sampling_fraction
    E[c((num_hidden_states * trait + 1), num_hidden_states + trait * num_hidden_states)] <- 1 - sampling_fraction

  }


  if (mainland) {
    DM3[c((num_hidden_states * trait_mainland_ancestor + 1),
           num_hidden_states + trait_mainland_ancestor * num_hidden_states)] <- 1
  }


  return( c(DE, DM3, E, DA3))
  }


#' Likelihood calculation including hidden traits
#' @title Using hidden traits
#'
#' @inheritParams default_params_doc
#'
#' @param rhs_func ll function
#' @export
loglik_R_tree <- function(parameter,
                          phy,
                          traits,
                          sampling_fraction,
                          num_hidden_states,
                          mainland = FALSE,
                          trait_mainland_ancestor = NULL,
                          atol = 1e-15,
                          rtol = 1e-15,
                          methode = "ode45",
                          rhs_func = loglik_hidden_rhs) {

  number_of_lineages <- length(phy$tip.label)
  num_unique_states <- length(parameter[[1]])
  states <- matrix(nrow = number_of_lineages + phy$Nnode,
                   ncol = 3 * length(parameter[[1]]) + 1,
                   data = NA)

  for (i in 1:length(traits)) {
    states[i, ] <- calc_init_state_hidden(trait = traits[i],
                                          sampling_fraction = sampling_fraction[ 1 + traits[i] ],
                                          num_unique_states = num_unique_states,
                                          num_hidden_states = num_hidden_states,
                                          mainland = mainland,
                                          trait_mainland_ancestor = trait_mainland_ancestor)
  }

  phy$node.label <- NULL
  split_times <- sort(secsse::event_times(phy), decreasing = FALSE)
  ances <- as.numeric(names(split_times))
  forTime <- cbind(phy$edge, phy$edge.length)

  loglik <- 0

  for (i in 1:length(ances)) {
    calcul <- calcThruNodes_hidden(ances = ances[i],
                                   states = states,
                                   loglik = loglik,
                                   forTime = forTime,
                                   parameter = parameter,
                                   methode = methode,
                                   phy = phy,
                                   rhs_func = rhs_func,
                                   reltol = rtol,
                                   abstol = atol)
    states <- calcul$states
    loglik <- calcul$loglik
    nodeN <- calcul$nodeN
  }

  prob_states <- calcul$combined_state
  prob_states <- matrix(prob_states, nrow = 1)
  return(prob_states)
}




#' Likelihood calculation including hidden traits
#' @title Using hidden traits
#'
#' @inheritParams default_params_doc
#'
#' @export
loglik_cpp_tree <- function(parameter,
                            phy,
                            traits,
                            sampling_fraction,
                            num_hidden_states,
                            mainland = FALSE,
                            trait_mainland_ancestor = NULL,
                            atol = 1e-15,
                            rtol = 1e-15,
                            method = "odeint::bulirsch_stoer",
                            use_normalization = TRUE,
                            num_threads = 1) {

  number_of_lineages <- length(phy$tip.label)
  num_unique_states <- length(parameter[[1]])
  states <- matrix(nrow = number_of_lineages + phy$Nnode,
                   ncol = 3 * length(parameter[[1]]) + 1,
                   data = NA)
  # sf = sampling fraction
  for (i in 1:length(traits)) {
    states[i, ] <- calc_init_state_hidden(traits[i],
                                          sampling_fraction[ 1 + traits[i] ],
                                          num_unique_states,
                                          num_hidden_states,
                                          mainland = mainland,
                                          trait_mainland_ancestor = trait_mainland_ancestor)
  }

  phy$node.label <- NULL
  split_times <- sort(secsse::event_times(phy), decreasing = FALSE)
  ances <- as.numeric(names(split_times))
  forTime <- cbind(phy$edge, phy$edge.length)

  lambda_c <- parameter[[1]]
  mus      <- parameter[[2]]
  gammas   <- parameter[[3]]
  lambda_a <- parameter[[4]]
  q_matrix       <- parameter[[5]]
  p_value       <- parameter[[6]]

  RcppParallel::setThreadOptions(numThreads = num_threads)

  calcul <- calc_ll_cpp(ances = ances,
                        states = states,
                        forTime = forTime,
                        lambda_cs = lambda_c,
                        lambda_as = lambda_a,
                        mus = mus,
                        gammas = gammas,
                        qs = q_matrix,
                        p = p_value,
                        method = method,
                        atol = atol,
                        rtol = rtol,
                        see_states = TRUE,
                        use_normalization = use_normalization)


  prob_states <- calcul$merge_branch
  prob_states <- matrix(prob_states, nrow = 1)
  return(prob_states)
}

