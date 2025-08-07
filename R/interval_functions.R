#' @keywords internal
interval2 <- function(t, state, parameter) {
  with(as.list(c(state, parameter)), {
    lambdac <- parameter[[1]]
    mu      <- parameter[[2]]
    gamma   <- parameter[[3]]
    lambdaa <- parameter[[4]]
    q       <- parameter[[5]]
    p       <- parameter[[6]]

    n <- (length(state) - 1) / 4

    dDE     <- numeric(n)
    dDM2    <- numeric(n)
    dDM3    <- numeric(n)
    dE      <- numeric(n)

    t_vec <- rowSums(q)

    DE  <- state[1:n]
    DM2 <- state[(n + 1):(n + n)]
    DM3 <- state[(n + n + 1):(n + n + n)]
    E   <- state[(n + n + n + 1):(n + n + n + n)]
    DA3 <- state[length(state)]

    gamma_matrix <- matrix(parameter[[3]], nrow = n, ncol = length(parameter[[3]]), byrow = TRUE)



    q_mult_E   <- t(q %*% E)
    q_mult_DE  <- t(q %*% DE)
    q_mult_DM2 <- t(q %*% DM2)
    q_mult_DM3 <- t(q %*% DM3)


    if (trait_mainland_ancestor == FALSE) {
      gamma <- gamma  # do nothing
      dist_gamma <- gamma/n
    } else {
      s <- ((trait_mainland_ancestor * num_hidden_states) + 1) :
        ((trait_mainland_ancestor + 1) * num_hidden_states)
      gamma <- gamma[s]
      dist_gamma <- gamma/num_hidden_states
    }


    dDE <- -(lambdac + mu + t_vec) * DE +
      2 * lambdac * DE * E +
      q_mult_DE

    dDM2 <- -(lambdac + mu + sum(dist_gamma) + lambdaa + t_vec) * DM2 +
      (lambdaa * DE + 2 * lambdac * DE * E + p * q_mult_DE) * DA3 +
      (1 - p) * q_mult_DM2

    dDM3 <- -(lambdac + mu + sum(dist_gamma) + lambdaa + t_vec) * DM3 +
      (mu + lambdaa * E + lambdac * E * E + p * q_mult_E) * DA3 +
      (1 - p) * q_mult_DM3 + sum(dist_gamma * DM3)

    dE <- mu - (mu + lambdac + t_vec) * E +
      lambdac * E * E +
      q_mult_E

    dDA3 <- -sum(dist_gamma) * DA3 + sum(dist_gamma * DM3)
    return(list(c(dDE, dDM2, dDM3, dE, dDA3)))
  })
}

#' @keywords internal
interval3 <- function(t, state, parameter) {
  with(as.list(c(state, parameter)), {
    lambdac <- parameter[[1]]
    mu      <- parameter[[2]]
    gamma   <- parameter[[3]]
    lambdaa <- parameter[[4]]
    q       <- parameter[[5]]
    p       <- parameter[[6]]

    n <- (length(state) - 2) / 5

    dDE     <- numeric(n)
    dDM1    <- numeric(n)
    dDM2    <- numeric(n)
    dDM3    <- numeric(n)
    dE      <- numeric(n)

    t_vec <- rowSums(q)

    DE  <- state[1:n]
    DM1 <- state[(n + 1):(n + n)]
    DM2 <- state[(n + n + 1):(n + n + n)]
    DM3 <- state[(n + n + n + 1):(n + n + n + n)]
    E   <- state[(n + n + n + n + 1):(n + n + n + n + n)]
    DA2 <- state[length(state) - 1]
    DA3 <- state[length(state)]

    gamma_matrix <- matrix(parameter[[3]], nrow = n, ncol = length(parameter[[3]]), byrow = TRUE)



    q_mult_E   <- t(q %*% E)
    q_mult_DE  <- t(q %*% DE)
    q_mult_DM1 <- t(q %*% DM1)
    q_mult_DM2 <- t(q %*% DM2)
    q_mult_DM3 <- t(q %*% DM3)

    if (trait_mainland_ancestor == FALSE) {
      gamma <- gamma  # do nothing
      dist_gamma <- gamma/n
    } else {
      s <- ((trait_mainland_ancestor * num_hidden_states) + 1) :
        ((trait_mainland_ancestor + 1) * num_hidden_states)
      gamma <- gamma[s]
      dist_gamma <- gamma/num_hidden_states
    }


    dDE <- -(lambdac + mu + t_vec) * DE +
      2 * lambdac * DE * E +
      q_mult_DE

    dDM1 <- -(lambdac + mu + sum(dist_gamma) + lambdaa + t_vec) * DM1 +
      (mu + lambdaa * E + lambdac * E * E + p * q_mult_E) * DA2 +
      (1 - p) * q_mult_DM1 + sum(dist_gamma * DM2)

    dDM2 <- -(lambdac + mu + sum(dist_gamma) + lambdaa + t_vec) * DM2 +
      (mu + lambdaa * E + lambdac * E * E + p * q_mult_E) * DA2 +
      (lambdaa * DE + 2 * lambdac * DE + p * q_mult_DE) * DA3 +
      (1 - p) * q_mult_DM2 + sum(dist_gamma * DM2)

    dDM3 <- -(lambdac + mu + sum(dist_gamma) + lambdaa + t_vec) * DM3 +
      (mu + lambdaa * E + lambdac * E * E + p * q_mult_E) * DA3 +
      (1 - p) * q_mult_DM3 + sum(dist_gamma * DM3)

    dE <- mu - (mu + lambdac + t_vec) * E +
      lambdac * E * E +
      q_mult_E

    dDA2 <- -sum(dist_gamma) * DA2 + sum(dist_gamma * DM2)
    dDA3 <- -sum(dist_gamma) * DA3 + sum(dist_gamma * DM3)

    return(list(c(dDE, dDM1, dDM2, dDM3, dE, dDA2, dDA3)))
  })
}

#' @keywords internal
interval4 <- function(t, state, parameter) {
  with(as.list(c(state, parameter)), {
    lambdac <- parameter[[1]]
    mu      <- parameter[[2]]
    gamma   <- parameter[[3]]
    lambdaa <- parameter[[4]]
    q       <- parameter[[5]]
    p       <- parameter[[6]]

    n <- (length(state) - 1) / 2

    dDM1 <- numeric(n)
    dDE  <- numeric(n)

    t_vec <- rowSums(q)

    DM1 <- state[1:n]
    E   <- state[(n + 1):(n + n)]
    DA1 <- state[length(state)]

    gamma_matrix <- matrix(parameter[[3]], nrow = n, ncol = length(parameter[[3]]), byrow = TRUE)



    q_mult_E   <- t(q %*% E)
    q_mult_DM1 <- t(q %*% DM1)


    if (trait_mainland_ancestor == FALSE) {
      gamma <- gamma  # do nothing
      dist_gamma <- gamma/n
    } else {
      s <- ((trait_mainland_ancestor * num_hidden_states) + 1) :
        ((trait_mainland_ancestor + 1) * num_hidden_states)
      gamma <- gamma[s]
      dist_gamma <- gamma/num_hidden_states
    }


    dDM1 <- -(lambdac + mu + sum(dist_gamma) + lambdaa + t_vec) * DM1 +
      (mu + lambdaa * E + lambdac * E * E + p * q_mult_E) * DA1 +
      (1 - p) * q_mult_DM1  + sum(dist_gamma * DM1)

    dE <- mu - (mu + lambdac + t_vec) * E +
      lambdac * E * E +
      q_mult_E

    dDA1 <- - sum(dist_gamma) * DA1 + sum(dist_gamma * DM1)

    return(list(c(dDM1, dE, dDA1)))
  })
}

