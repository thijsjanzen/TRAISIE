
# initial conditions system of equation interval 2
#' @keywords internal
get_initial_conditions2 <- function(status,
                                    res,
                                    trait,
                                    num_observed_states,
                                    num_hidden_states,
                                    brts,
                                    sampling_fraction,
                                    trait_mainland_ancestor) {
  n <- num_observed_states * num_hidden_states
  num_unique_states <- n

  initial_conditions2 <- c()

  DE  <- rep(0, num_unique_states)
  DM2 <- rep(0, num_unique_states)
  DM3 <- rep(0, num_unique_states)
  E   <- rep(0, num_unique_states)
  DA3 <- 1

  if (status == 2 && length(brts) > 2) {
    initial_conditions2 <- c(res[1:n],                      ## DE
                             (res[1:n]) * res[length(res)], ## DM2
                             res[(n + 1):(n + n)],          ## DM3
                             res[(n + n + 1):(n + n + n)],  ## E
                             res[length(res)])              ## DA3
    # pre-emptive return because this one is constructed differently.
    return(matrix(initial_conditions2, nrow = 1))
  } else if (status == 2 && length(brts) == 2) {

    if (trait == "FALSE") {
      DE[c(1, n)] <- sampling_fraction
      E[c(1, n)]  <- 1 - sampling_fraction
    }

    else if (trait == trait) {
      DE[c((num_hidden_states * trait + 1), num_hidden_states + trait * num_hidden_states)] <- sampling_fraction
      E[c((num_hidden_states * trait + 1), num_hidden_states + trait * num_hidden_states)]  <- 1 - sampling_fraction
    }

  }
  else if (status == 3 && length(brts) == 2 ) {
    DE[c((num_hidden_states*trait + 1), num_hidden_states + trait * num_hidden_states)] <- sampling_fraction
    E[c((num_hidden_states*trait + 1), num_hidden_states + trait * num_hidden_states)] <- 1 - sampling_fraction
    DM3[c((num_hidden_states * trait_mainland_ancestor + 1),
          num_hidden_states + trait_mainland_ancestor * num_hidden_states)] <- 1
  }
  else if (status == 4) {

    if (trait == "FALSE") {
      DM2[c(1, n)] <- 1
    }

    else if (trait == trait) {
      DM2[c((num_hidden_states * trait + 1), num_hidden_states + trait * num_hidden_states)] <- 1
    }

  }
  else if(status == 8) {
    if (trait == "FALSE") {
      DM2[c(1, n)] <- 1
    }

    else if (trait == trait) {
      DM2[c((num_hidden_states * trait + 1), num_hidden_states + trait * num_hidden_states)] <- 1
    }
  }
  else if(status == 9)  {
    if (trait == "FALSE") {
      DE[c(1, n)] <- sampling_fraction
      E[c(1, n)]  <- 1 - sampling_fraction
    }

    else if (trait == trait) {
      DE[c((num_hidden_states * trait + 1), num_hidden_states + trait * num_hidden_states)] <- sampling_fraction
      E[c((num_hidden_states * trait + 1), num_hidden_states + trait * num_hidden_states)]  <- 1 - sampling_fraction
    }
  }

  initial_conditions2 <- ( c(DE, DM2, DM3, E, DA3))
  return(matrix(initial_conditions2, nrow = 1))
}


# initial conditions system of equation interval 3
#' @keywords internal
get_initial_conditions3 <- function(status,
                                    res,
                                    num_observed_states,
                                    num_hidden_states,
                                    trait,
                                    solution) {

  n <- num_observed_states * num_hidden_states
  num_unique_states <- n

  initial_conditions3 <- c()
  DE  <- rep(0, num_unique_states)
  DM2 <- rep(0, num_unique_states)
  DM1 <- rep(0, num_unique_states)
  DM3 <- rep(0, num_unique_states)
  E   <- rep(0, num_unique_states)
  DA2 <- 0
  DA3 <- 1

  if (status == 1) {
    if (trait == "FALSE") {
      DM2[c(1, n)] <- 1
    }

    else if (trait == trait) {
      DM2[c((num_hidden_states * trait + 1), num_hidden_states + trait * num_hidden_states)] <- 1
    }

    initial_conditions3 <- c(DE, DM1, DM2, DM3, E, DA2, DA3)
  } else if (status == 6) {
    initial_conditions3 <- c(res[1:n],                                              ## DE
                             rep(0, n),                                             ## DM1
                             (res[1:n]) * res[length(res)],                         ## DM2
                             res[(n + 1):(n + n)],                                  ## DM3
                             res[(n + n + 1):(n + n + n)],                          ## E
                             0,                                                     ## DA2
                             res[length(res)])                                      ## DA3
  }
  else if (status == 5) {
    if (trait == "FALSE") {
      DE[c(1, n)] <- sampling_fraction
      E[c(1, n)]  <- 1 - sampling_fraction
    }

    else if (trait == trait) {
      DE[c((num_hidden_states * trait + 1), num_hidden_states + trait * num_hidden_states)] <- sampling_fraction
      E[c((num_hidden_states * trait + 1), num_hidden_states + trait * num_hidden_states)]  <- 1 - sampling_fraction
    }

    initial_conditions3 <- c(DE, DM1, DM2, DM3, E, DA2, DA3)
  } else if (status == 8 || status == 9) {
    initial_conditions3 <- c(solution[2, ][1:n],
                             rep(0, n),                                             ### DE: select DE in solution2
                             solution[2, ][(n + 1):(n + n)],                        ### DM2: select DM2 in solution2
                             solution[2, ][(n + n + 1):(n + n + n)],                ### DM3: select DM3 in solution2
                             solution[2, ][(n + n + n + 1):(n + n + n + n)],        ### E: select E in solution2
                             0,
                             solution[2, ][length(solution[2, ])])                 ### DA3: select DA3 in solution2
  }

  return(matrix(initial_conditions3, nrow = 1))
}


# initial conditions system of equation interval 4
#' @keywords internal
get_initial_conditions4 <- function(status,
                                    solution,
                                    parameter,
                                    trait_mainland_ancestor,
                                    num_observed_states,
                                    num_hidden_states) {
  n <- num_observed_states * num_hidden_states
  num_unique_states <- n

  if (status == 2 || status == 3 || status == 4) {
    if (trait_mainland_ancestor == "FALSE") {
      initial_conditions4 <- c(rep( sum((parameter[[3]]/n) * (solution[2,][(n + 1):(n + n)])), n), ### DM1: select DM2 in solution2
                               solution[2,][(n + n + n + 1):(n + n + n + n)],         ### E: select E in solution2
                                    sum((parameter[[3]]/n) * (solution[2,][(n + 1):(n + n)])))          ### DA1: select DM2 in solution2
    }
    #if the trait state of the species at the stem is known
    else if (trait_mainland_ancestor == trait_mainland_ancestor) {
      pos <- c((num_hidden_states*trait_mainland_ancestor + 1),
               num_hidden_states + trait_mainland_ancestor * num_hidden_states)
      initial_conditions4 <- c(rep(sum(parameter[[3]][pos] * (solution[2, ][(n + 1):(n + n)])[pos]) / num_hidden_states, n ), ### DM1: select DM2 in solution2
                               solution[2,][(n + n + n + 1):(n + n + n + n)],                                                                       ### E: select E in solution2
                               sum(parameter[[3]][pos] * (solution[2,][(n + 1):(n + n)])[pos]/num_hidden_states))          ### DA1: select DM2 in solution2
    }
  } else if (status == 1 || status == 5 || status == 6) {
    initial_conditions4 <- c(solution[2,][(n + 1):(n + n)],                                 ### DM1: select DM1 in solution1
                             solution[2,][(n + n + n + n + 1):(n + n + n + n + n)],         ### E: select E in solution1
                             solution[2,][length(solution[2, ]) - 1])                        ### DA1: select DA2 in solution1

  } else if (status == 8 || status == 9) {
    initial_conditions4 <- c(solution[2,][(n + 1):(n + n)],                                 ### DM1: select DM2 in solution3
                             solution[2,][(n + n + n + n + 1):(n + n + n + n + n)],         ### E: select E in solution3
                             solution[2,][length(solution[2, ]) - 1])                        ### DA1: select DA2 in solution3

  }
  return(matrix(initial_conditions4, nrow = 1))
}

