#' @keywords internal
solve_branch_cpp <- function(chosen_func,
                             initial_conditions,
                             time,
                             parameter,
                             methode = "odeint::bulirsch_stoer",
                             atol = 1e-15,
                             rtol = 1e-15) {

  lambda_c <- parameter[[1]]
  mus      <- parameter[[2]]
  gammas   <- parameter[[3]]
  lambda_a <- parameter[[4]]
  q_matrix       <- parameter[[5]]
  p_value       <- parameter[[6]]



  solution <- cpp_solve(lambda_c,
                        lambda_a,
                        mus,
                        gammas,
                        q_matrix,
                        p_value,
                        chosen_func,
                        methode,
                        initial_conditions,
                        time,
                        atol,
                        rtol)

  res <- matrix(data = NA, nrow = 2, ncol = length(solution$states))
  res[1, ] <- initial_conditions
  res[2, ] <- solution$states

  return(res)
}
