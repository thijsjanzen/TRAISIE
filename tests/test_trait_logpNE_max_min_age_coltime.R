test_that("logpNE_max_min_age_coltime", {



    brts <- c(4, 3, 2.5)

    parameter <- list(2.546591, 2.678781, 0.009326754, 1.008583, matrix(c(0), nrow = 1), 0 )

    res1 <-  DAISIE_DE_trait_logpNE_max_min_age_hidden(brts                  = brts,
                                                       trait                 = 0,
                                                       status                = 8,
                                                       parameter             = parameter,
                                                       num_observed_states   = 1,
                                                       num_hidden_states     = 1,
                                                       atol                  = 1e-15,
                                                       rtol                  = 1e-15,
                                                       methode               = "ode45"
                                                     )

    pars1 <- c(2.546591, 2.678781, Inf, 0.009326754, 1.008583)
    res2  <-  DAISIE:::DAISIE_loglik_CS_choice(pars1 = pars1,
                                               pars2 = c(100, 11, 0, 2),
                                               brts = brts,
                                               stac = 8,
                                               missnumspec = 0,
                                               datalist = datalist)

    testthat::expect_equal(res1, res2, tolerance = 0.01)
  }
)
