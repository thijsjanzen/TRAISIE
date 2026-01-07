require(treeLL)
data("Galapagos_datalist", package = "DAISIE")
datalist <- Galapagos_datalist
for (i in 2:length(datalist)) {
  datalist[[i]]$traits <- sample(c(0, 1), size = length(datalist[[i]]$branching_times),
                                 replace = TRUE)
  datalist[[i]]$sampling_fraction <- rep(1, 2)
  datalist[[i]]$phylogeny <- DDD::brts2phylo(datalist[[i]]$branching_times[-c(1, 2)])
  datalist[[i]]$root_state <- c(1, 1)
}

parameter <- list(2.546591, 2.678781, 0.009326754, 1.008583, matrix(c(0), nrow = 1), 0, 0 )

res1 <-  DAISIE_DE_trait_logp0(
  datalist            = datalist,
  weight_method = "mainland_stationary_weights",
  parameter           = parameter,
  num_observed_states     = 1,
  num_hidden_states       = 1,
  atol                = 1e-10,
  rtol                = 1e-10,
  methode             = "lsodes"
)

parameter <- list()
parameter[[1]] <- c(2.546591, 2.546591)
parameter[[2]] <- c(2.678781, 2.678781)
parameter[[3]] <- c(0.009326754, 0.009326754)
parameter[[4]] <- c(1.008583, 1.008583)
parameter[[5]] <- matrix(c(0), nrow = 2, ncol = 2)
parameter[[6]] <- 0


res2 <- DAISIE_DE_trait_loglik_CS(parameter = parameter,
                                  weight_method = "mainland_stationary_weights",
                                  datalist = datalist,
                                  num_observed_states = 2,
                                  num_hidden_states = 1,
                                  use_Rcpp = 2)



# but now we want to optimize

idparslist <- list()
num_observed_states <- 2
num_hidden_states <- 2

# we need in the idparslist the following
# lambda_c
# mu
# gamma
# lambda_a
# q

# let's first do a Constant Rates model
# 0A, 0B, 1A, 1B
# c(1, 1, 1, 1) # CR
# c(1, 1, 2, 2) # ETD
# c(1, 2, 1, 2) # CTD
# c(1, 2, 3, 4) # fancy



idparslist[[1]] <- c(1, 1, 1, 1) # lambda_c
idparslist[[2]] <- c(2, 2, 2, 2) # mu
idparslist[[3]] <- c(3, 3, 3, 3) # gamma
idparslist[[4]] <- c(4, 4, 4, 4) # lambda_a

idparslist[[5]] <- matrix(0, 4, 4)
# observed state transitions
idparslist[[5]][1, 3] <- 5 # 0 -> 1
idparslist[[5]][2, 4] <- 5 # 0 -> 1
idparslist[[5]][3, 1] <- 6 # 1 -> 0
idparslist[[5]][4, 2] <- 6 # 1 -> 0
colnames(idparslist[[5]]) <- c("0A", "0B", "1A", "1B")
rownames(idparslist[[5]]) <- colnames(idparslist[[5]])
# hidden state transitions
idparslist[[5]][1, 2] <- 7 # A -> B
idparslist[[5]][2, 1] <- 8 # B -> A
idparslist[[5]][3, 4] <- 7 # A -> B
idparslist[[5]][4, 3] <- 8 # B -> A

idparslist[[6]] <- 0 #p

idparsopt <- sort(unique(unlist(idparslist)))
idparsopt <- idparsopt[idparsopt > 0]

initvals <- c(2.546591, 2.678781, 0.009326754, 1.008583, runif(4))

# first, let's test we have an init ll
initparsopt <- initvals
parsfix <- c(0)
trparsopt <- initparsopt / (1 + initparsopt)
trparsopt[which(initparsopt == Inf)] <- 1
trparsfix <- parsfix / (1 + parsfix)
trparsfix[which(parsfix == Inf)] <- 1

initloglik <- treeLL::loglik_choosepar(trparsopt = trparsopt,
                                       weight_method = "mainland_stationary_weights",
                               trparsfix = trparsfix,
                               idparsopt = idparsopt,
                               idparsfix = c(0),
                               idparslist = idparslist,
                               datalist = datalist,
                               num_observed_states = num_observed_states,
                               num_hidden_states = num_hidden_states,
                               cond = 1,
                               atol = 1e-9,
                               rtol = 1e-9,
                               methode = "ode45",
                               verbose = TRUE,
                               use_Rcpp = FALSE)


initloglik

t0 <- Sys.time()
ml_res <- treeLL::calc_ml(  datalist,
                            weight_method = "mainland_stationary_weights",
                            num_observed_states = 2,
                            num_hidden_states = 2,
                            idparslist = idparslist,
                            idparsopt = idparsopt,
                            initparsopt = initvals,
                            idparsfix = c(0),
                            parsfix = c(0),
                            verbose = FALSE,
                            atol = 1e-15,
                            rtol = 1e-15,
                            use_Rcpp = TRUE)
t1 <- Sys.time()
ml_res2 <- treeLL::calc_ml( datalist,
                            weight_method = "mainland_stationary_weights",
                            num_observed_states = 2,
                            num_hidden_states = 2,
                            idparslist = idparslist,
                            idparsopt = idparsopt,
                            initparsopt = initvals,
                            idparsfix = c(0),
                            parsfix = c(0),
                            atol = 1e-15,
                            rtol = 1e-15,
                            verbose = FALSE,
                            use_Rcpp = FALSE)
t2 <- Sys.time()
difftime(t1, t0)
difftime(t2, t1)


