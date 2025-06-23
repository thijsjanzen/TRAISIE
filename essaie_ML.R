#load("C:/Users/P306709/Desktop/data_chapter3/datalist_test.RData")

#datalist <- datalist_test[1:11]

load("/Users/thijsjanzen/Downloads/data_list1.RData")
#load("C:/Users/P306709/Desktop/data_chapter3/data_list1.RData")

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



idparslist[[1]] <- c(1, 1, 2, 2) # lambda_c
idparslist[[2]] <- c(3, 3, 4, 4) # mu
idparslist[[3]] <- c(5, 5, 6, 6) # gamma
idparslist[[4]] <- c(7, 7, 8, 8) # lambda_a

idparslist[[5]] <- matrix(0, 4, 4)
# observed state transitions
idparslist[[5]][1, 3] <- 9 # 0 -> 1
idparslist[[5]][2, 4] <- 9 # 0 -> 1
idparslist[[5]][3, 1] <- 10 # 1 -> 0
idparslist[[5]][4, 2] <- 10 # 1 -> 0
colnames(idparslist[[5]]) <- c("0A", "0B", "1A", "1B")
rownames(idparslist[[5]]) <- colnames(idparslist[[5]])
# hidden state transitions
idparslist[[5]][1, 2] <- 11 # A -> B
idparslist[[5]][2, 1] <- 12 # B -> A
idparslist[[5]][3, 4] <- 11 # A -> B
idparslist[[5]][4, 3] <- 12 # B -> A

idparslist[[6]] <- 0 #p

idparsopt <- sort(unique(unlist(idparslist)))
idparsopt <- idparsopt[idparsopt > 0]

initvals <- c(0.546591, 0.678781, 0.9326754, 0.8583, 0.009326754, 0.008583, runif(6))

# first, let's test we have an init ll
initparsopt <- initvals
parsfix <- c(0)
trparsopt <- initparsopt / (1 + initparsopt)
trparsopt[which(initparsopt == Inf)] <- 1
trparsfix <- parsfix / (1 + parsfix)
trparsfix[which(parsfix == Inf)] <- 1



 parameter <- list(
   c(2.546591, 1.2, 1, 0.2),  # lambdac
   c(2.678781, 2, 1.9, 3),    # mu
   c(0.009326754, 0.003, 0.002, 0.2), #gamma
   c(1.008583, 1, 2, 1.5),  # lambda a
   matrix(c(                  # qmatrix
     0,    1,    0.5,  0,
     0,    0,    0.002, 0.005,
     runif(8)
   ), nrow = 4, ncol = 4, byrow = TRUE),
   0,                        # p
   c(1, 0)                   # trait mainland ancestor
 )


res <- DAISIE_DE_trait_loglik_CS (parameter,
                           datalist = data_list1[1:3],
                           methode = "lsodes",
                           rcpp_methode = "odeint::bulirsch_stoer",
                           atol = 1e-15,
                           rtol = 1e-15,
                           num_observed_states,
                           num_hidden_states,
                           cond = 1,
                           num_threads = 1,
                           verbose = FALSE,
                           use_Rcpp = 2)
if (1 == 2) {
initloglik <- treeLL::loglik_choosepar(trparsopt = trparsopt,
                                       trparsfix = trparsfix,
                                       idparsopt = idparsopt,
                                       idparsfix = c(0),
                                       idparslist = idparslist,
                                       datalist = data_list1,
                                       num_observed_states = num_observed_states,
                                       num_hidden_states = num_hidden_states,
                                       atol = 1e-15,
                                       rtol = 1e-15,
                                       methode = "ode45",
                                       rcpp_methode = "odeint::bulirsch_stoer",
                                       verbose = TRUE,
                                       use_Rcpp = 0,
                                       num_threads = 8)
initloglik

datalist1[1:10]


t0 <- Sys.time()
ml_res <- treeLL::calc_ml(  data_list1[1:40],
                            num_observed_states = 2,
                            num_hidden_states = 2,
                            idparslist = idparslist,
                            idparsopt = idparsopt,
                            initparsopt = initvals,
                            idparsfix = c(0),
                            parsfix = c(0),
                            atol = 1e-11,
                            rtol = 1e-11,
                            num_threads = 8,
                            verbose = TRUE,
                            use_Rcpp = 0)
t1 <- Sys.time()
ml_res2 <- treeLL::calc_ml(  data_list1,
                             num_observed_states = 2,
                             num_hidden_states = 2,
                             idparslist = idparslist,
                             idparsopt = idparsopt,
                             initparsopt = initvals,
                             idparsfix = c(0),
                             parsfix = c(0),
                             atol = 1e-11,
                             rtol = 1e-11,
                             num_threads = 8,
                             verbose = TRUE,
                             use_Rcpp = 2)
t2 <- Sys.time()
difftime(t1, t0)
difftime(t2, t1)
}