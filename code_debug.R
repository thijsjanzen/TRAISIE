load("/Users/thijsjanzen/Documents/GitHub/treeLL/height_data_list1.RData")

idparslist <- list()
idparslist[[1]] <- c(1) # lambda_c
idparslist[[2]] <- c(2) # mu
idparslist[[3]] <- c(3) # gamma
idparslist[[4]] <- c(4) # lambda_a

idparslist[[5]] <- matrix(0, 1, 1)
colnames(idparslist[[5]]) <- c("0")
rownames(idparslist[[5]]) <- colnames(idparslist[[5]])
# hidden state transitions
idparslist[[5]][1, 1] <- 5 #

idparslist[[6]] <- c(6)#p

idparsopt <- c(1:4)
idparsopt <- idparsopt[idparsopt > 0]
initvals2 <- c(1, 1, 0.005, 0.005)



#### example of starting values.
#### better to try different starting values otherwise risk of falling on local minima instead of global




# first, let's test we have an init ll
initparsopt <- initvals2
idparsfix <- c(5, 6)
parsfix = c(0.000, 0)
trparsopt <- initparsopt / (1 + initparsopt)
trparsopt[which(initparsopt == Inf)] <- 1
trparsfix <- parsfix / (1 + parsfix)
trparsfix[which(parsfix == Inf)] <- 1



ml_results11 <- treeLL:::calc_ml(data_list1,
                                 num_observed_states = 1,
                                 num_hidden_states = 1,
                                 idparslist = idparslist,
                                 maxiter = 2000 * round((1.25) ^ length(idparsopt)),
                                 idparsopt = idparsopt,
                                 initparsopt = initvals2,
                                 idparsfix = idparsfix,
                                 parsfix = parsfix,
                                 atol = 1e-15,
                                 rtol = 1e-15,
                                 num_threads = 8,
                                 verbose = TRUE,
                                 use_Rcpp = 2)
