####
# Scenario 1. No-purchase is not observed 
#   --- by EM algorithm
####

Sys.setenv(LANG = "en")
setwd("~/Dropbox/RCode/Choice_Gibbs.git/src/M1.L1")


source(file="Metropolis-Hastings.R")



scenarioName <- "MNL.M1.L1_Scenario1.EM"

## Load simulated choice data (NEED TO RUN MNL.M1.L1_InitData.R TO GENERATE THE DATA FIRST)
load(file="MNL.M1.L1_InitData.RData")

# final observation consists of the Demand
observation1 <- list(demand=Demand)



# negative log-likelihood
negLogLikelihood.beta <- function(beta, data) {
    
    beta.coef <- beta[1:L]
    beta.const <- beta[(L+1):(L+M)]
    
    score <- rbind(t(exp(beta.const + X_Mat %*% beta.coef)), 1)
    choice.prob <- apply(score, 2, function(x) x/sum(x))
    
    logLikelihood <- data*log(choice.prob)
    
    return(-sum(logLikelihood))
    
}


prob.nopurchase.given.demand <- function(nopurchase, demand, lambda, beta, k) {
    
    if (any(nopurchase<0))
        return(0)
    else {
        beta.coef <- beta[1:L]
        beta.const <- beta[(L+1):(L+M)]
        
        score <- rbind(exp(beta.const + X_Mat[k,] %*% beta.coef), 1)
        choice.prob <- score/sum(score)
        
        
        return(dmultinom(x=c(demand, nopurchase), prob=choice.prob) * dpois(demand+nopurchase, lambda))
        
    }
}


sample = function(data, parameters, nrun=100) {

    demand <- data$demand
    
    # initialize arrays to save samples
    nopurchases <- array(0, dim=c(nrun, K))
    lambdas <- array(0, dim=c(nrun, 1))
    betas <- array(0, dim=c(nrun, L+M))
    
    # initial parameters
    lambda1 <- parameters$lambda
    beta1 <- parameters$beta
    
    
    # start EM loop
    for(i in 1:nrun) {

        # Expectation Step for each nopurchase in period 1, ..., K
        nopurchase2 <- rep(0,K)
        np.max <- 200 # upper limit used in integration

        for (j in 1:K) {
            
            # conditional distribution of nopurchase^k
            prob <- rep(0, np.max+1)
            for (np in 0:np.max) {
                prob[np+1] <- prob.nopurchase.given.demand(nopurchase=np, demand=demand[j], lambda=lambda1, beta=beta1, k=j) 
            }
            prob <- prob / sum(prob)
            
            # calculate expectation of nopurchase^k
            nopurchase2[j] = (0:np.max) %*% prob
        }
        
        
        
        # MLE of lambda2 (this is simple because of Poisson distribution)
        lambda2 <- (sum(demand, nopurchase2)) / K
        
        

        # MLE of beta2 (this need to be searched over the likelihood function)
        opt <- optim(par=c(-1,1), fn=negLogLikelihood.beta, data=rbind(demand, nopurchase2), method="BFGS")
        beta2 <- opt$par
        
        
        
        # save the samples obtained in the current iteration
        lambdas[i,] <- lambda2
        betas[i,] <- beta2
        nopurchases[i,] <- nopurchase2
        
        cat("Run:", i, "\tlambda:", lambda2, "\tbeta2:", beta2, "\tnopurchase[1]:", nopurchase2[1], "\tloglikelihood:", -opt$value, "\toptim.converge?", opt$convergence, "\n", sep=" ")
        
        lambda1 <- lambda2
        beta1 <- beta2
    }
    
    # results
    return(list(lambdas=lambdas, betas=betas, nopurchases=nopurchases))

} # end function 'sample'




## initial sampling input
param0 <- list(beta=c(-1, 1), lambda=40)
z1.EM <- sample(data=observation1, parameters=param0, nrun=500)

save(z1.EM, observation1, file=paste0(scenarioName, ".RData"))



