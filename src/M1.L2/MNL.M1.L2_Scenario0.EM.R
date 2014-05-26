####
# Scenario 0. All choices are observed 
#   --- by EM algorithm
####

Sys.setenv(LANG = "en")
setwd("~/Dropbox/RCode/Choice_Gibbs.git/src/M1.L2")


source(file="Metropolis-Hastings.R")



scenarioName <- "MNL.M1.L2_Scenario0.EM"

## Load simulated choice data (NEED TO RUN MNL.M1.L2_InitData.R TO GENERATE THE DATA FIRST)
load(file="MNL.M1.L2_InitData.RData")

# final observation consists of Demand, and NoPurchase
observation0 <- list(demand=Demand, nopurchase=NoPurchase)



# negative log-likelihood
negLogLikelihood.beta <- function(beta, data) {
    
    beta.coef <- beta[1:L]
    beta.const <- beta[(L+1):(L+M)]
    
    score <- rbind(t(exp(beta.const + X_Mat %*% beta.coef)), 1)
    choice.prob <- apply(score, 2, function(x) x/sum(x))
    
    logLikelihood <- data*log(choice.prob)
    
    return(-sum(logLikelihood))
    
}


sample = function(data, parameters, nrun=100) {

    demand <- data$demand
    nopurchase <- data$nopurchase
    
    # initialize arrays to save samples
    lambdas <- array(0, dim=c(nrun, 1))
    betas <- array(0, dim=c(nrun, L+M))
    
    # initial parameters
    lambda <- parameters$lambda
    beta <- parameters$beta
    
    
    # start EM loop
    for(i in 1:nrun) {

        # Expectation Step for each nopurchase in period 1, ..., K
        # nothing to E here ...     
        
        
        
        # MLE of lambda (this is simple because of Poisson distribution)
        lambda <- (sum(demand, nopurchase)) / K
        
        
        
        # MLE of beta (this need to be searched over the likelihood function)
        opt <- optim(par=c(-1,1,1), fn=negLogLikelihood.beta, data=rbind(demand, nopurchase), method="BFGS")
        beta <- opt$par
        
        
        
        # save the samples obtained in the current iteration
        lambdas[i,] <- lambda
        betas[i,] <- beta
        
        cat("Run:", i, "\tlambda:", lambda, "\tbeta:", beta, "\tloglikelihood:", -opt$value, "\toptim.converge?", opt$convergence, "\n", sep=" ")
        
    }
    
    # results
    return(list(lambdas=lambdas, betas=betas))

} # end function 'sample'




## initial sampling input
param0 <- list(beta=c(-1, 1, 1), lambda=30)

z0.EM <- sample(data=observation0, parameters=param0, nrun=5)

save(z0.EM, observation0, file=paste0(scenarioName, ".RData"))



