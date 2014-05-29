####
# Scenario 1c. No-purchase is not observed 
#   --- but demand is censored by inventory
#   --- by EM algorithm
####

Sys.setenv(LANG = "en")
setwd("~/Dropbox/RCode/Choice_Gibbs.git/src/M1.L2")



scenarioName <- "MNL.M1.L2_Scenario1c.EM"

## Load simulated choice data (NEED TO RUN MNL.M1.L2_InitData.R TO GENERATE THE DATA FIRST)
load(file="MNL.M1.L2_InitData.RData")

# final observation consists of the Sales and Stockout status
observation1c <- list(sales=Sales, stockout=Stockout)



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


prob.nopurchase_demand.given.sales <- function(nopurchase, demand, sales, lambda, beta, k) {
    
    if (any(nopurchase<0) | any(demand<sales))
        return(0)
    else {
        beta.coef <- beta[1:L]
        beta.const <- beta[(L+1):(L+M)]
        
        score <- rbind(exp(beta.const + X_Mat[k,] %*% beta.coef), 1)
        choice.prob <- score/sum(score)
        
        dd <- c(demand, nopurchase)
        nn <- sum(dd)
        
        return(dmultinom(x=dd, prob=choice.prob) * dpois(nn, lambda))
    }
}


sample = function(data, parameters, nrun=100) {

    sales <- data$sales
    stockout <- data$stockout
    
    # initialize arrays to save samples
    nopurchases <- array(0, dim=c(nrun, K))
    demands <- array(0, dim=c(nrun, K))
    lambdas <- array(0, dim=c(nrun, 1))
    betas <- array(0, dim=c(nrun, L+M))
    
    # initial parameters
    lambda <- parameters$lambda
    beta <- parameters$beta
    demand <- sales
    
    
    # start EM loop
    for(i in 1:nrun) {
        
        # Expectation Step for each nopurchase/demand in period 1, ..., K
        nopurchase <- rep(0,K)
        
        for (j in 1:K) {
            if (!stockout[j]) {
                # conditional distribution of nopurchase^k
                prob <- rep(0, np.max+1)
                for (np in 0:np.max) {
                    prob[np+1] <- prob.nopurchase.given.demand(nopurchase=np, demand=demand[j], lambda=lambda, beta=beta, k=j) 
                }
                prob <- prob / sum(prob)
                
                # calculate expectation of nopurchase^k
                nopurchase[j] <- (0:np.max) %*% prob
            } else {
                # conditional distribution of nopurchase^k
                prob <- matrix(0, np.max+1, d.max+1)
                for (np in 0:np.max) {
                    for (d in 0:d.max) {
                        prob[np+1, d+1] <- prob.nopurchase_demand.given.sales(nopurchase=np, demand=d, sales=sales[j], lambda=lambda, beta=beta, k=j) 
                    }
                }
                prob <- prob / sum(prob)
                
                # calculate expectation of nopurchase^k
                nopurchase[j] <- (0:np.max) %*% rowSums(prob)
                demand[j] <- (0:d.max) %*% colSums(prob)
            }
        }
        
        
        
        # MLE of lambda (this is simple because of Poisson distribution)
        lambda <- (sum(demand, nopurchase)) / K
        
        
        
        # MLE of beta (this need to be searched over the likelihood function)
        opt <- optim(par=c(-1,1,1), fn=negLogLikelihood.beta, data=rbind(demand, nopurchase), method="BFGS")
        beta <- opt$par
        
        
        
        # save the samples obtained in the current iteration
        nopurchases[i,] <- nopurchase
        demands[i,] <- demand
        lambdas[i,] <- lambda
        betas[i,] <- beta
        
        cat("Run:", i, "\tlambda:", lambda, "\tbeta:", beta, "\tnopurchase[1]:", nopurchase[1], "\tdemand[1]:", demand[1], "\tloglikelihood:", -opt$value, "\toptim.converge?", opt$convergence, "\n", sep=" ")
        
    }
    
    # results
    return(list(lambdas=lambdas, betas=betas, nopurchases=nopurchases, demands=demands))

} # end function 'sample'




## initial sampling input
np.max <- 100 # upper limit used in integration
d.max <- 100 # upper limit used in integration
param0 <- list(beta=c(-1, 1, 1), lambda=30)

z1c.EM <- sample(data=observation1c, parameters=param0, nrun=5000)

save(z1c.EM, observation1c, file=paste0(scenarioName, ".RData"))



