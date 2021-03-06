####
# Scenario 0c. All choices are observed 
#   --- but demand is censored by inventory
#   --- by EM algorithm
####

Sys.setenv(LANG = "en")
setwd("~/Dropbox/RCode/Choice_Gibbs.git/src/M1.L2")



scenarioName <- "MNL.M1.L2_Scenario0c.EM"

## Load simulated choice data (NEED TO RUN MNL.M1.L2_InitData.R TO GENERATE THE DATA FIRST)
load(file="MNL.M1.L2_InitData.RData")

# final observation consists of Sales, Stockout status, and Rest
observation0c <- list(sales=Sales, stockout=Stockout, rest=NoPurchase+LostSales)



# negative log-likelihood
negLogLikelihood.beta <- function(beta, data) {
    
    beta.coef <- beta[1:L]
    beta.const <- beta[(L+1):(L+M)]
    
    score <- rbind(t(exp(beta.const + X_Mat %*% beta.coef)), 1)
    choice.prob <- apply(score, 2, function(x) x/sum(x))
    
    logLikelihood <- data*log(choice.prob)
    
    return(-sum(logLikelihood))
    
}


prob.nopurchase.given.demand <- function(nopurchase, sales, rest, lambda, beta, k) {
    
    if (any(nopurchase<0) | any(nopurchase>rest))
        return(0)
    else {
        beta.coef <- beta[1:L]
        beta.const <- beta[(L+1):(L+M)]
        
        score <- rbind(exp(beta.const + X_Mat[k,] %*% beta.coef), 1)
        choice.prob <- score/sum(score)
        
        dd <- c(sales + rest - nopurchase, nopurchase)
        nn <- sum(dd)
        
        
        return(dmultinom(x=dd, prob=choice.prob) * dpois(nn, lambda))
    }
}


sample = function(data, parameters, nrun=100) {

    sales <- data$sales
    stockout <- data$stockout
    rest <- data$rest
    
    # initialize arrays to save samples
    nopurchases <- array(0, dim=c(nrun, K))
    lambdas <- array(0, dim=c(nrun, 1))
    betas <- array(0, dim=c(nrun, L+M))
    
    # initial parameters
    nopurchase <- rest
    lambda <- parameters$lambda
    beta <- parameters$beta
    
    
    # start EM loop
    for(i in 1:nrun) {

        # Expectation Step for each nopurchase in period 1, ..., K with stock-out
        for (j in which(stockout)) {
            
            # conditional distribution of nopurchase^k
            prob <- rep(0, np.max+1)
            for (np in 0:np.max) {
                prob[np+1] <- prob.nopurchase.given.demand(nopurchase=np, sales=sales[j], rest=rest[j], lambda=lambda, beta=beta, k=j) 
            }
            prob <- prob / sum(prob)
            
            # calculate expectation of nopurchase^k
            #nopurchase[j] <- round((0:np.max) %*% prob)
            nopurchase[j] <- (0:np.max) %*% prob
        }
        
        
        
        # MLE of lambda (this is simple because of Poisson distribution)
        lambda <- (sum(sales, rest)) / K
        
        
        
        # MLE of beta (this need to be searched over the likelihood function)
        demand <- sales + rest - nopurchase
        opt <- optim(par=c(-1,1,1), fn=negLogLikelihood.beta, data=rbind(demand, nopurchase), method="BFGS")
        beta <- opt$par
        
        
        
        # save the samples obtained in the current iteration
        nopurchases[i,] <- nopurchase
        lambdas[i,] <- lambda
        betas[i,] <- beta
        
        cat("Run:", i, "\tlambda:", lambda, "\tbeta:", beta, "\tnopurchase[1]:", nopurchase[1], "\tloglikelihood:", -opt$value, "\toptim.converge?", opt$convergence, "\n", sep=" ")
        
    }
    
    # results
    return(list(lambdas=lambdas, betas=betas, nopurchases=nopurchases))

} # end function 'sample'




## initial sampling input
np.max <- 100 # upper limit used in integration
param0 <- list(beta=c(-1, 1, 1), lambda=30)

z0c.EM <- sample(data=observation0c, parameters=param0, nrun=100)

save(z0c.EM, observation0c, file=paste0(scenarioName, ".RData"))



