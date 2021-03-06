####
# Scenario 3c. A noisy observation of No-purchase is observed 
# --- Noisy observation (say, traffic flow) is modeled as T = exp^epsilon1 * N + exp^epsilon2, where N is total realized number of potential customers and epsilon1 and epsilon2 are noise terms ~ N(epsilon1.mean, epsilon1.sd) and N(epsilon2.mean, epsilon2.sd)
# --- Case 1: only multiplicative noise (only epsilon1, epsilon2=-Inf)
#   --- but demand is censored by inventory
#   --- by EM algorithm
####

Sys.setenv(LANG = "en")
setwd("~/Dropbox/RCode/Choice_Gibbs.git/src/M1.L1")



scenarioName <- "MNL.M1.L1_Scenario3Mc.m.EM"

## Load simulated choice data (NEED TO RUN MNL.M1.L1_InitData.R TO GENERATE THE DATA FIRST)
load(file="MNL.M1.L1_InitData.RData")

# final observation consists of the Sales, Stockout status and the Traffic flow
observation3Mc.m <- list(sales=Sales, stockout=Stockout, traffic=TrafficM.m)



# negative log-likelihood
negLogLikelihood.beta <- function(beta, data) {
    
    beta.coef <- beta[1:L]
    beta.const <- beta[(L+1):(L+M)]
    
    score <- rbind(t(exp(beta.const + X_Mat %*% beta.coef)), 1)
    choice.prob <- apply(score, 2, function(x) x/sum(x))
    
    logLikelihood <- data*log(choice.prob)
    
    return(-sum(logLikelihood))
    
}


prob.nopurchase.given.demand <- function(nopurchase, demand, traffic, lambda, beta, k, eps1.mu, eps1.sd) {
    
    if (any(nopurchase<0))
        return(0)
    else {
        
        beta.coef <- beta[1:L]
        beta.const <- beta[(L+1):(L+M)]
        
        score <- rbind(exp(beta.const + X_Mat[k,] %*% beta.coef), 1)
        choice.prob <- score/sum(score)
    
        return(dmultinom(x=c(demand, nopurchase), prob=choice.prob) * dpois(demand+nopurchase, lambda) * dnorm(x=log(traffic/(demand+nopurchase)), mean=eps1.mu, sd=eps1.sd))
        
    }
}


prob.nopurchase_demand.given.sales <- function(nopurchase, demand, sales, traffic, lambda, beta, k, eps1.mu, eps1.sd) {
    
    if (any(nopurchase<0) | any(demand<sales))
        return(0)
    else {
        beta.coef <- beta[1:L]
        beta.const <- beta[(L+1):(L+M)]
        
        score <- rbind(exp(beta.const + X_Mat[k,] %*% beta.coef), 1)
        choice.prob <- score/sum(score)
        
        dd <- c(demand, nopurchase)
        nn <- sum(dd)
        
        return(dmultinom(x=dd, prob=choice.prob) * dpois(nn, lambda) * dnorm(x=log(traffic/(demand+nopurchase)), mean=eps1.mu, sd=eps1.sd))
    }
}


sample = function(data, parameters, nrun=100) {

    sales <- data$sales
    stockout <- data$stockout
    traffic <- data$traffic
    
    # initialize arrays to save samples
    nopurchases <- array(0, dim=c(nrun, K))
    demands <- array(0, dim=c(nrun, K))
    lambdas <- array(0, dim=c(nrun, 1))
    betas <- array(0, dim=c(nrun, L+M))
    eps1.sds <- array(0, dim=c(nrun, 1))
    
    # initial parameters
    lambda <- parameters$lambda
    beta <- parameters$beta
    eps1.mu <- parameters$eps1.mu
    eps1.sd <- parameters$eps1.sd
    demand <- sales
    
    # start EM loop
    for(i in 1:nrun) {

        # Expectation Step for each nopurchase in period 1, ..., K
        nopurchase <- rep(0,K)
        
        for (j in 1:K) {
            if (!stockout[j]) {
                # conditional distribution of nopurchase^k
                prob <- rep(0, np.max+1)
                for (np in 0:np.max) {
                    prob[np+1] <- prob.nopurchase.given.demand(nopurchase=np, demand=demand[j], traffic=traffic[j], lambda=lambda, beta=beta, k=j, eps1.mu=eps1.mu, eps1.sd=eps1.sd)
                }
                prob <- prob / sum(prob)
                
                # calculate expectation of nopurchase^k
                nopurchase[j] <- (0:np.max) %*% prob
            } else {
                # conditional distribution of nopurchase^k
                prob <- matrix(0, np.max+1, d.max+1)
                for (np in 0:np.max) {
                    for (d in 0:d.max) {
                        prob[np+1, d+1] <- prob.nopurchase_demand.given.sales(nopurchase=np, demand=d, sales=sales[j], traffic=traffic[j], lambda=lambda, beta=beta, k=j, eps1.mu=eps1.mu, eps1.sd=eps1.sd)
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
        opt <- optim(par=c(-1,1), fn=negLogLikelihood.beta, data=rbind(demand, nopurchase), method="BFGS")
        beta <- opt$par
        
        
        
        # MLE of eps1.sd (simply sample stdev)
        eps1 <- log(traffic) - log(demand + nopurchase)
        #eps1.mu <- mean(eps1)
        eps1.sd <- sqrt( sum((eps1 - eps1.mu)^2) / K )
        
        
        
        # save the samples obtained in the current iteration
        nopurchases[i,] <- nopurchase
        demands[i,] <- demand
        lambdas[i,] <- lambda
        betas[i,] <- beta
        eps1.sds[i,] <- eps1.sd
        
        cat("Run:", i, "\tlambda:", lambda, "\tbeta:", beta, "\tnopurchase[1]:", nopurchase[1], "\tdemand[1]:", demand[1], "\teps1.sd:", eps1.sd, "\teps1.mu:", eps1.mu, "\tloglikelihood:", -opt$value, "\toptim.converge?", opt$convergence, "\n", sep=" ")
        
    }
    
    # results
    return(list(lambdas=lambdas, betas=betas, eps1.sds=eps1.sds, nopurchases=nopurchases, demands=demands))

} # end function 'sample'




## initial sampling input
np.max <- 100 # upper limit used in integration
d.max <- 100 # upper limit used in integration
param0 <- list(beta=c(-1, 1), lambda=30, eps1.mu=epsilon1.mean, eps1.sd=1)

z3Mc.m.EM <- sample(data=observation3Mc.m, parameters=param0, nrun=100)

save(z3Mc.m.EM, observation3Mc.m, file=paste0(scenarioName, ".RData"))



