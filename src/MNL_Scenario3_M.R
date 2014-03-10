####
# Scenario 3. A noisy observation of No-purchase is observed 
# --- Noisy observation (say, traffic flow) is modeled as T = exp^epsilon1 * N + exp^epsilon2, where N is total realized number of potential customers and epsilon1 and epsilon2 are noise terms ~ N(epsilon1.mean, epsilon1.sd) and N(epsilon2.mean, epsilon2.sd)
# --- Case 1: only multiplicative noise (only epsilon1, epsilon2=-Inf)
####

Sys.setenv(LANG = "en")
setwd("~/Dropbox/RCode/Choice_Gibbs.git/src")

require("MCMCpack")
require("mvtnorm")


## Load simulated choice data (NEED TO RUN MNL_InitData.R TO GENERATE THE DATA FIRST)
load(file="MNL_InitData.RData")




# parameters for the noise terms epsilon1 and epsilon2
epsilon1.mean <- 3
epsilon1.sd <- 0.5
epsilon2 <- -Inf



### simulate traffic observation
epsilon1 <- rnorm(K, mean=epsilon1.mean, sd=epsilon1.sd)
traffic <- exp(epsilon1) * colSums(choice.mat) + exp(epsilon2)

# final observation consists of the sales of alternative 1 and 2 and the traffic flow
observation3 <- list(sales = choice.mat[1:M-1,], traffic=traffic)



#log-posterior of beta, to be called by M-H algorithm
logpost.beta <- function(beta, data) {

    if (any(beta<=0))
        return(-1.0e99)
    else {
        score <- apply(X_Mat, 1, function (x) exp(matrix(c(x,0,0), nrow=M, ncol=L, byrow=TRUE) %*% beta))
        choice.prob <- apply(score, 2, function(x) x/sum(x)) # M*N matrix
        
        logLikelihood <- data*log(choice.prob)
        
        logprior <- dmvnorm(log(beta), mean=beta.mu, sigma=beta.sg, log=TRUE)
            
        return(sum(logLikelihood) + logprior)
    }
}


logpost.d0 <- function(d0, data, lambda, beta, k, eps1.mu, eps1.sd) {
    
    score <-  exp(matrix(c(X_Mat[k,],0,0), nrow=M, ncol=L, byrow=TRUE) %*% beta)
    choice.prob <- score/sum(score) # M*N matrix
    
    dd <- c(data$sales,d0)
    nn <- sum(dd)

    return(dmultinom(x=dd, prob=choice.prob, log=TRUE) + dpois(nn, lambda, log=TRUE) + dnorm(x=log(data$traffic/nn), mean=eps1.mu, sd=eps1.sd, log=TRUE))
}

#implementation of discrete Metropolic-Hastings with Negative Binomial proposal
discreteMH <-function (logpost, proposal, start, m, ...) 
{
    pb = length(start)
    Mpar = array(0, c(m, pb))
    b = matrix(t(start))
    lb = logpost(start, ...)
    
    size = proposal$size
    
    accept = 0
    for (i in 1:m) {
        #proposed next point is Negative Binomial using the current position as its mean, variance is controled by $size$
        bc <- rnbinom(pb, size=size, mu=b)
        
        lbc = logpost(t(bc), ...) 
        #since proposal is asymmetric, need to adjust the jumping probability accordingly
        prob = exp(lbc + dnbinom(b, size=size, mu=bc, log=TRUE) - lb - dnbinom(bc, size=size, mu=b, log=TRUE) )
        
        if (is.na(prob) == FALSE) {
            if (runif(1) < prob) {
                lb = lbc
                b = bc
                accept = accept + 1
            }
        }
        Mpar[i, ] = b
    }
    accept = accept/m
    stuff = list(par = Mpar, accept = accept)
    return(stuff)
}




sample = function(data, parameters, nrun=1000) {
    
    sales <- data$sales
    traffic <- data$traffic

    
    # initialize arrays to save samples
    lambdas = array(0, dim=c(nrun, 1))
    d0s = array(0, dim=c(nrun, K))
    betas = array(0, dim=c(nrun, L))
    eps1.mus = array(0, dim=c(nrun, 1))
    eps1.sds = array(0, dim=c(nrun, 1))
    
    # initial parameters
    #lambda1 <- parameters$lambda
    d01 <- parameters$d0
    beta1 <- parameters$beta
    #eps1.mu1 <- parameters$eps1.mu
    eps1.sd1 <- parameters$eps1.sd
    
    
    # start sampling loop
    for(i in 1:nrun) {

        #update posterior of lambda by conjugacy
        alpha2 <- lambda.alpha + sum(sales) + sum(d01)
        beta2 <- lambda.beta + K
        
        #simulate lambda2
        lambda2 <- rpois(1, alpha2/beta2)
        
        
        
        #simulate beta2 by Metropolis-Hastings
        beta2 <- MCMCmetrop1R(logpost.beta, theta.init=beta1,
                         data=rbind(sales,d01),
                         thin=1, mcmc=1, burnin=500, tune=0.014,
                         verbose=0,  V=matrix(c(1,0,0,1),2,2))[1,]
                         #optim.lower=1e-6, optim.method="L-BFGS-B")[1,]


        #update and simulate eps1.mu and eps1.sd
        eps1 <- log(traffic) - log(colSums(sales)+d01)
        eps1.mu2 <- rnorm(1, mean=mean(eps1), sd=eps1.sd1/sqrt(K))
        eps1.pr2 <- rgamma(1, shape=eps1.pr.alpha0+K/2, rate=eps1.pr.beta0+sum((eps1-eps1.mu2)^2)/2)
        eps1.sd2 <- 1/sqrt(eps1.pr2)
        
        
        
        #simulate d0 by discrete Metropolis-Hastings
        d02 <- rep(0, K)
        d0.accept <- rep(0, K)
        for (j in 1:K) {
            dataj <- list(sales=sales[,j], traffic=traffic[j])
            sim <- discreteMH(logpost.d0, proposal=list(size=2), start=d01[j], m=500, 
                       data=dataj, lambda=lambda2, beta=beta2, k=j, eps1.mu=eps1.mu2, eps1.sd=eps1.sd2)
            d02[j] <- sim$par[500]
            d0.accept[j] <- sim$accept
        }
        cat("discreteMH acceptance rate was ", mean(d0.accept), "\n\n")
        

        # save the samples obtained in the current iteration
        lambdas[i,] = lambda2
        d0s[i,] = d02
        betas[i,] = beta2
        eps1.mus[i,] = eps1.mu2
        eps1.sds[i,] = eps1.sd2

        cat("Run:", i, "\tlambda:", lambda2, "\tbeta2:", beta2, "\teps1.mu2:", eps1.mu2, "\teps1.sd2:", eps1.sd2, "\n", sep=" ")
        
        lambda1 <- lambda2
        d01 <- d02
        beta1 <- beta2
        eps1.mu1 <- eps1.mu2
        eps1.sd1 <- eps1.sd2
    }
    
    # results
    return(list(lambdas=lambdas, d0s=d0s, betas=betas, eps1.mus=eps1.mus, eps1.sds=eps1.sds))

} # end function 'sample'




### initialize input before sampling
# beta prior ~ logN(beta.mu, beta.sg)
beta.mu <- c(-2.5, -2.5)
beta.sg <- matrix(c(0.5, 0, 0, 0.5), 2, 2)

# lambda ~ Gamma(lambda.alpha, lambda.beta)
lambda.alpha <- 0.5
lambda.beta <- 0.01

# prior of the precision of epsilon
eps1.pr.alpha0 <- 0.2
eps1.pr.beta0 <- 0.1

## initial sampling input
param0 <- list(beta=rep(0.1, L), d0=rep(10,K), eps1.sd=eps1.pr.alpha0/eps1.pr.beta0)
nrun <- 5000



### sample
z <- sample(data=observation3, parameters=param0, nrun=nrun)

save.image(file="MNL_Scenario3_M.RData")

#stopCluster(cl)



### Visualize results
burnin <- 0.1*nrun

#plot lambda
samples.lambda <- z$lambdas
plot(samples.lambda, type="l")


samples.lambda.truncated <- samples.lambda[burnin:nrun,]
quantile(samples.lambda.truncated, c(.025,.5,.975))
mean(samples.lambda.truncated)
hist(samples.lambda.truncated)



#plot beta
samples.beta <- data.frame(z$betas)
plot(samples.beta$X1, type="l")
plot(samples.beta$X2, type="l")

plot(c(1,nrun), c(-0.05, 0.2), type="n", xlab="Samples", ylab="a", xaxt="n", yaxt="n")
lines(samples.beta$X1, col="GREEN")
lines(samples.beta$X2, col="BLUE")
axis(side=1)
axis(side=2)

samples.beta.truncated <- samples.beta[burnin:nrun,]
quantile(samples.beta.truncated$X1, c(.025,.5,.975))
quantile(samples.beta.truncated$X2, c(.025,.5,.975))
colMeans(samples.beta.truncated)
hist(samples.beta.truncated$X1)
hist(samples.beta.truncated$X2)



#plot d0[1]
samples.d0 <- data.frame(z$d0s)
plot(samples.d0$X1, type="l")
plot(samples.d0$X2, type="l")
plot(samples.d0$X3, type="l")
plot(samples.d0$X4, type="l")
plot(samples.d0$X5, type="l")


samples.d0.truncated <- samples.d0[burnin:nrun,]
quantile(samples.d0.truncated$X1, c(.025,.5,.975))
quantile(samples.d0.truncated$X2, c(.025,.5,.975))
quantile(samples.d0.truncated$X3, c(.025,.5,.975))
quantile(samples.d0.truncated$X4, c(.025,.5,.975))
quantile(samples.d0.truncated$X5, c(.025,.5,.975))
colMeans(samples.d0.truncated)
hist(samples.d0.truncated$X1)
hist(samples.d0.truncated$X2)
hist(samples.d0.truncated$X3)
hist(samples.d0.truncated$X4)
hist(samples.d0.truncated$X5)



#plot eps1.mu and eps1.sd
samples.eps1 <- data.frame(mu=z$eps1.mus, sd=z$eps1.sds)
plot(samples.eps1$mu, type="l")
plot(samples.eps1$sd, type="l")

samples.eps1.truncated <- samples.eps1[burnin:nrun,]
quantile(samples.eps1.truncated$mu, c(.025,.5,.975))
quantile(samples.eps1.truncated$sd, c(.025,.5,.975))
colMeans(samples.eps1.truncated)
hist(samples.eps1.truncated$mu)
hist(samples.eps1.truncated$sd)

