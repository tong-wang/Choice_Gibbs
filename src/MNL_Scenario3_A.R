####
# Scenario 3. A noisy observation of No-purchase is observed 
# --- Noisy observation (say, traffic flow) is modeled as T = exp^epsilon1 * N + exp^epsilon2, where N is total realized number of potential customers and epsilon1 and epsilon2 are noise terms ~ N(epsilon1.mean, epsilon1.sd) and N(epsilon2.mean, epsilon2.sd)
# --- Case 2: only additive noise (only epsilon2, epsilon1=0)
####

Sys.setenv(LANG = "en")
setwd("~/Dropbox/RCode/Choice_Gibbs.git/src")

require("MCMCpack")
require("mvtnorm")



## Load simulated choice data (NEED TO RUN MNL_InitData.R TO GENERATE THE DATA FIRST)
load(file="MNL_InitData.RData")



# parameters for the noise terms epsilon1 and epsilon2
epsilon1 <- 0
epsilon2.mean <- 5
epsilon2.sd <- 1



### simulate traffic observation
epsilon2 <- rnorm(K, mean=epsilon2.mean, sd=epsilon2.sd)
traffic <- exp(epsilon1) * colSums(choice.mat) + exp(epsilon2)

# final observation consists of the sales of alternative 1 and 2 and the traffic flow
observation3 <- list(sales = choice.mat[1:M-1,], traffic=traffic)



#log-posterior of beta, to be called by M-H algorithm
logpost.beta <- function(beta, data) {

    if (any(beta<0))
        return(-Inf)
    else {
        score <- apply(X_Mat, 1, function (x) exp(matrix(c(x,0,0), nrow=M, ncol=L, byrow=TRUE) %*% beta))
        choice.prob <- apply(score, 2, function(x) x/sum(x)) # M*N matrix
        
        logLikelihood <- data*log(choice.prob)
        
        logprior <- dmvnorm(log(beta), mean=beta.mu, sigma=beta.sg, log=TRUE)
            
        return(sum(logLikelihood) + logprior)
    }
}


logpost.d0 <- function(d0, data, lambda, beta, k, eps2.mu, eps2.sd) {
    
    if (any(d0<0))
        return(-Inf)
    else {
        score <-  exp(matrix(c(X_Mat[k,],0,0), nrow=M, ncol=L, byrow=TRUE) %*% beta)
        choice.prob <- score/sum(score) # M*N matrix
        
        dd <- c(data$sales,d0)
        nn <- sum(dd)
    
        return(dmultinom(x=dd, prob=choice.prob, log=TRUE) + dpois(nn, lambda, log=TRUE) + dnorm(x=log(data$traffic-nn), mean=eps2.mu, sd=eps2.sd, log=TRUE))
    }
}

# implementation of discrete Metropolis-Hastings algorithm with discretized symmetric univariate Normal proposal
# tested for one- and high-dimensional discrete distribution
discreteMH.norm <-function (logpost, start, scale, nrun, ...) 
{
    dim = length(start)
    MC = array(0, c(nrun, dim))
    b1 = start
    ll.b1 = logpost(start, ...)
    
    
    accept = 0
    for (i in 1:nrun) {
        #proposed next point is discrete Multivariate Normal using the current position as its mean, variance is controled by $scale$
        b2 <- round(rnorm(dim, mean=b1, sd=scale))
        
        ll.b2 = logpost(b2, ...) 
        ll.ratio = exp(ll.b2 - ll.b1)
        
        if (!is.na(ll.ratio)) {
            if (runif(1) <= ll.ratio) {
                ll.b1 = ll.b2
                b1 = b2
                accept = accept + 1
            }
        }
        MC[i, ] = b1
        
    }
    accept = accept/nrun
    
    list(MC = MC, accept = accept)
}




sample = function(data, parameters, nrun=1000) {
    
    sales <- data$sales
    traffic <- data$traffic

    
    # initialize arrays to save samples
    lambdas = array(0, dim=c(nrun, 1))
    d0s = array(0, dim=c(nrun, K))
    betas = array(0, dim=c(nrun, L))
    #eps1.mus = array(0, dim=c(nrun, 1))
    #eps1.sds = array(0, dim=c(nrun, 1))
    eps2.mus = array(0, dim=c(nrun, 1))
    eps2.sds = array(0, dim=c(nrun, 1))
    
    # initial parameters
    #lambda1 <- parameters$lambda
    d01 <- parameters$d0
    beta1 <- parameters$beta
    #eps1.mu1 <- parameters$eps1.mu
    #eps1.sd1 <- parameters$eps1.sd
    #eps2.mu1 <- parameters$eps2.mu
    eps2.sd1 <- parameters$eps2.sd
    
    
    # start sampling loop
    for(i in 1:nrun) {

        #update posterior of lambda by conjugacy
        alpha2 <- lambda.alpha + sum(sales) + sum(d01)
        beta2 <- lambda.beta + K
        
        #simulate lambda2
        lambda2 <- rgamma(1, shape=alpha2, rate=beta2)
        
        
        
        #simulate beta2 by Metropolis-Hastings
        beta2 <- MCMCmetrop1R(logpost.beta, theta.init=beta1,
                         data=rbind(sales,d01),
                         thin=1, mcmc=1, burnin=10, tune=0.015,
                         verbose=0,  V=matrix(c(1,0,0,1),2,2))[1,]
                         #optim.lower=1e-6, optim.method="L-BFGS-B")[1,]


        #update and simulate eps1.mu and eps1.sd
        eps2 <- log(traffic - colSums(sales) - d01)
        eps2.mu2 <- rnorm(1, mean=mean(eps2), sd=eps2.sd1/sqrt(K))
        eps2.pr2 <- rgamma(1, shape=eps2.pr.alpha0+K/2, rate=eps2.pr.beta0+sum((eps2-eps2.mu2)^2)/2)
        eps2.sd2 <- 1/sqrt(eps2.pr2)
        
        
        
        #simulate d0 by discrete Metropolis-Hastings
        d02 <- d01
        d0.accept <- rep(0, K)
        for (j in 1:K) {
            dataj <- list(sales=sales[,j], traffic=traffic[j])
            sim <- discreteMH.norm(logpost.d0, start=d01[j], scale=8, nrun=10, 
                                   data=dataj, lambda=lambda2, beta=beta2, k=j, eps2.mu=eps2.mu2, eps2.sd=eps2.sd2)
            d02[j] <- sim$MC[10]
            d0.accept[j] <- sim$accept
        }
        cat("discreteMH acceptance rate was ", mean(d0.accept), "\n\n")

        

        # save the samples obtained in the current iteration
        lambdas[i,] = lambda2
        d0s[i,] = d02
        betas[i,] = beta2
        eps2.mus[i,] = eps2.mu2
        eps2.sds[i,] = eps2.sd2

        cat("Run:", i, "\tlambda:", lambda2, "\tbeta2:", beta2, "\teps2.mu2:", eps2.mu2, "\teps2.sd2:", eps2.sd2, "\n", sep=" ")
        
        lambda1 <- lambda2
        d01 <- d02
        beta1 <- beta2
        eps2.mu1 <- eps2.mu2
        eps2.sd1 <- eps2.sd2
    }
    
    # results
    return(list(lambdas=lambdas, d0s=d0s, betas=betas, eps2.mus=eps2.mus, eps2.sds=eps2.sds))

} # end function 'sample'




### initialize input before sampling
# beta prior ~ logN(beta.mu, beta.sg)
beta.mu <- c(-2.5, -2.5)
beta.sg <- matrix(c(0.5, 0, 0, 0.5), 2, 2)

# lambda ~ Gamma(lambda.alpha, lambda.beta)
lambda.alpha <- 0.5
lambda.beta <- 0.01

# prior of the precision of epsilon
eps2.pr.alpha0 <- 0.1
eps2.pr.beta0 <- 0.1

## initial sampling input
param0 <- list(beta=rep(0.1, L), d0=rep(10,K), eps2.sd=eps2.pr.alpha0/eps2.pr.beta0)
nrun <- 5000



### sample
z <- sample(data=observation3, parameters=param0, nrun=nrun)

save.image(file="MNL_Scenario3_A.RData")




### Visualize results
burnin <- 0.2*nrun

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

### save plots
require(ggplot2)

pdf('MNL_Scenario3_A.beta1.pdf', width = 8, height = 8)
ggplot(data=samples.beta.truncated) + geom_density(aes(x=X1), color="black")
dev.off()
pdf('MNL_Scenario3_A.beta2.pdf', width = 8, height = 8)
ggplot(data=samples.beta.truncated) + geom_density(aes(x=X2), color="black")
dev.off()


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



#plot eps2.mu and eps2.sd
samples.eps2 <- data.frame(mu=z$eps2.mus, sd=z$eps2.sds)
plot(samples.eps2$mu, type="l")
plot(samples.eps2$sd, type="l")

samples.eps2.truncated <- samples.eps2[burnin:nrun,]
quantile(samples.eps2.truncated$mu, c(.025,.5,.975))
quantile(samples.eps2.truncated$sd, c(.025,.5,.975))
colMeans(samples.eps2.truncated)
hist(samples.eps2.truncated$mu)
hist(samples.eps2.truncated$sd)

