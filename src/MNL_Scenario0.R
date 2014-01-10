####
# Scenario 0. All choices are observed 
#
####

Sys.setenv(LANG = "en")
setwd("~/Dropbox/RCode/Choice_Gibbs.git/src")

require("MCMCpack")
require("mvtnorm")

#for %dopar%
library(doParallel)
cl <- makeCluster(detectCores(logical=FALSE))
registerDoParallel(cl)


### True parameters
# Number of choices
M <- 2

#XMAT is the observed attributes of the alternatives.We consider two products,each product has two attributes. Therefore, XMAT is a 2X2 matrix.
#by row: [X11 X12, X21 X22]
XMAT <- matrix(c(0.5, 0.3, 0.4, 0.1), nrow=M, byrow=TRUE);

#true coefficient of beta~mvN(b,W)
b <- c(-1, 3);
W <- matrix(c(4, 0.5, 0.5, 1), nrow=M); 



### simulate data
N <- 1000
beta <- rmvnorm(N, mean=b, sigma=W)
#score of choice 1 and 2 (col) by users (row) is exp(beta*X')
score <- exp(beta %*% t(XMAT))
#score of no-purchase is exp(0)=1
score <- cbind(score, 1)
#choice probabilities
choice.prob <- t(apply(score, 1, function(x) x/sum(x)))
#simulate actual choice
choice <-  t(apply(choice.prob, 1, function(x) rmultinom(1, 1, x)))
choice <- cbind(choice, choice %*% (1:(M+1)))

colMeans(choice)




##M-H log-posterior function
#log-posterior of beta
logpost.beta <- function(beta, data, bb, WW) {
    
    logScore <- cbind(beta %*% t(XMAT), 0)
    loglikelihood <- logScore[data] - log(sum(exp(logScore)))

    logprior <- dmvnorm(beta, mean=bb, sigma=WW, log=TRUE)
    
    return(loglikelihood + logprior)
}


##sampling
#one-step sampling function
sample1 = function(data, parameters) {
    
    mm <- length(parameters$b)
    nn <- length(data)
    
    b1 <- parameters$b
    W1 <- parameters$W
    beta1 <- parameters$beta
    
    
    #update posterior of W
    nu2 <- W.nu + nn
    Phi2 <- W.Phi + (t(beta1)-b1) %*% t((t(beta1)-b1))
    
    #simulate W2
    W2 <- riwish(nu2, Phi2)
    
    
    #update posterior of b
    m.beta <- colMeans(beta1)
    v.beta.i <- ginv(cov(beta1))
    v.i <- ginv(b.v)
    v2 <- ginv(v.i + nn * v.beta.i)
    m2 <- v2 %*% (v.i %*% b.m + nn * v.beta.i %*% m.beta)
    
    #simulate b2
    b2 <- as.vector(rmvnorm(1, mean=m2, sigma=v2))
    
    

    
    #simulate beta2 by Metropolic-Hasting
    #parallel
    beta2 <- foreach(betai=iter(beta1, by="row"), datai=iter(data), .combine="rbind", .packages=c("MCMCpack", "mvtnorm"), .export=c("logpost.beta", "XMAT"))  %dopar%  {
        MCMCmetrop1R(logpost.beta, theta.init=betai,
                     data=datai, bb=b2, WW=W2,
                     thin=1, mcmc=1, burnin=0, tune=2,
                     verbose=0, logfun=TRUE)[1,]
    }
    
    
    #return samples
    return(list(b=b2, W=W2, beta=beta2))
}


samplen = function(data, parameters, nrun=1000) {
    mm <- length(parameters$b)
    nn <- length(data)
    
    bs = array(0, dim=c(nrun, mm))
    Ws = array(0, dim=c(nrun, mm, mm))
    betas = array(0, dim=c(nrun, nn, mm))
    
    
    for(i in 1:nrun) {
        z = sample1(data, parameters)
        bs[i,] = z$b
        Ws[i,,] = z$W
        betas[i,,] = z$beta
        
        parameters = list(b=z$b, W=z$W, beta=z$beta)
    }
    
    # results
    output <- list(bs=bs, Ws=Ws, betas=betas)
    return(output)
} # end function 'samplen'




### initialize input before sampling
## beta ~ mvN(b, W)
# b ~ mvN(b.m, b.v)
b.m <- c(0, 0)
b.v <- matrix(c(100, 0, 0, 100), nrow=M);

# W ~ iWishart(nu, Phi)
W.nu <- 2
W.Phi <- matrix(c(100, 0, 0, 100), nrow=M)



#initial input (i.e., mu)
beta0 <- rmvnorm(N, mean=b.m, sigma=W.Phi)
param <- list(b=b.m, W=W.Phi, beta=beta0)
nrun <- 1000



#sample
z <- samplen(data=choice[,4], parameters=param, nrun=nrun)

stopCluster(cl)


## Visualize results
burnin <- 0.1*nrun

#plot b
#samples.b <- data.frame(matrix(unlist(z$bs), ncol=ncol(beta), byrow=T))
samples.b <- data.frame(z$bs)
plot(samples.b$X1, type="l")
plot(samples.b$X2, type="l")

plot(c(1,nrun), c(-2, 4), type="n", xlab="Samples", ylab="a", xaxt="n", yaxt="n")
lines(samples.b$X1, col="GREEN")
lines(samples.b$X2, col="RED")
axis(side=1)
axis(side=2)

samples.b.truncated <- samples.b[burnin:nrun,]
quantile(samples.b.truncated$X1, c(.025,.5,.975))
quantile(samples.b.truncated$X2, c(.025,.5,.975))
colMeans(samples.b.truncated)
hist(samples.b.truncated$X1)
hist(samples.b.truncated$X2)


#plot W
#samples.W <- data.frame(matrix(unlist(z$Ws), ncol=4, byrow=T))
samples.W <- data.frame(matrix(z$Ws, ncol=4))
plot(samples.W$X1, type="l")
plot(samples.W$X4, type="l")
plot(samples.W$X2, type="l")

plot(c(1,nrun), c(0, 5), type="n", xlab="Samples", ylab="a", xaxt="n", yaxt="n")
lines(samples.W$X1, col="GREEN")
lines(samples.W$X4, col="RED")
lines(samples.W$X2, col="BLUE")
axis(side=1)
axis(side=2)

samples.W.truncated <- samples.W[burnin:nrun,]
quantile(samples.W.truncated$X1, c(.025,.5,.975))
quantile(samples.W.truncated$X2, c(.025,.5,.975))
quantile(samples.W.truncated$X4, c(.025,.5,.975))
colMeans(samples.W.truncated)
hist(samples.W.truncated$X1)
hist(samples.W.truncated$X4)
hist(samples.W.truncated$X2)




#plot beta
samples.beta4 <- data.frame(z$betas[,4,])

plot(c(1,nrun), c(-5, 5), type="n", xlab="Samples", ylab="a", xaxt="n", yaxt="n")
lines(samples.beta4$X1, col="GREEN")
lines(samples.beta4$X2, col="BLUE")
axis(side=1)
axis(side=2)

samples.beta4.truncated <- samples.beta4[burnin:nrun,]
quantile(samples.beta4.truncated$X1, c(.025,.5,.975))
quantile(samples.beta4.truncated$X2, c(.025,.5,.975))
colMeans(samples.beta4.truncated)
hist(samples.beta4.truncated$X1)
hist(samples.beta4.truncated$X2)

