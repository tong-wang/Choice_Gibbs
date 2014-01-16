####
# Mixed MNL
# -- This is a UNIDENTIFIABLE model!
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
M <- 3 # number of alternatives (alternative 3 is dummy for no-purchase)
L <- 2 # number of covariates

#XMAT is the attributes of the alternatives; [Xij] is an M*L matrix, i=1...M, j=1...L.
#by col: [X11 X12; X21 X22; 0, 0]
XMAT <- matrix(c(0.4, 0.3, 
                 0.6, 0.1,
                 0, 0), 
               nrow=M, ncol=L, byrow=TRUE);

#true coefficient of beta~mvN(b,W)
b <- c(0.6, 0.3); # L-dimensional
W <- matrix(c(0.04, -0.01, -0.01, 0.01), nrow=L, ncol=L); # L*L matrix


### simulate data
N <- 10000 # number of data points
beta <- t(rmvnorm(N, mean=b, sigma=W)) # L*N matrix
#score of choice 1 and 2 (col) by users (row) is exp(X*beta)
score <- exp(XMAT %*% beta) # M*N matrix
#choice probabilities
choice.prob <- apply(score, 2, function(x) x/sum(x)) # M*N matrix
#simulate actual choice
choice.mat <- apply(choice.prob, 2, function(x) rmultinom(1, 1, x)) # M*N matrix
choice.vec <- t(1:M) %*%  choice.mat   # N-dimensional vector

rowMeans(choice.mat)




#log-posterior of beta, to be called by M-H algorithm
logpost.beta <- function(beta, data, bb, WW) {
    
    logScore <- XMAT %*% beta
    logLikelihood <- logScore[data] - log(sum(exp(logScore)))

    logPrior <- dmvnorm(beta, mean=bb, sigma=WW, log=TRUE)
    
    return(logLikelihood + logPrior)
}



sample = function(data, parameters, nrun=1000) {

    # get dimension parameters
    ll <- length(parameters$b)
    nn <- length(data)
    
    # initialize arrays to save samples
    bs = array(0, dim=c(nrun, ll))
    Ws = array(0, dim=c(nrun, ll,ll))
    betas = array(0, dim=c(nrun, ll, nn))
    
    # initial parameters
    b1 <- parameters$b
    W1 <- parameters$W
    beta1 <- parameters$beta
    
    
    # start sampling loop
    for(i in 1:nrun) {

        #update posterior of W by conjugacy
        nu2 <- W.nu + nn
        Psi2 <- W.Psi + (beta1-b1) %*% t(beta1-b1)
        
        #simulate W2
        W2 <- riwish(nu2, Psi2)
        
        
        ## update posterior of b
        # sample mean of beta
        m.beta <- rowMeans(beta1)
        
        # assume b has diffuse prior, so b'' ~ N(m.beta, W2/N)
        #simulate b2
        b2 <- as.vector(rmvnorm(1, mean=m.beta, sigma=W2/nn))
        
        
        #simulate beta2 by Metropolic-Hasting
        #parallel
        beta2 <- foreach(betai=iter(beta1, by="column"), datai=iter(data, by="column"), .combine="cbind", .packages=c("MCMCpack", "mvtnorm"), .export=c("logpost.beta", "XMAT"))  %dopar%  {
            #print(betai)
            #print(datai)
            #print(logpost.beta(as.vector(betai), as.vector(datai), b2, W2))
            
            MCMCmetrop1R(logpost.beta, theta.init=betai,
                         data=datai, bb=b2, WW=W2,
                         thin=1, mcmc=1, burnin=0, tune=2,
                         verbose=0, logfun=TRUE)[1,]
        }
        
        

        
        # save the samples obtained in the current iteration
        bs[i,] = b2
        Ws[i,,] = W2
        betas[i,,] = beta2
        
        cat("Run:", i, "\tb:", b2, "\tW:", W2, "\n", sep=" ")
        
        b1 <- b2
        W1 <- W2
        beta1 <- beta2
    }
    
    # results
    return(list(bs=bs, Ws=Ws, betas=betas))

} # end function 'sample'




### initialize input before sampling
## beta ~ mvN(b, W)
# b has a diffuse prior
# W ~ iWishart(nu, Psi)
W.nu <- L
W.Psi <- diag(L) * L


## initial sampling input
param0 <- list(b=c(0, 0), W=W.Psi, beta=t(rmvnorm(N, mean=c(0, 0), sigma=W.Psi/W.nu)))
nrun <- 1000



### sample
z <- sample(data=choice.vec, parameters=param0, nrun=nrun)

stopCluster(cl)


### test results
b2 <- c(-1.835909, 1.347521)
W2 <-  matrix(c(0.1995626, 0.2219391, 0.2219391, 0.2571687), nrow=L, ncol=L)
beta2 <- t(rmvnorm(N, mean=b2, sigma=W2)) # L*N matrix
score2 <- exp(XMAT %*% beta2) # M*N matrix
choice.prob2 <- apply(score2, 2, function(x) x/sum(x)) # M*N matrix
choice.mat2 <- apply(choice.prob2, 2, function(x) rmultinom(1, 1, x)) # M*N matrix
rowMeans(choice.mat2)
rowMeans(choice.mat)




### Visualize results
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
samples.beta4 <- data.frame(z$betas[,,4])
plot(samples.beta4$X1, type="l")
plot(samples.beta4$X2, type="l")

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

