Sys.setenv(LANG = "en")
setwd("~/Dropbox/RCode/")

require("MCMCpack")
require("mvtnorm")


#true coefficient of beta~mvN(b,W)
b <- c(0.6, 0.3);
W <- matrix(c(0.04, -0.01, -0.01, 0.01), nrow=2, ncol=2); 


### simulate data
N <- 1000
beta <- t(rmvnorm(N, mean=b, sigma=W))
hist(beta[1,])
hist(beta[2,])



#one-step sampling function
sample1 = function(data, parameters) {
    
    nn <- ncol(data)
    
    b1 <- parameters$b
    W1 <- parameters$W
    
    
    #update posterior of W
    nu2 <- W.nu + nn
    Psi2 <- W.Psi + (data-b1) %*% t(data-b1)
    
    #simulate W2
    W2 <- riwish(nu2, Psi2)
    
    
    #update posterior of b using Normal-Normal conjugate
    m.data <- rowMeans(data)
    #v.data.i <- ginv(W2)
    #v.i <- ginv(b.v)
    #v2 <- ginv(v.i + nn * v.data.i)
    #m2 <- v2 %*% (v.i %*% b.m + nn * v.data.i %*% m.data)
    
    #simulate b2
    #b2 <- as.vector(rmvnorm(1, mean=m2, sigma=v2))

    # Alternatively, assume b has diffuse prior, so b'' ~ N(m.beta, W2/N)
    b2 <- as.vector(rmvnorm(1, mean=m.data, sigma=W2/nn))
    
    
    #return samples
    return(list(b=b2, W=W2))
}


samplen = function(data, parameters, nrun=1000) {
    bs = array(0, dim=c(nrun, 2))
    Ws = array(0, dim=c(nrun, 2, 2))
    
    
    for(i in 1:nrun) {
        z = sample1(data, parameters)
        bs[i,] = z$b
        Ws[i,,] = z$W
        
        parameters = z
    }
    
    # results
    output <- list(bs=bs, Ws=Ws)
    return(output)
} # end function 'samplen'




### initialize input before sampling
# beta ~ mvN(b, W)
# b ~ mvN(b.m, b.v)
b.m <- c(0, 0)
#b.v <- matrix(c(100, 0, 0, 100), nrow=2);

# W ~ iWishart(nu, Phi)
W.nu <- 2
W.Psi <-  diag(2) * 2

#initial input (i.e., mu)
param <- list(b=b.m, W=W.Psi)
nrun <- 1000

#sample
z <- samplen(data=beta, parameters=param, nrun=nrun)



## Visualize results
burnin <- 0.1*nrun

#plot b
samples.b <- data.frame(z$bs)
plot(samples.b$X1, type="l")
plot(samples.b$X2, type="l")

plot(c(1,nrun), c(0.2, 0.7), type="n", xlab="Samples", ylab="a", xaxt="n", yaxt="n")
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

samples.W <- data.frame(matrix(z$Ws, ncol=4))
plot(samples.W$X1, type="l")
plot(samples.W$X4, type="l")
plot(samples.W$X2, type="l")

plot(c(1,nrun), c(-0.02, 0.05), type="n", xlab="Samples", ylab="a", xaxt="n", yaxt="n")
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





### Theoretical result using Normal-iWishart
# prior of beta ~ Normal-iWishart(mu, k, nu, Phi)

beta.mu0 <- c(0, 0)
beta.k0 <- 2
beta.nu0 <- 2
beta.Phi0 <- matrix(c(2, 0, 0, 2), nrow=2)

data.mu <- rowMeans(beta)
data.SS <- cov(t(beta))*(N-1)
beta.mu2 <- (beta.k0 * beta.mu0 + N * data.mu)/(beta.k0 + N)
beta.k2 <- beta.k0 + N
beta.nu2 <- beta.nu0 + N
beta.Phi2 <- beta.Phi0 + data.SS + (beta.k0 * N) * (data.mu - beta.mu0) %*% t(data.mu-beta.mu0) / (beta.k0 + N)

beta.mu2
beta.Phi2/(beta.nu2 - 2 -1)

