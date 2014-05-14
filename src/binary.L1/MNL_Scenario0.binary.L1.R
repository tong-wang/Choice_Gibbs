####
# Scenario 0. All choices are observed 
#
####

Sys.setenv(LANG = "en")
setwd("~/Dropbox/RCode/Choice_Gibbs.git/src/binary.L1")

require("mvtnorm")
source(file="Metropolis-Hastings.R")


## Load simulated choice data (NEED TO RUN MNL_InitData.binary.R TO GENERATE THE DATA FIRST)
load(file="MNL_InitData.binary.L1.RData")



observation0 <- choice.mat



### Estimating beta by M-H sampling
# log-posterior of beta, to be called by M-H algorithm
logpost.beta <- function(beta, data) {
    
    score <- rbind(t(exp(X_Mat %*% beta)), 1)
    choice.prob <- apply(score, 2, function(x) x/sum(x)) # M*N matrix
    
    logLikelihood <- data*log(choice.prob)
    
    logprior <- dmvnorm(beta, mean=beta.mu, sigma=beta.sg, log=TRUE)
    
    return(sum(logLikelihood) + logprior)
    
}



## initialize input before sampling
# beta prior ~ N(beta.mu, beta.sg)
beta.mu <- rep(0, L)
beta.sg <- 100*diag(L)


nrun <- 10000
burnin <- 0.5
start <- burnin*nrun+1


#direct sampling beta
MH <- MH.mvnorm(logpost.beta, sigma=diag(L), scale=0.08, start=rep(0.1, L), nrun = nrun, data=observation0)
cat("MH acceptance rate: ", MH$accept, "\n")

betas <- MH$MC
plot(betas, type="l")

beta.estm <- mean(betas[start:nrun,])


### Estimate lambda by conjugate prior

# lambda ~ Gamma(lambda.alpha, lambda.beta)
lambda.alpha <- 0.005
lambda.beta <- 0.0001

#update posterior of lambda by conjugacy
alpha2 <- lambda.alpha + sum(observation0)
beta2 <- lambda.beta + K

lambda.estm <- alpha2/beta2

score.estm <- exp(matrix(X_Mean, M-1,L) %*% beta.estm)
score.estm <- c(score.estm, 1)
choice.prob.estm <- score.estm / sum(score.estm)
demand.estm <- lambda.estm * choice.prob.estm




z0 <- list(lambdas=rgamma(nrun, shape=alpha2, rate=beta2), betas=betas)

save(z0, observation0, file="MNL_Scenario0.binary.L1.RData")



#plot lambda
samples.lambda <- z0$lambdas
plot(samples.lambda, type="l")


samples.lambda.truncated <- samples.lambda[start:nrun]
quantile(samples.lambda.truncated, c(.025,.5,.975))
mean(samples.lambda.truncated)
hist(samples.lambda.truncated)



#plot beta
samples.beta <- z0$betas
plot(samples.beta, type="l")


samples.beta.truncated <- samples.beta[start:nrun]
quantile(samples.beta.truncated, c(.025,.5,.975))
mean(samples.beta.truncated)
hist(samples.beta.truncated)


### save plots
require(ggplot2)

pdf('MNL_Scenario0.binary.L1.lambda.pdf', width = 8, height = 8)
ggplot(data=data.frame(samples.lambda.truncated)) + geom_density(aes(x=samples.lambda.truncated), color="black")
dev.off()

pdf('MNL_Scenario0.binary.L1.beta.pdf', width = 8, height = 8)
ggplot(data=data.frame(samples.beta.truncated)) + geom_density(aes(x=samples.beta.truncated), color="black")
dev.off()


