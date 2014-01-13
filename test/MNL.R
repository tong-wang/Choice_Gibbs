####
# Standard Multinomial Logit model
#   Estimation using Metropolic-Hasting algorithm
####

Sys.setenv(LANG = "en")
setwd("~/Dropbox/RCode/Choice_Gibbs.git/test")

require("MCMCpack")
require("mvtnorm")


### True parameters (M>=L for identifiability)
M <- 6 # number of alternatives (the last alternative is dummy for no-purchase)
L <- 3 # number of covariates

#XMAT is the attributes of the alternatives; [Xij] is an M*L matrix, i=1...M, j=1...L.
#by col: [X11 X12; X21 X22; 0, 0]
XMAT <- matrix(c(4, 3, 5, 
                 6, 1, 2, 
                 2, 2, 3,
                 5, 4, 2,
                 4, 6, 7,
                 1, 3, 2),
               nrow=M, ncol=L, byrow=TRUE);

#true coefficient of beta
beta <- c(0.06, 0.04, 0.03); # L-dimensional


### simulate data
N <- 50000 # number of data points
#score of choice 1 and 2 (col) by users (row) is exp(X*beta)
score <- exp(XMAT %*% beta)
#choice probabilities
choice.prob <- score / sum(score)
#simulate actual choice
choice.mat <- rmultinom(N, 1, choice.prob)
choice.vec <- t(1:M) %*%  choice.mat   # N-dimensional vector

rowMeans(choice.mat)
data <- rowSums(choice.mat)




### M-H sampling
#log-posterior of beta, to be called by M-H algorithm
logpost.beta <- function(beta, data) {
    
    score <- exp(XMAT %*% beta)
    choice.prob <- score / sum(score)

    logLikelihood <- 0
    for (i in 1:length(choice.prob)) {
        logLikelihood <- logLikelihood + data[i]*log(choice.prob[i])
    }
    
    return(logLikelihood)
}



#sampling
z <- MCMCmetrop1R(logpost.beta, theta.init=rep(0,L),
             data=data,
             thin=10, mcmc=100000, burnin=1000, tune=1.2,
             verbose=10000, logfun=TRUE)


summary(z)
plot(z)

beta.estm <- colMeans(z)
score.estm <- exp(XMAT %*% beta.estm)
choice.prob.estm <- score.estm / sum(score.estm)
choice.prob.estm
