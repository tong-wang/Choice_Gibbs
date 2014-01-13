################################################################
# Standard Multinomial Logit
#   Estimation using Metropolic-Hasting algorithm
#
# Setup:
#   -- M alternatives
#   -- Utility of alternative i is U_i = S_i + epsilon
#   -- Task: estimate S (M-dimensional vector)
#
# Note:
#   -- The last alternative has zero X (this is a dummy for the no-purchase option), which helps to normalize the utilities.
################################################################

Sys.setenv(LANG = "en")
setwd("~/Dropbox/RCode/Choice_Gibbs.git/test")

require("MCMCpack")


### True parameters
M <- 6 # number of alternatives (the last alternative is dummy for no-purchase)

#true coefficient of S
S <- c(1, 2, 3, 1.6, 2.3, 0); # L-dimensional


### simulate data
N <- 50000 # number of data points
#score of choice 1 and 2 (col) by users (row) is exp(X*beta)
score <- exp(S)
#choice probabilities
choice.prob <- score / sum(score)
#simulate actual choice
data <- rmultinom(1, N, choice.prob)




### M-H sampling
# log-posterior of beta, to be called by M-H algorithm
# assuming diffuse prior
logpost.beta <- function(s, data) {
    
    score <- exp(c(s,0))
    choice.prob <- score / sum(score)

    logLikelihood <- 0
    for (i in 1:length(choice.prob)) {
        logLikelihood <- logLikelihood + data[i]*log(choice.prob[i])
    }
    
    return(logLikelihood)
}



#sampling
z <- MCMCmetrop1R(logpost.beta, theta.init=rep(0, M-1),
             data=data,
             thin=10, mcmc=100000, burnin=10000, tune=1,
             verbose=10000, logfun=TRUE)


summary(z)
plot(z)

S.estm <- colMeans(z)
S.estm
