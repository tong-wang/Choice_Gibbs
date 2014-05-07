####
# Initialize choice data
# --- binary choice case
####

Sys.setenv(LANG = "en")
setwd("~/Dropbox/RCode/Choice_Gibbs.git/src")

require("mvtnorm")



### Known parameters
M <- 2 # number of alternatives (alternative 3 is dummy for no-purchase)
L <- 2 # number of covariates
K <- 360 # number of periods

#X is the attributes of the alternatives; in each period, [Xij] is an (M-1)*L matrix, i=1...M-1, j=1...L.
#by row: [X11 X12; X21 X22]
X_Mean <- c(4, 3)

#for each of the K periods, generate an X matrix
#dim of X_Mat is K, (M-1)*L
X_Mat <- rmvnorm(K, mean=X_Mean, sigma=diag(0.5, length(X_Mean)))


## true values of the parameters to be estimated
# beta is the MNL coefficient
beta <- c(0.06, 0.03); # L-dimensional
# lambda is the poisson demand rate in each perid
lambda <- 100


### simulate data
# simulate number of demand per period, ~Poisson(lambda)
N <- rpois(K, lambda) 

#in each period, score of choice 1 and 2 (col) is exp(X*beta)
#dim of score is M by K
score <- apply(X_Mat, 1, function (x) exp(matrix(c(x,0,0), nrow=M, ncol=L, byrow=TRUE) %*% beta))
#choice probabilities in each period (M by K)
choice.prob <- apply(score, 2, function(x) x/sum(x)) # M*N matrix
#simulate actual choice in each period (M by K)
choice.mat <- matrix(0, nrow=M, ncol=K)
for (k in 1:K) {
    choice.mat[, k] <- rmultinom(1, N[k], choice.prob[, k])    
}    

rowMeans(choice.mat)

rm(k)

save.image(file="MNL_InitData.binary.RData")

