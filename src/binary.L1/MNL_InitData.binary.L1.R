####
# Initialize choice data
# --- binary choice case
# --- with one-dimensional covaraite
# --- with inventory and stock-out info
####

#clear environment
rm(list=ls(all.names=TRUE))


Sys.setenv(LANG = "en")
setwd("~/Dropbox/RCode/Choice_Gibbs.git/src/binary.L1")

require("mvtnorm")



### Known parameters
M <- 2 # number of alternatives (the last alternative is dummy for no-purchase)
L <- 1 # number of covariates
K <- 90 # number of periods

#X is the attributes of the alternatives; in each period, [Xij] is an (M-1)*L matrix, i=1...M-1, j=1...L.
X_Mean <- rep(0, (M-1)*L)

#for each of the K periods, generate an X matrix
#dim of X_Mat is K, (M-1)*L
X_Mat <- rmvnorm(K, mean=X_Mean)


## true values of the parameters to be estimated
# beta is the MNL coefficient
beta <- c(0.3); # L-dimensional
# lambda is the poisson demand rate in each perid
lambda <- 50


### simulate data
# simulate number of demand per period, ~Poisson(lambda)
N <- rpois(K, lambda) 

#in each period, score of choice 1 and 2 (col) is exp(X*beta)
#dim of score is M by K
score <- rbind(t(exp(X_Mat %*% beta)), 1)
#choice probabilities in each period (M by K)
choice.prob <- apply(score, 2, function(x) x/sum(x)) # M*N matrix
#simulate actual choice in each period (M by K)
choice.mat <- matrix(0, nrow=M, ncol=K)
for (k in 1:K) {
    choice.mat[, k] <- rmultinom(1, N[k], choice.prob[, k])    
}    

rowMeans(choice.mat)


### organize data for model use
Inventory <- round(runif(K, min=20, max=30))
Demand <- choice.mat[1:M-1,]
Sales <- pmin(Demand, Inventory)
Stockout <- (Demand >= Inventory)
LostSales <- pmax(Demand-Inventory, 0)
NoPurchase <- choice.mat[M,]


### simulate traffic data with different accuracy
epsilon1.mean <- 3
epsilon1.xl <- rnorm(K, mean=epsilon1.mean, sd=1)
epsilon1.l <- rnorm(K, mean=epsilon1.mean, sd=0.5)
epsilon1.m <- rnorm(K, mean=epsilon1.mean, sd=0.1)
epsilon1.h <- rnorm(K, mean=epsilon1.mean, sd=0.05)
epsilon1.xh <- rnorm(K, mean=epsilon1.mean, sd=0.01)

Traffic.xl <- exp(epsilon1.xl) * colSums(choice.mat) 
Traffic.l <- exp(epsilon1.l) * colSums(choice.mat) 
Traffic.m <- exp(epsilon1.m) * colSums(choice.mat) 
Traffic.h <- exp(epsilon1.h) * colSums(choice.mat) 
Traffic.xh <- exp(epsilon1.xh) * colSums(choice.mat) 



rm(k)

save.image(file="MNL_InitData.binary.L1.RData")

