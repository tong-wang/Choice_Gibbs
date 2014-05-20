####
# Initialize choice data
# --- binary choice case
# --- with inventory and stock-out info
####

#clear environment
rm(list=ls(all.names=TRUE))


Sys.setenv(LANG = "en")
setwd("~/Dropbox/RCode/Choice_Gibbs.git/src/binary")

require("mvtnorm")



### Known parameters
M <- 2 # number of alternatives (the last alternative is the dummy for no-purchase)
L <- 3 # number of covariates (the last covariate takes constant value 1, which is the dummy for the constant beta)
K <- 60 # number of periods

#X is the attributes of the alternatives; in each period, [Xij] is an (M-1)*L matrix, i=1...M-1, j=1...L.
#for each of the K periods, generate an X matrix
#dim of X_Mat is K, (M-1)*L
X1 <- rnorm(K, mean=2, sd=0.3) # price covariate ~ N(2, 0.3)
X2 <- sample(x=1:5, size=K, replace=TRUE) # quality covariate ~ DiscreteUniform(1,5)
X3 <- rep(1, K) # dummy covariate = 1
X_Mat <- cbind(X1, X2, X3)


## true values of the parameters to be estimated
# beta is the MNL coefficient
beta <- c(-3, 1, 3); # L-dimensional
betaT <- c(beta[1], beta[2:L]/beta[1])


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

save.image(file="MNL_InitData.binary.RData")

