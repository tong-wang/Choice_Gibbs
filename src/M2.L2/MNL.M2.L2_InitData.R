####
# Initialize choice data
# --- with inventory and stock-out info
####

#clear environment
rm(list=ls(all.names=TRUE))


Sys.setenv(LANG = "en")
setwd("~/Dropbox/RCode/Choice_Gibbs.git/src/M2.L2")

require("mvtnorm")



### Known parameters
M <- 2 # number of alternatives (not including the no-purchase option)
L <- 2 # number of covariates (not including dummy for constant)
K <- 90 # number of periods


## Generate X_Mat, the covariates matrix
# In each period k, covariates X[i,j] is an M*L matrix, i=1...M, j=1...L. X[i,j] is coerced into a row vector
# for each of the K periods, generate the X_Mat matrix with K rows of X, dimension = K * (M*L)
X1 <- rmvnorm(K, mean=c(2, 2.4), sigma=matrix(c(0.09, 0.05, 0.05, 0.16), 2, 2)) # price covariate ~ Normal
X2 <- cbind(sample(x=1:3, size=K, replace=TRUE), sample(x=1:3, size=K, replace=TRUE)) # quality covariate ~ DiscreteUniform
X_Mat <- cbind(X1, X2)


### true values of the parameters to be estimated
# lambda is the poisson demand rate in each perid
lambda <- 50

# beta is the MNL coefficient: utility of alternative m = beta.const[m] + beta.coef * X[m]
beta.coef <- c(-2, 0.5); # coefficient terms, L-dimensional
beta.const <- c(3, 4); # constant terms, M-dimensional
# in learning, we estimate a transformation of beta: betaT = beta[1] , beta[2]/beta[1], ..., beta[L]/beta[1]
beta <- c(beta.coef, beta.const)
betaT <- c(beta[1], beta[2:(L+M)] / beta[1])


### simulate data
# simulate number of demand per period, ~Poisson(lambda)
N <- rpois(K, lambda) 

#in each period, score of choice 1 and 2 is exp(utility), score of no-purchase (choice M+1) is exp(0)
#dim of score is M+1 by K
score <- rbind(apply(X_Mat, 1, function (x) exp(beta.const + matrix(x, nrow=M, ncol=L) %*% beta.coef)), 1)
#choice probabilities in each period (M+1 by K)
choice.prob <- apply(score, 2, function(x) x/sum(x))
#simulate actual choice in each period (M+1 by K)
choice.mat <- matrix(0, nrow=M+1, ncol=K)
for (k in 1:K) {
    choice.mat[, k] <- rmultinom(1, N[k], choice.prob[, k])    
}    

rowMeans(choice.mat)


Inventory <- matrix(round(runif((M-1)*K, min=10, max=20)), M, K)
Demand <- choice.mat[1:M,]
Sales <- pmin(Demand, Inventory)
Stockout <- (Demand >= Inventory)
LostSales <- pmax(Demand-Inventory, 0)
NoPurchase <- choice.mat[M+1,]


### simulate traffic data with different accuracy
epsilon1.mean <- 3
epsilon1.xl <- rnorm(K, mean=epsilon1.mean, sd=1)
epsilon1.l <- rnorm(K, mean=epsilon1.mean, sd=0.5)
epsilon1.m <- rnorm(K, mean=epsilon1.mean, sd=0.1)
epsilon1.h <- rnorm(K, mean=epsilon1.mean, sd=0.05)
epsilon1.xh <- rnorm(K, mean=epsilon1.mean, sd=0.01)

TrafficM.xl <- exp(epsilon1.xl) * colSums(choice.mat) 
TrafficM.l <- exp(epsilon1.l) * colSums(choice.mat) 
TrafficM.m <- exp(epsilon1.m) * colSums(choice.mat) 
TrafficM.h <- exp(epsilon1.h) * colSums(choice.mat) 
TrafficM.xh <- exp(epsilon1.xh) * colSums(choice.mat) 



rm(k)

save.image(file="MNL.M2.L2_InitData.RData")

