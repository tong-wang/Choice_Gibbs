####
# Initialize choice data
# --- binary choice case
# --- with one-dimensional covaraite
# --- with inventory and stock-out info
####

#clear environment
rm(list=ls(all.names=TRUE))


Sys.setenv(LANG = "en")
setwd("~/Dropbox/RCode/Choice_Gibbs.git/src/binary.L2")



### Known parameters
M <- 2 # number of alternatives + no-purchase option
L <- 2 # number of covariates + dummy for constant
K <- 90 # number of periods

## Generate X_Mat, the covariates matrix
# In each period k, covariates X[i,j] is an (M-1)*L matrix, i=1...M-1, j=1...L. X[i,j] is coerced into a row vector
# for each of the K periods, generate the X_Mat matrix with K rows of X, dimension = K, (M-1)*L
X1 <- rnorm(K, mean=2, sd=0.3) # price covariate ~ N(2, 0.3)
XL <- rep(1, K) # dummy covariate fixed at 1
X_Mat <- cbind(X1, XL)


## true values of the parameters to be estimated
# lambda is the poisson demand rate in each perid
lambda <- 50
# beta is the MNL coefficient: score of alternative m = beta * X[m, ]
beta <- c(-3, 6); # L-dimensional
# in learning, we estimate a transformation of beta: betaT = beta[1] , beta[2]/beta[1], ..., beta[L]/beta[1]
betaT <- c(beta[1], beta[2:L]/beta[1])

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

TrafficM.xl <- exp(epsilon1.xl) * colSums(choice.mat) 
TrafficM.l <- exp(epsilon1.l) * colSums(choice.mat) 
TrafficM.m <- exp(epsilon1.m) * colSums(choice.mat) 
TrafficM.h <- exp(epsilon1.h) * colSums(choice.mat) 
TrafficM.xh <- exp(epsilon1.xh) * colSums(choice.mat) 



rm(k)

save.image(file="MNL.binary.L2_InitData.RData")

