require("plyr")
require("ggplot2")

setwd("~/Dropbox/RCode/Choice_Gibbs.git/src/M1.L1/output/K90L50")

load("MNL.M1.L1_Scenario0.RData")
load("MNL.M1.L1_Scenario1.RData")
load("MNL.M1.L1_Scenario3M.xl.RData")
load("MNL.M1.L1_Scenario3M.l.RData")
load("MNL.M1.L1_Scenario3M.m.RData")
load("MNL.M1.L1_Scenario3M.h.RData")
load("MNL.M1.L1_Scenario3M.xh.RData")


load("MNL.M1.L1_Scenario0c.RData")
load("MNL.M1.L1_Scenario1c.RData")
load("MNL.M1.L1_Scenario3Mc.xl.RData")
load("MNL.M1.L1_Scenario3Mc.l.RData")
load("MNL.M1.L1_Scenario3Mc.m.RData")
load("MNL.M1.L1_Scenario3Mc.h.RData")
load("MNL.M1.L1_Scenario3Mc.xh.RData")


nrun <- 50000
burnin <- 0.5
start <- burnin*nrun+1

# generate a sample index (remove burn-in and thinning by a factor of 10)
thin <- 10
sample.index <- ((1:nrun) >= start) & ((1:nrun) %% thin == 0)



##merge estimation into a big dataframe
data <- data.frame(
    lam0=z0$lambdas[sample.index],
    lam1=z1$lambdas[sample.index,], 
    lam3.xl=z3M.xl$lambdas[sample.index,], 
    lam3.l=z3M.l$lambdas[sample.index,], 
    lam3.m=z3M.m$lambdas[sample.index,], 
    lam3.h=z3M.h$lambdas[sample.index,],
    lam3.xh=z3M.xh$lambdas[sample.index,],
    lam0c=z0c$lambdas[sample.index],
    lam1c=z1c$lambdas[sample.index,], 
    lam3c.xl=z3Mc.xl$lambdas[sample.index,], 
    lam3c.l=z3Mc.l$lambdas[sample.index,], 
    lam3c.m=z3Mc.m$lambdas[sample.index,], 
    lam3c.h=z3Mc.h$lambdas[sample.index,],
    lam3c.xh=z3Mc.xh$lambdas[sample.index,],
    beta0=z0$betas[sample.index,],
    beta1=z1$betas[sample.index,], 
    beta3.xl=z3M.xl$betas[sample.index,], 
    beta3.l=z3M.l$betas[sample.index,], 
    beta3.m=z3M.m$betas[sample.index,], 
    beta3.h=z3M.h$betas[sample.index,],
    beta3.xh=z3M.xh$betas[sample.index,],
    beta0c=z0c$betas[sample.index,],
    beta1c=z1c$betas[sample.index,], 
    beta3c.xl=z3Mc.xl$betas[sample.index,], 
    beta3c.l=z3Mc.l$betas[sample.index,], 
    beta3c.m=z3Mc.m$betas[sample.index,], 
    beta3c.h=z3Mc.h$betas[sample.index,],
    beta3c.xh=z3Mc.xh$betas[sample.index,]
)
data$id <- 1:nrow(data)
colMeans(data)



# generate a new set of covaraites for prediction
#X_new <- c(-2, -1, 0, 1, 2)
X_new <- 1
Cost_new <- 0.7

##generate lambda samples
sample.lambda <- function (row) {
    
    row$lambda0 <- row$lam0 * (exp(X_new * row$beta0.1 + row$beta0.2) / (1 + exp(X_new * row$beta0.1 + row$beta0.2)))
    row$lambda0c <- row$lam0c * (exp(X_new * row$beta0c.1 + row$beta0c.2) / (1 + exp(X_new * row$beta0c.1 + row$beta0c.2)))
    row$lambda1 <- row$lam1 * (exp(X_new * row$beta1.1 + row$beta1.2) / (1 + exp(X_new * row$beta1.1 + row$beta1.2)))
    row$lambda1c <- row$lam1c * (exp(X_new * row$beta1c.1 + row$beta1c.2) / (1 + exp(X_new * row$beta1c.1 + row$beta1c.2)))
    row$lambda3.m <- row$lam3.m * (exp(X_new * row$beta3.m.1 + row$beta3.m.2) / (1 + exp(X_new * row$beta3.m.1 + row$beta3.m.2)))
    row$lambda3c.m <- row$lam3c.m * (exp(X_new * row$beta3c.m.1 + row$beta3c.m.2) / (1 + exp(X_new * row$beta3c.m.1 + row$beta3c.m.2)))
    
    return(row)
}

data <- ddply(data, .(id), sample.lambda)
summary(data)

#plot lambda distribution
ggplot(data=data) + 
    geom_density(aes(x=lambda0), color="black") + 
    geom_density(aes(x=lambda1), color="grey") +
    geom_density(aes(x=lambda3.m), color="blue") +
    geom_density(aes(x=lambda0c), color="black", linetype="dashed") + 
    geom_density(aes(x=lambda1c), color="grey", linetype="dashed") +
    geom_density(aes(x=lambda3c.m), color="blue", linetype="dashed")




### overall mean and sd of profits under different inventory
##generate demand samples
num <- 100
data2 <- data.frame(
    demand0 = rpois(n=num*nrow(data), lambda=data$lambda0),
    demand0c = rpois(n=num*nrow(data), lambda=data$lambda0c),
    demand1 = rpois(n=num*nrow(data), lambda=data$lambda1),
    demand1c = rpois(n=num*nrow(data), lambda=data$lambda1c),
    demand3.m = rpois(n=num*nrow(data), lambda=data$lambda3.m),
    demand3c.m = rpois(n=num*nrow(data), lambda=data$lambda3c.m)
)
summary(data2)

#plot demand distribution
ggplot(data=data2) + 
    geom_density(aes(x=demand0), color="black",) + 
    geom_density(aes(x=demand1), color="grey") +
    geom_density(aes(x=demand3.m), color="blue") +
    geom_density(aes(x=demand0c), color="black", linetype="dashed") + 
    geom_density(aes(x=demand1c), color="grey", linetype="dashed") +
    geom_density(aes(x=demand3c.m), color="blue", linetype="dashed")


#optimal inventory
quantile(data2$demand0, c(.1,.3,.5,.7,.9))
quantile(data2$demand1, c(.1,.3,.5,.7,.9))
quantile(data2$demand1c, c(.1,.3,.5,.7,.9))
quantile(data2$demand3.m, c(.1,.3,.5,.7,.9))
quantile(data2$demand3c.m, c(.1,.3,.5,.7,.9))


#profit
profit0 <- NULL
profit1 <- NULL
profit1c <- NULL
profit3.m <- NULL

for (inv in 21:50) {
    profit0 <- cbind(profit0, X_new * pmin(data2$demand0, inv) - Cost_new * inv)
    profit1 <- cbind(profit1, X_new * pmin(data2$demand1, inv) - Cost_new * inv)
    profit1c <- cbind(profit1c, X_new * pmin(data2$demand1c, inv) - Cost_new * inv)
    profit3.m <- cbind(profit3.m, X_new * pmin(data2$demand3.m, inv) - Cost_new * inv)
}

#mean-variance
mean0 <- colMeans(profit0)
mean1 <- colMeans(profit1)
mean1c <- colMeans(profit1c)
mean3.m <- colMeans(profit3.m)
sd0 <- apply(profit0, 2, sd)
sd1 <- apply(profit1, 2, sd)
sd1c <- apply(profit1c, 2, sd)
sd3.m <- apply(profit3.m, 2, sd)

#plot efficient frontier
plot(mean0, sd0, type="l")
lines(mean1, sd1, type="l", col="grey")
lines(mean1c, sd1c, type="l", col="darkgrey")
lines(mean3.m, sd3.m, type="l", col="red")

#expectation-maximization inv
20 + which.max(mean0)
20 + which.max(mean1)
20 + which.max(mean1c)
20 + which.max(mean3.m)




### profit 
N.exp <- 1000

profit.given.lambda <- function(row) {
    d0 <- rpois(N.exp, row$lambda0)
    profit0 <- X_new * pmin(d0, 47) - Cost_new * 47
    row$mean0 <- mean(profit0)
    row$sd0 <- sd(profit0)
    
    d1 <- rpois(N.exp, row$lambda1)
    profit1 <- X_new * pmin(d1, 45) - Cost_new * 45
    row$mean1 <- mean(profit1)
    row$sd1 <- sd(profit1)
    
    
    return(row)
}


data3 <- data[,c("lambda0", "lambda0c", "lambda1", "lambda1c", "lambda3.m", "lambda3c.m")]

data3 <- ddply(data3, .(lambda0), profit.given.lambda)


plot(data3$mean0, data3$sd0)
plot(data3$mean1, data3$sd1)

ggplot(data=data3) + 
    geom_density(aes(x=mean0), color="black",) + 
    geom_density(aes(x=mean1), color="grey") #+
    geom_density(aes(x=demand3.m), color="blue") +
    geom_density(aes(x=demand0c), color="black", linetype="dashed") + 
    geom_density(aes(x=demand1c), color="grey", linetype="dashed") +
    geom_density(aes(x=demand3c.m), color="blue", linetype="dashed")
