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
X_new <- 0.5


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



##generate demand samples
num <- 10
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



ggplot() +
    #geom_density(aes(x=data$lam0), color="black") +
    #geom_density(aes(x=data$lam1c), color="grey") #+
    geom_density(aes(x=rpois(50000, lambda=data$lam0)), color="black") +
    geom_density(aes(x=rpois(50000, lambda=data$lam1c)), color="grey") + 
    geom_density(aes(x=rpois(50000, lambda=50)), color="red")
