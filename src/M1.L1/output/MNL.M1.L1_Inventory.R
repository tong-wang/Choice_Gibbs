####
# Inventory optimization
####

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



# merge estimation into a big dataframe
posteriors <- data.frame(
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
posteriors$id <- 1:nrow(posteriors)
colMeans(posteriors)



# generate effective lambda (lambda * choice.prob) with given covariates
Cost <- 1.2
Price <- 1.5 * Cost

# generate effective lambda samples
calc.elambda <- function (row, scenario) {
    
    out <- data.frame(id=row$id)
    out$elambda <- row[, paste0("lam", scenario)]  / (1 + exp(- Price * row[, paste0("beta",scenario, ".1")] - row[, paste0("beta",scenario, ".2")]))
    
    return(out)
}

elambda0 <- ddply(posteriors, .(id), calc.elambda, scenario="0")
elambda1 <- ddply(posteriors, .(id), calc.elambda, scenario="1")
elambda3.m <- ddply(posteriors, .(id), calc.elambda, scenario="3.m")
elambda0c <- ddply(posteriors, .(id), calc.elambda, scenario="0c")
elambda1c <- ddply(posteriors, .(id), calc.elambda, scenario="1c")
elambda3c.m <- ddply(posteriors, .(id), calc.elambda, scenario="3c.m")
head(elambda1)

# plot effective lambda distribution
ggplot() + 
    geom_density(aes(x=elambda0$elambda), color="black") + 
    geom_density(aes(x=elambda1$elambda), color="grey") +
    geom_density(aes(x=elambda3.m$elambda), color="blue") +
    geom_density(aes(x=elambda0c$elambda), color="black", linetype="dashed") + 
    geom_density(aes(x=elambda1c$elambda), color="grey", linetype="dashed") +
    geom_density(aes(x=elambda3c.m$elambda), color="blue", linetype="dashed")




### mean and sd of profits under different inventory

## generate demand samples
num <- 100

generate.demand <- function(row) {
    
    out <- data.frame(
        id = rep(row$id, num),
        #elambda = rep(row$elambda, num),
        demand = rpois(num, row$elambda)
    )
    
    return(out)
}

demand0 <- ddply(elambda0, .(id), generate.demand)
demand1 <- ddply(elambda1, .(id), generate.demand)
demand3.m <- ddply(elambda3.m, .(id), generate.demand)
demand0c <- ddply(elambda0c, .(id), generate.demand)
demand1c <- ddply(elambda1c, .(id), generate.demand)
demand3c.m <- ddply(elambda3c.m, .(id), generate.demand)
head(demand0)

#plot demand distribution
ggplot() + 
    geom_density(aes(x=demand0$demand), color="black",) + 
    geom_density(aes(x=demand1$demand), color="grey") +
    geom_density(aes(x=demand3.m$demand), color="blue") +
    geom_density(aes(x=demand0c$demand), color="black", linetype="dashed") + 
    geom_density(aes(x=demand1c$demand), color="grey", linetype="dashed") +
    geom_density(aes(x=demand3c.m$demand), color="blue", linetype="dashed")


#optimal inventory, aka quantiles
quantile(demand0$demand, c(.1,.3,.5,.7,.9))
quantile(demand1$demand, c(.1,.3,.5,.7,.9))
quantile(demand3.m$demand, c(.1,.3,.5,.7,.9))
quantile(demand0c$demand, c(.1,.3,.5,.7,.9))
quantile(demand1c$demand, c(.1,.3,.5,.7,.9))
quantile(demand3c.m$demand, c(.1,.3,.5,.7,.9))


#inventory candidates
inv.range <- 21:50

#calculate realized profit
profit0 <- demand0
profit1 <- demand1
profit3.m <- demand3.m
profit0c <- demand0c
profit1c <- demand1c
profit3c.m <- demand3c.m

for (inv in inv.range) {
    profit0[, paste0("profit.",inv)] <-  Price * pmin(profit0$demand, inv) - Cost * inv
    profit1[, paste0("profit.",inv)] <-  Price * pmin(profit1$demand, inv) - Cost * inv
    profit3.m[, paste0("profit.",inv)] <-  Price * pmin(profit3.m$demand, inv) - Cost * inv
    profit0c[, paste0("profit.",inv)] <-  Price * pmin(profit0c$demand, inv) - Cost * inv
    profit1c[, paste0("profit.",inv)] <-  Price * pmin(profit1c$demand, inv) - Cost * inv
    profit3c.m[, paste0("profit.",inv)] <-  Price * pmin(profit3c.m$demand, inv) - Cost * inv
}


head(profit0)



# conditional mean and variance (conditioning on elambda)
cmean0 <- ddply(profit0, .(id), colMeans)
cmean1 <- ddply(profit1, .(id), colMeans)
cmean3.m <- ddply(profit3.m, .(id), colMeans)
cmean0c <- ddply(profit0c, .(id), colMeans)
cmean1c <- ddply(profit1c, .(id), colMeans)
cmean3c.m <- ddply(profit3c.m, .(id), colMeans)

cvar0 <- ddply(profit0, .(id), function(x) {apply(x, 2, var)})
cvar1 <- ddply(profit1, .(id), function(x) {apply(x, 2, var)})
cvar3.m <- ddply(profit3.m, .(id), function(x) {apply(x, 2, var)})
cvar0c <- ddply(profit0c, .(id), function(x) {apply(x, 2, var)})
cvar1c <- ddply(profit1c, .(id), function(x) {apply(x, 2, var)})
cvar3c.m <- ddply(profit3c.m, .(id), function(x) {apply(x, 2, var)})

#plot cmean distribution
ggplot() + 
    geom_density(aes(x=cmean0$profit.27), color="black",) + 
    geom_density(aes(x=cmean1$profit.27), color="grey") +
    geom_density(aes(x=cmean3.m$profit.27), color="blue") 
#plot cvar distribution
ggplot() + 
    geom_density(aes(x=cvar0$profit.27), color="black",) + 
    geom_density(aes(x=cvar1$profit.27), color="grey") +
    geom_density(aes(x=cvar3.m$profit.27), color="blue") 




# total mean and var
mean0 <- colMeans(profit0)
mean1 <- colMeans(profit1)
mean3.m <- colMeans(profit3.m)
mean0c <- colMeans(profit0c)
mean1c <- colMeans(profit1c)
mean3c.m <- colMeans(profit3c.m)


var0 <- apply(profit0, 2, var)
var1 <- apply(profit1, 2, var)
var3.m <- apply(profit3.m, 2, var)
var0c <- apply(profit0c, 2, var)
var1c <- apply(profit1c, 2, var)
var3c.m <- apply(profit3c.m, 2, var)



#plot efficient frontier of total mean-variance
plot(mean0[3:16], var0[3:16], type="l")
lines(mean1[3:16], var1[3:16], type="l", col="grey")
lines(mean3.m[3:16], var3.m[3:16], type="l", col="red")

plot(mean0c[3:16], var0c[3:16], type="l")
lines(mean1c[3:16], var1c[3:16], type="l", col="darkgrey")
lines(mean3c.m[3:16], var3c.m[3:16], type="l", col="blue")


#expectation-maximization inv
inv.range[1] - 1 + which.max(mean0[-(1:2)])
inv.range[1] - 1 +  which.max(mean1[-(1:2)])
inv.range[1] - 1 +  which.max(mean3.m[-(1:2)])
inv.range[1] - 1 + which.max(mean0c[-(1:2)])
inv.range[1] - 1 +  which.max(mean1c[-(1:2)])
inv.range[1] - 1 +  which.max(mean3c.m[-(1:2)])

