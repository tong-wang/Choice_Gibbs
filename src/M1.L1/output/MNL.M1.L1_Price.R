####
# Price optimization
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

#plot the posteriors
ggplot(data=posteriors) + geom_density(aes(x=lam0), color="black") + 
    geom_density(aes(x=lam1), color="grey") +
    geom_density(aes(x=lam3.m), color="blue") 
ggplot(data=posteriors) + geom_density(aes(x=beta0.1), color="black") + 
    geom_density(aes(x=beta1.1), color="grey") +
    geom_density(aes(x=beta3.m.1), color="blue")
ggplot(data=posteriors) + geom_density(aes(x=beta0.2), color="black") + 
    geom_density(aes(x=beta1.2), color="grey") +
    geom_density(aes(x=beta3.m.2), color="blue")

ggplot(data=posteriors) + geom_density(aes(x=lam0c), color="black") + 
    geom_density(aes(x=lam1c), color="grey") +
    geom_density(aes(x=lam3c.m), color="blue") 
ggplot(data=posteriors) + geom_density(aes(x=beta0c.1), color="black") + 
    geom_density(aes(x=beta1c.1), color="grey") +
    geom_density(aes(x=beta3c.m.1), color="blue")
ggplot(data=posteriors) + geom_density(aes(x=beta0c.2), color="black") + 
    geom_density(aes(x=beta1c.2), color="grey") +
    geom_density(aes(x=beta3c.m.2), color="blue")




#price candidates
Cost <- 1
price.range <- seq(1, 3, 0.025) * Cost
price.index <- 1:length(price.range)



## generate profit samples
num <- 100

generate.profit <- function(row, scenario, n) {
    
    out <- data.frame(id=rep(row$id, n))
    for (i in price.index) {
        out[, paste0("profit.",i)] <- (price.range[i] - Cost) * rpois(n, row[, paste0("lam", scenario)]  / (1 + exp(- price.range[i] * row[, paste0("beta",scenario, ".1")] - row[, paste0("beta",scenario, ".2")])))
    }
    
    return(out)
}

profit0 <- ddply(posteriors, .(id), generate.profit, scenario="0", n=num)
profit1 <- ddply(posteriors, .(id), generate.profit, scenario="1", n=num)
profit3.m <- ddply(posteriors, .(id), generate.profit, scenario="3.m", n=num)
profit0c <- ddply(posteriors, .(id), generate.profit, scenario="0c", n=num)
profit1c <- ddply(posteriors, .(id), generate.profit, scenario="1c", n=num)
profit3c.m <- ddply(posteriors, .(id), generate.profit, scenario="3c.m", n=num)

head(profit0)


#plot profit distribution for a given price level
ggplot() + 
    geom_density(aes(x=profit0$profit.10), color="black",) + 
    geom_density(aes(x=profit1$profit.10), color="grey") +
    geom_density(aes(x=profit3.m$profit.10), color="blue") +
    geom_density(aes(x=profit0c$profit.10), color="black", linetype="dashed") + 
    geom_density(aes(x=profit1c$profit.10), color="grey", linetype="dashed") +
    geom_density(aes(x=profit3c.m$profit.10), color="blue", linetype="dashed")




####### add estimates from EM methods ##########

load("MNL.M1.L1_Scenario0.EM.RData")
#load("MNL.M1.L1_Scenario1c.EM.RData")
load("MNL.M1.L1_Scenario3Mc.m.EM.RData")

lambda0.EM <- tail(z0.EM$lambdas, n=1)
beta0.EM <- tail(z0.EM$betas, n=1)

lambda3c.m.EM <- tail(z3Mc.m.EM$lambdas, n=1)
beta3c.m.EM <- tail(z3Mc.m.EM$betas, n=1)


num.EM <- 100000

profit0.EM <- data.frame(id=rep(1, num.EM))
profit3c.m.EM <- data.frame(id=rep(1, num.EM))
for (i in price.index) {
    profit0.EM[, paste0("profit.",i)] <- (price.range[i] - Cost) * rpois(num.EM, lambda0.EM  / (1 + exp(- price.range[i] * beta0.EM[1] - beta0.EM[2])))
    profit3c.m.EM[, paste0("profit.",i)] <- (price.range[i] - Cost) * rpois(num.EM, lambda3c.m.EM  / (1 + exp(- price.range[i] * beta3c.m.EM[1] - beta3c.m.EM[2])))
}

################################################




# conditional mean and variance (conditioning on elambda)
cmean0 <- ddply(profit0, .(id), colMeans)[,-1]
cmean1 <- ddply(profit1, .(id), colMeans)[,-1]
cmean3.m <- ddply(profit3.m, .(id), colMeans)[,-1]
cmean0c <- ddply(profit0c, .(id), colMeans)[,-1]
cmean1c <- ddply(profit1c, .(id), colMeans)[,-1]
cmean3c.m <- ddply(profit3c.m, .(id), colMeans)[,-1]
cmean0.EM <- ddply(profit0.EM, .(id), colMeans)[,-1]
cmean3c.m.EM <- ddply(profit3c.m.EM, .(id), colMeans)[,-1]

cvar0 <- ddply(profit0, .(id), function(x) {apply(x, 2, var)})[,-1]
cvar1 <- ddply(profit1, .(id), function(x) {apply(x, 2, var)})[,-1]
cvar3.m <- ddply(profit3.m, .(id), function(x) {apply(x, 2, var)})[,-1]
cvar0c <- ddply(profit0c, .(id), function(x) {apply(x, 2, var)})[,-1]
cvar1c <- ddply(profit1c, .(id), function(x) {apply(x, 2, var)})[,-1]
cvar3c.m <- ddply(profit3c.m, .(id), function(x) {apply(x, 2, var)})[,-1]


var.cmean0 <- apply(cmean0, 2, var)
var.cmean1 <- apply(cmean1, 2, var)
var.cmean3.m <- apply(cmean3.m, 2, var)
var.cmean0c <- apply(cmean0c, 2, var)
var.cmean1c <- apply(cmean1c, 2, var)
var.cmean3c.m <- apply(cmean3c.m, 2, var)


#plot cmean distribution
ggplot() + 
    geom_density(aes(x=cmean0$profit.10), color="black",) + 
    geom_density(aes(x=cmean1$profit.10), color="grey") +
    geom_density(aes(x=cmean3.m$profit.10), color="blue") 
ggplot() + 
    geom_density(aes(x=cmean0c$profit.10), color="black",) + 
    geom_density(aes(x=cmean1c$profit.10), color="grey") +
    geom_density(aes(x=cmean3c.m$profit.10), color="blue") 
#plot cvar distribution
ggplot() + 
    geom_density(aes(x=cvar0$profit.10), color="black",) + 
    geom_density(aes(x=cvar1$profit.10), color="grey") +
    geom_density(aes(x=cvar3.m$profit.10), color="blue") 



# total mean, var, quantile
mean0 <- colMeans(cmean0)
mean1 <- colMeans(cmean1)
mean3.m <- colMeans(cmean3.m)
mean0c <- colMeans(cmean0c)
mean1c <- colMeans(cmean1c)
mean3c.m <- colMeans(cmean3c.m)
mean0.EM <- colMeans(cmean0.EM)
mean3c.m.EM <- colMeans(cmean3c.m.EM)

var0 <- apply(profit0, 2, var)[-1]
var1 <- apply(profit1, 2, var)[-1]
var3.m <- apply(profit3.m, 2, var)[-1]
var0c <- apply(profit0c, 2, var)[-1]
var1c <- apply(profit1c, 2, var)[-1]
var3c.m <- apply(profit3c.m, 2, var)[-1]
var0.EM <- apply(profit0.EM, 2, var)[-1]
var3c.m.EM <- apply(profit3c.m.EM, 2, var)[-1]



#plot efficient frontier of total mean-variance
ggplot() +
    geom_path(aes(x=mean0, y=var0)) + 
    #geom_path(aes(x=mean1, y=var1), color="grey") + 
    #geom_path(aes(x=mean3.m, y=var3.m), color="red") + 
    #geom_path(aes(x=mean0c, y=var0c), color="yellow") + 
    geom_path(aes(x=mean1c, y=var1c), color="darkgrey") + 
    geom_path(aes(x=mean3c.m, y=var3c.m), color="blue") +
    #geom_path(aes(x=mean0.EM, y=var0.EM), color="green") 
    geom_path(aes(x=mean3c.m.EM, y=var3c.m.EM), color="green") 
    


#optimal price (risk neutral)
price.range[which.max(mean0)]
price.range[which.max(mean1)]
price.range[which.max(mean3.m)]
price.range[which.max(mean0c)]
price.range[which.max(mean1c)]
price.range[which.max(mean3c.m)]
price.range[which.max(mean3c.m.EM)]

#optimal expected profit
max(mean0)
max(mean1)
max(mean3.m)
max(mean0c)
max(mean1c)
max(mean3c.m)
max(mean3c.m.EM)



#optimal price (risk averse: mean >= b*var)
b <- 1
inv.range[which.max(mean0[mean0>=b*var0])]
inv.range[which.max(mean1c[mean1c>=b*var1c])]
inv.range[which.max(mean3c.m[mean3c.m>=b*var3c.m])]
inv.range[which.max(mean3c.m.EM[mean3c.m.EM>=b*var3c.m.EM])]
