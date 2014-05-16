setwd("~/Dropbox/RCode/Choice_Gibbs.git/src/binary.L1/output/K90L50Beta0.3")

load("MNL_Scenario0.binary.L1.RData")
load("MNL_Scenario1.binary.L1.RData")
load("MNL_Scenario3_M.binary.L1.xl.RData")
load("MNL_Scenario3_M.binary.L1.l.RData")
load("MNL_Scenario3_M.binary.L1.m.RData")
load("MNL_Scenario3_M.binary.L1.h.RData")
load("MNL_Scenario3_M.binary.L1.xh.RData")


load("MNL_Scenario0c.binary.L1.RData")
load("MNL_Scenario1c.binary.L1.RData")
load("MNL_Scenario3c_M.binary.L1.xl.RData")
load("MNL_Scenario3c_M.binary.L1.l.RData")
load("MNL_Scenario3c_M.binary.L1.m.RData")
load("MNL_Scenario3c_M.binary.L1.h.RData")
load("MNL_Scenario3c_M.binary.L1.xh.RData")


nrun <- 50000
burnin <- 0.5
start <- burnin*nrun+1






### Visualize results
require(ggplot2)


#plot lambda (true lambda=50)
lam <- data.frame(
    lam0=z0$lambdas[start:nrun],
    lam1=z1$lambdas[start:nrun,], 
    lam3.xl=z3.xl$lambdas[start:nrun,], 
    lam3.l=z3.l$lambdas[start:nrun,], 
    lam3.m=z3.m$lambdas[start:nrun,], 
    lam3.h=z3.h$lambdas[start:nrun,],
    lam3.xh=z3.xh$lambdas[start:nrun,],
    lam0c=z0c$lambdas[start:nrun],
    lam1c=z1c$lambdas[start:nrun,], 
    lam3c.xl=z3c.xl$lambdas[start:nrun,], 
    lam3c.l=z3c.l$lambdas[start:nrun,], 
    lam3c.m=z3c.m$lambdas[start:nrun,], 
    lam3c.h=z3c.h$lambdas[start:nrun,],
    lam3c.xh=z3c.xh$lambdas[start:nrun,]
)

colMeans(lam)



plot(lam$lam0, type="l")
plot(lam$lam1, type="l")
plot(lam$lam3.xl, type="l")
plot(lam$lam3.l, type="l")
plot(lam$lam3.m, type="l")
plot(lam$lam3.h, type="l")
plot(lam$lam3.xh, type="l")

plot(lam$lam0c, type="l")
plot(lam$lam1c, type="l")
plot(lam$lam3c.xl, type="l")
plot(lam$lam3c.l, type="l")
plot(lam$lam3c.m, type="l")
plot(lam$lam3c.h, type="l")
plot(lam$lam3c.xh, type="l")


#non censored
ggplot(data=lam) + geom_density(aes(x=lam0), color="black") + 
    geom_density(aes(x=lam1), color="grey") +
    geom_density(aes(x=lam3.xl), color="red") +
    geom_density(aes(x=lam3.l), color="green") +
    geom_density(aes(x=lam3.m), color="blue") +
    geom_density(aes(x=lam3.h), color="yellow") +
    geom_density(aes(x=lam3.xh), color="purple") #+ #
#censored
ggplot(data=lam) + 
    geom_density(aes(x=lam0c), color="black", linetype="dashed") + 
    geom_density(aes(x=lam1c), color="grey", linetype="dashed") +
    geom_density(aes(x=lam3c.xl), color="red", linetype="dashed") +
    geom_density(aes(x=lam3c.l), color="green", linetype="dashed") +
    geom_density(aes(x=lam3c.m), color="blue", linetype="dashed") +
    geom_density(aes(x=lam3c.h), color="yellow", linetype="dashed") +
    geom_density(aes(x=lam3c.xh), color="purple", linetype="dashed")


#by scenario
ggplot(data=lam) + geom_density(aes(x=lam0), color="black") + geom_density(aes(x=lam0c), color="grey")
ggplot(data=lam) + geom_density(aes(x=lam1), color="black") + geom_density(aes(x=lam1c), color="grey")
ggplot(data=lam) + geom_density(aes(x=lam3.xl), color="black") + geom_density(aes(x=lam3c.xl), color="grey")
ggplot(data=lam) + geom_density(aes(x=lam3.l), color="black") + geom_density(aes(x=lam3c.l), color="grey")
ggplot(data=lam) + geom_density(aes(x=lam3.m), color="black") + geom_density(aes(x=lam3c.m), color="grey")
ggplot(data=lam) + geom_density(aes(x=lam3.h), color="black") + geom_density(aes(x=lam3c.h), color="grey")
ggplot(data=lam) + geom_density(aes(x=lam3.xh), color="black") + geom_density(aes(x=lam3c.xh), color="grey")


### Plot beta
beta <- data.frame(
    beta0=z0$betas[start:nrun,],
    beta1=z1$betas[start:nrun,], 
    beta3.xl=z3.xl$betas[start:nrun,], 
    beta3.l=z3.l$betas[start:nrun,], 
    beta3.m=z3.m$betas[start:nrun,], 
    beta3.h=z3.h$betas[start:nrun,],
    beta3.xh=z3.xh$betas[start:nrun,],
    beta0c=z0c$betas[start:nrun,],
    beta1c=z1c$betas[start:nrun,], 
    beta3c.xl=z3c.xl$betas[start:nrun,], 
    beta3c.l=z3c.l$betas[start:nrun,], 
    beta3c.m=z3c.m$betas[start:nrun,], 
    beta3c.h=z3c.h$betas[start:nrun,],
    beta3c.xh=z3c.xh$betas[start:nrun,]
)

colMeans(beta)


#plot beta1 (0.3)
plot(beta$beta0, type="l")
plot(beta$beta1, type="l")
plot(beta$beta3.xl, type="l")
plot(beta$beta3.l, type="l")
plot(beta$beta3.m, type="l")
plot(beta$beta3.h, type="l")
plot(beta$beta3.xh, type="l")

plot(beta$beta0c, type="l")
plot(beta$beta1c, type="l")
plot(beta$beta3c.xl, type="l")
plot(beta$beta3c.l, type="l")
plot(beta$beta3c.m, type="l")
plot(beta$beta3c.h, type="l")
plot(beta$beta3c.xh, type="l")


#non censored
ggplot(data=beta) + geom_density(aes(x=beta0), color="black") + 
    geom_density(aes(x=beta1), color="grey") +
    geom_density(aes(x=beta3.xl), color="red") +
    geom_density(aes(x=beta3.l), color="green") +
    geom_density(aes(x=beta3.m), color="blue") +
    geom_density(aes(x=beta3.h), color="yellow") +
    geom_density(aes(x=beta3.xh), color="purple") #+ #
#censored
ggplot(data=beta) + 
    geom_density(aes(x=beta0c), color="black", linetype="dashed") + 
    geom_density(aes(x=beta1c), color="grey", linetype="dashed") +
    geom_density(aes(x=beta3c.xl), color="red", linetype="dashed") +
    geom_density(aes(x=beta3c.l), color="green", linetype="dashed") +
    geom_density(aes(x=beta3c.m), color="blue", linetype="dashed") +
    geom_density(aes(x=beta3c.h), color="yellow", linetype="dashed") +
    geom_density(aes(x=beta3c.xh), color="purple", linetype="dashed")

#by scenario
ggplot(data=beta) + geom_density(aes(x=beta0), color="black") + geom_density(aes(x=beta0c), color="black", linetype="dashed")
ggplot(data=beta) + geom_density(aes(x=beta1), color="grey") + geom_density(aes(x=beta1c), color="grey", linetype="dashed")
ggplot(data=beta) + geom_density(aes(x=beta3.xl), color="red") + geom_density(aes(x=beta3c.xl), color="red", linetype="dashed")
ggplot(data=beta) + geom_density(aes(x=beta3.l), color="green") + geom_density(aes(x=beta3c.l), color="green", linetype="dashed")
ggplot(data=beta) + geom_density(aes(x=beta3.m), color="blue") + geom_density(aes(x=beta3c.m), color="blue", linetype="dashed")
ggplot(data=beta) + geom_density(aes(x=beta3.h), color="yellow") + geom_density(aes(x=beta3c.h), color="yellow", linetype="dashed")
ggplot(data=beta) + geom_density(aes(x=beta3.xh), color="purple") + geom_density(aes(x=beta3c.xh), color="purple", linetype="dashed")

