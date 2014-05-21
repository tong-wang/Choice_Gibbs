setwd("~/Dropbox/RCode/Choice_Gibbs.git/src/binary.L2/output/K90L50")

load("MNL.binary.L2_Scenario0.RData")
load("MNL.binary.L2_Scenario1.RData")
load("MNL.binary.L2_Scenario3M.xl.RData")
load("MNL.binary.L2_Scenario3M.l.RData")
load("MNL.binary.L2_Scenario3M.m.RData")
load("MNL.binary.L2_Scenario3M.h.RData")
load("MNL.binary.L2_Scenario3M.xh.RData")


load("MNL.binary.L2_Scenario0c.RData")
load("MNL.binary.L2_Scenario1c.RData")
load("MNL.binary.L2_Scenario3Mc.xl.RData")
load("MNL.binary.L2_Scenario3Mc.l.RData")
load("MNL.binary.L2_Scenario3Mc.m.RData")
load("MNL.binary.L2_Scenario3Mc.h.RData")
load("MNL.binary.L2_Scenario3Mc.xh.RData")


nrun <- 50000
burnin <- 0.5
start <- burnin*nrun+1






### Visualize results
require(ggplot2)


#plot lambda (true lambda=50)
lam <- data.frame(
    lam0=z0$lambdas[start:nrun],
    lam1=z1$lambdas[start:nrun,], 
    lam3.xl=z3M.xl$lambdas[start:nrun,], 
    lam3.l=z3M.l$lambdas[start:nrun,], 
    lam3.m=z3M.m$lambdas[start:nrun,], 
    lam3.h=z3M.h$lambdas[start:nrun,],
    lam3.xh=z3M.xh$lambdas[start:nrun,],
    lam0c=z0c$lambdas[start:nrun],
    lam1c=z1c$lambdas[start:nrun,], 
    lam3c.xl=z3Mc.xl$lambdas[start:nrun,], 
    lam3c.l=z3Mc.l$lambdas[start:nrun,], 
    lam3c.m=z3Mc.m$lambdas[start:nrun,], 
    lam3c.h=z3Mc.h$lambdas[start:nrun,],
    lam3c.xh=z3Mc.xh$lambdas[start:nrun,]
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
    beta3.xl=z3M.xl$betas[start:nrun,], 
    beta3.l=z3M.l$betas[start:nrun,], 
    beta3.m=z3M.m$betas[start:nrun,], 
    beta3.h=z3M.h$betas[start:nrun,],
    beta3.xh=z3M.xh$betas[start:nrun,],
    beta0c=z0c$betas[start:nrun,],
    beta1c=z1c$betas[start:nrun,], 
    beta3c.xl=z3Mc.xl$betas[start:nrun,], 
    beta3c.l=z3Mc.l$betas[start:nrun,], 
    beta3c.m=z3Mc.m$betas[start:nrun,], 
    beta3c.h=z3Mc.h$betas[start:nrun,],
    beta3c.xh=z3Mc.xh$betas[start:nrun,]
)

colMeans(beta)


#plot beta1 (-3)
plot(beta$beta0.1, type="l")
plot(beta$beta1.1, type="l")
plot(beta$beta3.xl.1, type="l")
plot(beta$beta3.l.1, type="l")
plot(beta$beta3.m.1, type="l")
plot(beta$beta3.h.1, type="l")
plot(beta$beta3.xh.1, type="l")

plot(beta$beta0c.1, type="l")
plot(beta$beta1c.1, type="l")
plot(beta$beta3c.xl.1, type="l")
plot(beta$beta3c.l.1, type="l")
plot(beta$beta3c.m.1, type="l")
plot(beta$beta3c.h.1, type="l")
plot(beta$beta3c.xh.1, type="l")


#non censored
ggplot(data=beta) + geom_density(aes(x=beta0.1), color="black") + 
    geom_density(aes(x=beta1.1), color="grey") +
    geom_density(aes(x=beta3.xl.1), color="red") +
    geom_density(aes(x=beta3.l.1), color="green") +
    geom_density(aes(x=beta3.m.1), color="blue") +
    geom_density(aes(x=beta3.h.1), color="yellow") +
    geom_density(aes(x=beta3.xh.1), color="purple") #+ #
#censored
ggplot(data=beta) + 
    geom_density(aes(x=beta0c.1), color="black", linetype="dashed") + 
    geom_density(aes(x=beta1c.1), color="grey", linetype="dashed") +
    geom_density(aes(x=beta3c.xl.1), color="red", linetype="dashed") +
    geom_density(aes(x=beta3c.l.1), color="green", linetype="dashed") +
    geom_density(aes(x=beta3c.m.1), color="blue", linetype="dashed") +
    geom_density(aes(x=beta3c.h.1), color="yellow", linetype="dashed") +
    geom_density(aes(x=beta3c.xh.1), color="purple", linetype="dashed")

#by scenario
ggplot(data=beta) + geom_density(aes(x=beta0.1), color="black") + geom_density(aes(x=beta0c.1), color="black", linetype="dashed")
ggplot(data=beta) + geom_density(aes(x=beta1.1), color="grey") + geom_density(aes(x=beta1c.1), color="grey", linetype="dashed")
ggplot(data=beta) + geom_density(aes(x=beta3.xl.1), color="red") + geom_density(aes(x=beta3c.xl.1), color="red", linetype="dashed")
ggplot(data=beta) + geom_density(aes(x=beta3.l.1), color="green") + geom_density(aes(x=beta3c.l.1), color="green", linetype="dashed")
ggplot(data=beta) + geom_density(aes(x=beta3.m.1), color="blue") + geom_density(aes(x=beta3c.m.1), color="blue", linetype="dashed")
ggplot(data=beta) + geom_density(aes(x=beta3.h.1), color="yellow") + geom_density(aes(x=beta3c.h.1), color="yellow", linetype="dashed")
ggplot(data=beta) + geom_density(aes(x=beta3.xh.1), color="purple") + geom_density(aes(x=beta3c.xh.1), color="purple", linetype="dashed")



#plot betaL (6)
plot(beta$beta0.2, type="l")
plot(beta$beta1.2, type="l")
plot(beta$beta3.xl.2, type="l")
plot(beta$beta3.l.2, type="l")
plot(beta$beta3.m.2, type="l")
plot(beta$beta3.h.2, type="l")
plot(beta$beta3.xh.2, type="l")

plot(beta$beta0c.2, type="l")
plot(beta$beta1c.2, type="l")
plot(beta$beta3c.xl.2, type="l")
plot(beta$beta3c.l.2, type="l")
plot(beta$beta3c.m.2, type="l")
plot(beta$beta3c.h.2, type="l")
plot(beta$beta3c.xh.2, type="l")


#non censored
ggplot(data=beta) + geom_density(aes(x=beta0.2), color="black") + 
    geom_density(aes(x=beta1.2), color="grey") +
    geom_density(aes(x=beta3.xl.2), color="red") +
    geom_density(aes(x=beta3.l.2), color="green") +
    geom_density(aes(x=beta3.m.2), color="blue") +
    geom_density(aes(x=beta3.h.2), color="yellow") +
    geom_density(aes(x=beta3.xh.2), color="purple") #+ #
#censored
ggplot(data=beta) + 
    geom_density(aes(x=beta0c.2), color="black", linetype="dashed") + 
    geom_density(aes(x=beta1c.2), color="grey", linetype="dashed") +
    geom_density(aes(x=beta3c.xl.2), color="red", linetype="dashed") +
    geom_density(aes(x=beta3c.l.2), color="green", linetype="dashed") +
    geom_density(aes(x=beta3c.m.2), color="blue", linetype="dashed") +
    geom_density(aes(x=beta3c.h.2), color="yellow", linetype="dashed") +
    geom_density(aes(x=beta3c.xh.2), color="purple", linetype="dashed")

#by scenario
ggplot(data=beta) + geom_density(aes(x=beta0.2), color="black") + geom_density(aes(x=beta0c.2), color="black", linetype="dashed")
ggplot(data=beta) + geom_density(aes(x=beta1.2), color="grey") + geom_density(aes(x=beta1c.2), color="grey", linetype="dashed")
ggplot(data=beta) + geom_density(aes(x=beta3.xl.2), color="red") + geom_density(aes(x=beta3c.xl.2), color="red", linetype="dashed")
ggplot(data=beta) + geom_density(aes(x=beta3.l.2), color="green") + geom_density(aes(x=beta3c.l.2), color="green", linetype="dashed")
ggplot(data=beta) + geom_density(aes(x=beta3.m.2), color="blue") + geom_density(aes(x=beta3c.m.2), color="blue", linetype="dashed")
ggplot(data=beta) + geom_density(aes(x=beta3.h.2), color="yellow") + geom_density(aes(x=beta3c.h.2), color="yellow", linetype="dashed")
ggplot(data=beta) + geom_density(aes(x=beta3.xh.2), color="purple") + geom_density(aes(x=beta3c.xh.2), color="purple", linetype="dashed")

