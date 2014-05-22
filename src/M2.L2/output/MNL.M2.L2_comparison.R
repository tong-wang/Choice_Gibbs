require("ggplot2")

setwd("~/Dropbox/RCode/Choice_Gibbs.git/src/M2.L2/output/K90L50")

load("MNL.M2.L2_Scenario0.RData")
load("MNL.M2.L2_Scenario1.RData")
load("MNL.M2.L2_Scenario3M.xl.RData")
load("MNL.M2.L2_Scenario3M.l.RData")
load("MNL.M2.L2_Scenario3M.m.RData")
load("MNL.M2.L2_Scenario3M.h.RData")
load("MNL.M2.L2_Scenario3M.xh.RData")



nrun <- 50000
burnin <- 0.5
start <- burnin*nrun+1




### Visualize results

#plot lambda (true lambda=50)
lam <- data.frame(
    lam0=z0$lambdas[start:nrun],
    lam1=z1$lambdas[start:nrun,], 
    lam3.xl=z3M.xl$lambdas[start:nrun,], 
    lam3.l=z3M.l$lambdas[start:nrun,], 
    lam3.m=z3M.m$lambdas[start:nrun,], 
    lam3.h=z3M.h$lambdas[start:nrun,],
    lam3.xh=z3M.xh$lambdas[start:nrun,]
)

colMeans(lam)



plot(lam$lam0, type="l")
plot(lam$lam1, type="l")
plot(lam$lam3.xl, type="l")
plot(lam$lam3.l, type="l")
plot(lam$lam3.m, type="l")
plot(lam$lam3.h, type="l")
plot(lam$lam3.xh, type="l")

ggplot(data=lam) + geom_density(aes(x=lam0), color="black") + 
    geom_density(aes(x=lam1), color="grey") +
    geom_density(aes(x=lam3.xl), color="red") +
    geom_density(aes(x=lam3.l), color="green") +
    geom_density(aes(x=lam3.m), color="blue") +
    geom_density(aes(x=lam3.h), color="yellow") +
    geom_density(aes(x=lam3.xh), color="purple")




### Plot beta
beta <- data.frame(
    beta0=z0$betas[start:nrun,],
    beta1=z1$betas[start:nrun,], 
    beta3.xl=z3M.xl$betas[start:nrun,], 
    beta3.l=z3M.l$betas[start:nrun,], 
    beta3.m=z3M.m$betas[start:nrun,], 
    beta3.h=z3M.h$betas[start:nrun,],
    beta3.xh=z3M.xh$betas[start:nrun,]
)

colMeans(beta)


#plot beta.coef1 (-2)
plot(beta$beta0.1, type="l")
plot(beta$beta1.1, type="l")
plot(beta$beta3.xl.1, type="l")
plot(beta$beta3.l.1, type="l")
plot(beta$beta3.m.1, type="l")
plot(beta$beta3.h.1, type="l")
plot(beta$beta3.xh.1, type="l")

ggplot(data=beta) + geom_density(aes(x=beta0.1), color="black") + 
    geom_density(aes(x=beta1.1), color="grey") +
    geom_density(aes(x=beta3.xl.1), color="red") +
    geom_density(aes(x=beta3.l.1), color="green") +
    geom_density(aes(x=beta3.m.1), color="blue") +
    geom_density(aes(x=beta3.h.1), color="yellow") +
    geom_density(aes(x=beta3.xh.1), color="purple")




#plot beta.coef2 (0.5)
plot(beta$beta0.2, type="l")
plot(beta$beta1.2, type="l")
plot(beta$beta3.xl.2, type="l")
plot(beta$beta3.l.2, type="l")
plot(beta$beta3.m.2, type="l")
plot(beta$beta3.h.2, type="l")
plot(beta$beta3.xh.2, type="l")

ggplot(data=beta) + geom_density(aes(x=beta0.2), color="black") + 
    geom_density(aes(x=beta1.2), color="grey") +
    geom_density(aes(x=beta3.xl.2), color="red") +
    geom_density(aes(x=beta3.l.2), color="green") +
    geom_density(aes(x=beta3.m.2), color="blue") +
    geom_density(aes(x=beta3.h.2), color="yellow") +
    geom_density(aes(x=beta3.xh.2), color="purple")



plot(beta$beta1.3, type="l")
plot(beta$beta3.xl.3, type="l")
plot(beta$beta3.l.3, type="l")
plot(beta$beta3.m.3, type="l")
plot(beta$beta3.h.3, type="l")
plot(beta$beta3.xh.3, type="l")

ggplot(data=beta) + geom_density(aes(x=beta0.3), color="black") + 
    geom_density(aes(x=beta1.3), color="grey") +
    geom_density(aes(x=beta3.xl.3), color="red") +
    geom_density(aes(x=beta3.l.3), color="green") +
    geom_density(aes(x=beta3.m.3), color="blue") +
    geom_density(aes(x=beta3.h.3), color="yellow") +
    geom_density(aes(x=beta3.xh.3), color="purple")


#plot beta.const2 (4)
plot(beta$beta0.4, type="l")
plot(beta$beta1.4, type="l")
plot(beta$beta3.xl.4, type="l")
plot(beta$beta3.l.4, type="l")
plot(beta$beta3.m.4, type="l")
plot(beta$beta3.h.4, type="l")
plot(beta$beta3.xh.4, type="l")

ggplot(data=beta) + geom_density(aes(x=beta0.4), color="black") + 
    geom_density(aes(x=beta1.4), color="grey") +
    geom_density(aes(x=beta3.xl.4), color="red") +
    geom_density(aes(x=beta3.l.4), color="green") +
    geom_density(aes(x=beta3.m.4), color="blue") +
    geom_density(aes(x=beta3.h.4), color="yellow") +
    geom_density(aes(x=beta3.xh.4), color="purple")
