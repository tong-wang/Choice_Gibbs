setwd("~/Dropbox/RCode/Choice_Gibbs.git/src/K360Lambda100/")

load("MNL_Scenario0.RData")

load("MNL_Scenario1.RData")

load("MNL_Scenario3_Ml.RData")
load("MNL_Scenario3_M.RData")
load("MNL_Scenario3_M1.RData")
load("MNL_Scenario3_M2.RData")
load("MNL_Scenario3_M3.RData")
load("MNL_Scenario3_Mh.RData")
load("MNL_Scenario3_M1h.RData")
load("MNL_Scenario3_M2h.RData")
load("MNL_Scenario3_M3h.RData")
load("MNL_Scenario3_Mxh.RData")


load("MNL_Scenario3_A.RData")
load("MNL_Scenario3_Al.RData")
load("MNL_Scenario3_Ah.RData")
load("MNL_Scenario3_Axh.RData")

load("MNL_Scenario3_Reg.RData")
load("MNL_Scenario3_Regl.RData")
load("MNL_Scenario3_Regh.RData")
load("MNL_Scenario3_Regxh.RData")



nrun <- 100000
burnin <- 0.5
start <- burnin*nrun+1






### Visualize results
require(ggplot2)




#plot lambda (true lambda=100)
lam <- data.frame(
    lam0=z0$lambdas,
    lam1=z1$lambdas[start:nrun,], 
    lam3ml=z3ml$lambdas[start:nrun,], 
    lam3m=z3m$lambdas[start:nrun,], 
    lam3m1=z3m1$lambdas[start:nrun,],
    lam3m2=z3m2$lambdas[start:nrun,],
    lam3m3=z3m3$lambdas[start:nrun,],
    lam3mh=z3mh$lambdas[start:nrun,],
    lam3mxh=z3mxh$lambdas[start:nrun,],
    lam3m1h=z3m1h$lambdas[start:nrun,],
    lam3m2h=z3m2h$lambdas[start:nrun,],
    lam3m3h=z3m3h$lambdas[start:nrun,]#,
    #lam3al=z3al$lambdas[start:nrun,], 
    #lam3a=z3a$lambdas[start:nrun,], 
    #lam3ah=z3ah$lambdas[start:nrun,],
    #lam3axh=z3axh$lambdas[start:nrun,]#,
    #lam3rl=z3rl$lambdas[start:nrun,], 
    #lam3r=z3r$lambdas[start:nrun,], 
    #lam3rh=z3rh$lambdas[start:nrun,],
    #lam3rxh=z3rxh$lambdas[start:nrun,]
)

colMeans(lam)



plot(lam$lam0, type="l")
plot(lam$lam1, type="l")
plot(z3ml$lambdas, type="l")
plot(z3m$lambdas, type="l")
plot(z3m1$lambdas, type="l")
plot(z3m2$lambdas, type="l")
plot(z3m3$lambdas, type="l")
plot(z3mh$lambdas, type="l")
plot(z3m1h$lambdas, type="l")
plot(z3m2h$lambdas, type="l")
plot(z3m3h$lambdas, type="l")
plot(z3mxh$lambdas, type="l")
plot(lam$lam3al, type="l")
plot(lam$lam3a, type="l")
plot(lam$lam3ah, type="l")
plot(lam$lam3axh, type="l")
plot(lam$lam3rl, type="l")
plot(lam$lam3r, type="l")
plot(lam$lam3rh, type="l")
plot(lam$lam3rxh, type="l")




ggplot(data=lam) + geom_density(aes(x=lam0), color="black") + 
    geom_density(aes(x=lam1), color="grey") +
    geom_density(aes(x=lam3ml), color="red") +
    geom_density(aes(x=lam3m), color="green") +
    geom_density(aes(x=lam3m1), color="grey10") +
    geom_density(aes(x=lam3m2), color="grey50") +
    geom_density(aes(x=lam3m3), color="grey80") +
    geom_density(aes(x=lam3mh), color="yellow") +
    geom_density(aes(x=lam3m1h), color="red") +
    geom_density(aes(x=lam3m2h), color="green") +
    geom_density(aes(x=lam3m3h), color="blue") +
    geom_density(aes(x=lam3mxh), color="purple")

    
ggplot(data=lam) + geom_density(aes(x=lam0), color="black") + 
    geom_density(aes(x=lam1), color="grey") +
    geom_density(aes(x=lam3al), color="red") +
    geom_density(aes(x=lam3a), color="green") +
    geom_density(aes(x=lam3ah), color="blue") +
    geom_density(aes(x=lam3axh), color="purple")


ggplot(data=lam) + geom_density(aes(x=lam0), color="black") + 
    geom_density(aes(x=lam1), color="grey") +
    geom_density(aes(x=lam3rl), color="red") +
    geom_density(aes(x=lam3r), color="green") +
    geom_density(aes(x=lam3rh), color="blue") +
    geom_density(aes(x=lam3rxh), color="purple")






### Plot beta

beta <- data.frame(
    beta0=z0$betas[,], 
    beta1=z1$betas[start:nrun,], 
    beta3ml=z3ml$betas[start:nrun,], 
    beta3m=z3m$betas[start:nrun,], 
    beta3m1=z3m1$betas[start:nrun,], 
    beta3m2=z3m2$betas[start:nrun,], 
    beta3m3=z3m3$betas[start:nrun,], 
    beta3mh=z3mh$betas[start:nrun,],
    beta3mxh=z3mxh$betas[start:nrun,],
    beta3m1h=z3m1h$betas[start:nrun,],
    beta3m2h=z3m2h$betas[start:nrun,],
    beta3m3h=z3m3h$betas[start:nrun,]#,
    #beta3al=z3al$betas[start:nrun,], 
    #beta3a=z3a$betas[start:nrun,], 
    #beta3ah=z3ah$betas[start:nrun,],
    #beta3axh=z3axh$betas[start:nrun,]#,
    #beta3rl=z3rl$betas[start:nrun,], 
    #beta3r=z3r$betas[start:nrun,], 
    #beta3rh=z3rh$betas[start:nrun,],
    #beta3rxh=z3rxh$betas[start:nrun,]
)


colMeans(beta)



#plot beta1 (0.06)
plot(beta$beta0.1, type="l")
plot(beta$beta1.1, type="l")
plot(beta$beta3ml.1, type="l")
plot(beta$beta3m.1, type="l")
plot(beta$beta3mh.1, type="l")
plot(beta$beta3m1h.1, type="l")
plot(beta$beta3m2h.1, type="l")
plot(beta$beta3m3h.1, type="l")
plot(beta$beta3mxh.1, type="l")
plot(beta$beta3a.1, type="l")
plot(beta$beta3rl.1, type="l")
plot(beta$beta3r.1, type="l")
plot(beta$beta3rh.1, type="l")
plot(beta$beta3rxh.1, type="l")

ggplot(data=beta) + geom_density(aes(x=beta0.1), color="black") +
    geom_density(aes(x=beta1.1), color="grey") +
    geom_density(aes(x=beta3ml.1), color="red") +
    geom_density(aes(x=beta3m.1), color="green") +
    geom_density(aes(x=beta3m1.1), color="red") +
    geom_density(aes(x=beta3m2.1), color="green") +
    geom_density(aes(x=beta3m3.1), color="blue") +
    geom_density(aes(x=beta3mh.1), color="yellow") +
    geom_density(aes(x=beta3m1h.1), color="red") +
    geom_density(aes(x=beta3m2h.1), color="green") +
    geom_density(aes(x=beta3m3h.1), color="blue") +
    geom_density(aes(x=beta3mxh.1), color="purple")

ggplot(data=beta) + geom_density(aes(x=beta0.1), color="black") +
    geom_density(aes(x=beta1.1), color="grey") +
    geom_density(aes(x=beta3al.1), color="red") +
    geom_density(aes(x=beta3a.1), color="green") +
    geom_density(aes(x=beta3ah.1), color="blue") +
    geom_density(aes(x=beta3axh.1), color="purple")


ggplot(data=beta) + geom_density(aes(x=beta0.1), color="black") +
    geom_density(aes(x=beta1.1), color="grey") +
    geom_density(aes(x=beta3rl.1), color="red") +
    geom_density(aes(x=beta3r.1), color="green") +
    geom_density(aes(x=beta3rh.1), color="blue") +
    geom_density(aes(x=beta3rxh.1), color="purple")



#plot beta2 (0.03)
plot(beta$beta0.2, type="l")
plot(beta$beta1.2, type="l")
plot(beta$beta3ml.2, type="l")
plot(beta$beta3m.2, type="l")
plot(beta$beta3mh.2, type="l")
plot(beta$beta3m1h.2, type="l")
plot(beta$beta3m2h.2, type="l")
plot(beta$beta3m3h.2, type="l")
plot(beta$beta3mxh.2, type="l")
plot(beta$beta3a.2, type="l")
plot(beta$beta3r.2, type="l")


ggplot(data=beta) + geom_density(aes(x=beta0.2), color="black") +
    geom_density(aes(x=beta1.2), color="grey") +
    geom_density(aes(x=beta3ml.2), color="red") +
    geom_density(aes(x=beta3m.2), color="green") +
    geom_density(aes(x=beta3m1.2), color="green") +
    geom_density(aes(x=beta3m2.2), color="green") +
    geom_density(aes(x=beta3m3.2), color="green") +
    geom_density(aes(x=beta3mh.2), color="yellow") +
    geom_density(aes(x=beta3m1h.2), color="red") +
    geom_density(aes(x=beta3m2h.2), color="green") +
    geom_density(aes(x=beta3m3h.2), color="blue") +
    geom_density(aes(x=beta3mxh.2), color="purple")

ggplot(data=beta) + geom_density(aes(x=beta0.2), color="black") +
    geom_density(aes(x=beta1.2), color="grey") +
    geom_density(aes(x=beta3al.2), color="red") +
    geom_density(aes(x=beta3a.2), color="green") +
    geom_density(aes(x=beta3ah.2), color="blue") +
    geom_density(aes(x=beta3axh.2), color="purple")


ggplot(data=beta) + geom_density(aes(x=beta0.2), color="black") +
    geom_density(aes(x=beta1.2), color="grey") +
    geom_density(aes(x=beta3rl.2), color="red") +
    geom_density(aes(x=beta3r.2), color="green") +
    geom_density(aes(x=beta3rh.2), color="blue") +
    geom_density(aes(x=beta3rxh.2), color="purple")




# plot eps1.mu and eps1.sd in the M model (3, (1.5, 0.5, 0.1))
eps1.mu <- data.frame(
    eps1.mul=z3ml$eps1.mus[start:nrun,], 
    eps1.mu=z3m$eps1.mus[start:nrun,], 
    eps1.muh=z3mh$eps1.mus[start:nrun,],
    eps1.muxh=z3mxh$eps1.mus[start:nrun,]
)
eps1.sd <- data.frame(
    eps1.sdl=z3ml$eps1.sds[start:nrun,], 
    eps1.sd=z3m$eps1.sds[start:nrun,], 
    eps1.sdh=z3mh$eps1.sds[start:nrun,],
    eps1.sdxh=z3mxh$eps1.sds[start:nrun,]
)


colMeans(eps1.mu)
colMeans(eps1.sd)

plot(eps1.mu$eps1.mul, type="l")
plot(eps1.mu$eps1.mu, type="l")
plot(eps1.mu$eps1.muh, type="l")
plot(eps1.mu$eps1.muxh, type="l")

plot(eps1.sd$eps1.sdl, type="l")
plot(eps1.sd$eps1.sd, type="l")
plot(eps1.sd$eps1.sdh, type="l")
plot(eps1.sd$eps1.sdxh, type="l")


ggplot(data=eps1.mu) +
    geom_density(aes(x=eps1.mul), color="red") +
    geom_density(aes(x=eps1.mu), color="green") +
    geom_density(aes(x=eps1.muh), color="blue") #+ 
    #geom_density(aes(x=eps1.muxh), color="purple")

ggplot(data=eps1.sd) +
    geom_density(aes(x=eps1.sdl), color="red") +
    geom_density(aes(x=eps1.sd), color="green") +
    geom_density(aes(x=eps1.sdh), color="blue")





# plot eps2.mu and eps2.sd in the A model (5, (3, 1, 0.1))
eps2.mu <- data.frame(
    eps2.mul=z3al$eps2.mus[start:nrun,], 
    eps2.mu=z3a$eps2.mus[start:nrun,], 
    eps2.muh=z3ah$eps2.mus[start:nrun,],
    eps2.muxh=z3axh$eps2.mus[start:nrun,]
)
eps2.sd <- data.frame(
    eps2.sdl=z3al$eps2.sds[start:nrun,], 
    eps2.sd=z3a$eps2.sds[start:nrun,], 
    eps2.sdh=z3ah$eps2.sds[start:nrun,],
    eps2.sdxh=z3axh$eps2.sds[start:nrun,]
)


colMeans(eps2.mu)
colMeans(eps2.sd)

plot(eps2.mu$eps2.mul, type="l")
plot(eps2.mu$eps2.mu, type="l")
plot(eps2.mu$eps2.muh, type="l")
plot(eps2.mu$eps2.muxh, type="l")

plot(eps2.sd$eps2.sdl, type="l")
plot(eps2.sd$eps2.sd, type="l")
plot(eps2.sd$eps2.sdh, type="l")
plot(eps2.sd$eps2.sdxh, type="l")

ggplot(data=eps2.mu) +
    geom_density(aes(x=eps2.mul), color="red") +
    geom_density(aes(x=eps2.mu), color="green") +
    geom_density(aes(x=eps2.muh), color="blue")

ggplot(data=eps2.sd) +
    geom_density(aes(x=eps2.sdl), color="red") +
    geom_density(aes(x=eps2.sd), color="green") +
    geom_density(aes(x=eps2.sdh), color="blue")




# plot t.b and eps.sd in the Reg model (500, 5, (250, 150, 50))
t.b <- data.frame(
    t.bl=z3rl$t.bs[start:nrun,], 
    t.b=z3r$t.bs[start:nrun,], 
    t.bh=z3rh$t.bs[start:nrun,],
    t.bxh=z3rxh$t.bs[start:nrun,]
)
eps.sd <- data.frame(
    eps.sdl=z3rl$eps.sds[start:nrun,], 
    eps.sd=z3r$eps.sds[start:nrun,], 
    eps.sdh=z3rh$eps.sds[start:nrun,],
    eps.sdxh=z3rxh$eps.sds[start:nrun,]
)


colMeans(t.b)
colMeans(eps.sd)


plot(t.b$t.bl.1, type="l")
plot(t.b$t.b.1, type="l")
plot(t.b$t.bh.1, type="l")
plot(t.b$t.bxh.1, type="l")

plot(t.b$t.bl.2, type="l")
plot(t.b$t.b.2, type="l")
plot(t.b$t.bh.2, type="l")
plot(t.b$t.bxh.2, type="l")


ggplot(data=t.b) +
    geom_density(aes(x=t.bl.1), color="red") +
    geom_density(aes(x=t.b.1), color="green") +
    geom_density(aes(x=t.bh.1), color="blue") #+
    #geom_density(aes(x=t.bxh.1), color="purple")

ggplot(data=t.b) +
    geom_density(aes(x=t.bl.2), color="red") +
    geom_density(aes(x=t.b.2), color="green") +
    geom_density(aes(x=t.bh.2), color="blue")

ggplot(data=eps.sd) +
    geom_density(aes(x=eps.sdl), color="red") +
    geom_density(aes(x=eps.sd), color="green") +
    geom_density(aes(x=eps.sdh), color="blue")





