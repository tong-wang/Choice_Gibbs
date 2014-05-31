####
# Joint Price & Inventory optimization
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
posteriors$posterior.id <- 1:nrow(posteriors)
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



#remove all objects except for "posterior"
rm(list=ls()[ls()!="posteriors"])



#price candidates
Cost <- 1.4
price.range <- seq(1.5, 1.7, 0.01) * Cost
inv.range <- 15:25



## generate profit samples
num <- 10

generate.profit <- function(row, scenario, num) {
    
    out1 <- expand.grid(price.range, 1:num)
    colnames(out1) <- c("price", "sample.id")
    
    elambda <- row[, paste0("lam", scenario)]  / (1 + exp(- out1$price * row[, paste0("beta",scenario, ".1")] - row[, paste0("beta",scenario, ".2")]))
    out1$demand <- rpois(nrow(out1), elambda)
    
    out <- data.frame()
    for (i in inv.range) {
        out <- rbind(out, cbind(out1, data.frame(inv=rep(i, nrow(out1)))))
    }
    
    out$profit <- out$price * pmin(out$demand, out$inv) - Cost * out$inv

    return(out)
}

profit0 <- ddply(posteriors, .(posterior.id), generate.profit, scenario="0", num=num)
profit1 <- ddply(posteriors, .(posterior.id), generate.profit, scenario="1", num=num)
profit3.m <- ddply(posteriors, .(posterior.id), generate.profit, scenario="3.m", num=num)
profit0c <- ddply(posteriors, .(posterior.id), generate.profit, scenario="0c", num=num)
profit1c <- ddply(posteriors, .(posterior.id), generate.profit, scenario="1c", num=num)
profit3c.m <- ddply(posteriors, .(posterior.id), generate.profit, scenario="3c.m", num=num)

head(profit0)


####### add estimates from EM methods ##########

load("MNL.M1.L1_Scenario0.EM.RData")
#load("MNL.M1.L1_Scenario1c.EM.RData")
load("MNL.M1.L1_Scenario3Mc.m.EM.RData")

posteriors.EM <- data.frame(
    posterior.id = 1,
    lam0 = tail(z0.EM$lambdas, n=1),
    lam3c.m = tail(z3Mc.m.EM$lambdas, n=1),
    beta0 = tail(z0.EM$betas, n=1),
    beta3c.m = tail(z3Mc.m.EM$betas, n=1)
)

num.EM <- 10000

profit0.EM <- ddply(posteriors.EM, .(posterior.id), generate.profit, scenario="0", num=num.EM)
profit3c.m.EM <- ddply(posteriors.EM, .(posterior.id), generate.profit, scenario="3c.m", num=num.EM)

################################################




# total mean, var, quantile
meanvar0 <- ddply(profit0, .(price, inv), summarize, mean=mean(profit), var=var(profit))
meanvar1 <- ddply(profit1, .(price, inv), summarize, mean=mean(profit), var=var(profit))
meanvar3.m <- ddply(profit3.m, .(price, inv), summarize, mean=mean(profit), var=var(profit))
meanvar0c <- ddply(profit0c, .(price, inv), summarize, mean=mean(profit), var=var(profit))
meanvar1c <- ddply(profit1c, .(price, inv), summarize, mean=mean(profit), var=var(profit))
meanvar3c.m <- ddply(profit3c.m, .(price, inv), summarize, mean=mean(profit), var=var(profit))
meanvar0.EM <- ddply(profit0.EM, .(price, inv), summarize, mean=mean(profit), var=var(profit))
meanvar3c.m.EM <- ddply(profit3c.m.EM, .(price, inv), summarize, mean=mean(profit), var=var(profit))


#plot efficient frontier of total mean-variance
ggplot() +
    geom_point(data=meanvar0, aes(x=mean, y=var)) + 
    #geom_point(data=meanvar1, aes(x=mean, y=var), color="grey") + 
    #geom_point(data=meanvar3.m, aes(x=mean, y=var), color="red") + 
    #geom_point(data=meanvar0c, aes(x=mean, y=var), color="yellow") + 
    geom_point(data=meanvar1c, aes(x=mean, y=var), color="darkgrey") + 
    geom_point(data=meanvar3c.m, aes(x=mean, y=var), color="blue")  +
    geom_point(data=meanvar3c.m.EM, aes(x=mean, y=var), color="green")  


#optimal price-inventory (risk neutral)
meanvar0[which.max(meanvar0$mean),]
meanvar1[which.max(meanvar1$mean),]
meanvar3.m[which.max(meanvar3.m$mean),]
meanvar0c[which.max(meanvar0c$mean),]
meanvar1c[which.max(meanvar1c$mean),]
meanvar3c.m[which.max(meanvar3c.m$mean),]
meanvar3c.m.EM[which.max(meanvar3c.m.EM$mean),]


#optimal price-inventory (risk averse: mean >= b*var)
b <- 1.5
meanvar0.ra <- meanvar0[meanvar0$mean >= b*meanvar0$var,]
meanvar1c.ra <- meanvar1c[meanvar1c$mean >= b*meanvar1c$var,]
meanvar3c.m.ra <- meanvar3c.m[meanvar3c.m$mean >= b*meanvar3c.m$var,]
meanvar3c.m.EM.ra <- meanvar3c.m.EM[meanvar3c.m.EM$mean >= b*meanvar3c.m.EM$var,]

meanvar0.ra[which.max(meanvar0.ra$mean),] 
meanvar1c.ra[which.max(meanvar1c.ra$mean),] 
meanvar3c.m.ra[which.max(meanvar3c.m.ra$mean),] 
meanvar3c.m.EM.ra[which.max(meanvar3c.m.EM.ra$mean),] 



# conditional mean and variance (conditioning on elambda)
cmeanvar0 <- ddply(profit0, .(posterior.id, price, inv), summarize, mean=mean(profit), var=var(profit))
#cmeanvar1 <- ddply(profit1, .(posterior.id, price, inv), summarize, mean=mean(profit), var=var(profit))
#cmeanvar3.m <- ddply(profit3.m, .(posterior.id, price, inv), summarize, mean=mean(profit), var=var(profit))
#cmeanvar0c <- ddply(profit0c, .(posterior.id, price, inv), summarize, mean=mean(profit), var=var(profit))
cmeanvar1c <- ddply(profit1c, .(posterior.id, price, inv), summarize, mean=mean(profit), var=var(profit))
cmeanvar3c.m <- ddply(profit3c.m, .(posterior.id, price, inv), summarize, mean=mean(profit), var=var(profit))
#cmeanvar0.EM <- ddply(profit0.EM, .(posterior.id, price, inv), summarize, mean=mean(profit), var=var(profit))
cmeanvar3c.m.EM <- ddply(profit3c.m.EM, .(posterior.id, price, inv), summarize, mean=mean(profit), var=var(profit))


# variace of the conditional mean
var.cmean0 <- ddply(cmeanvar0, .(price, inv), summarize, var=var(mean))
#var.cmean1 <- ddply(cmeanvar1, .(price, inv), summarize, var=var(mean))
#var.cmean3.m <- ddply(cmeanvar3.m, .(price, inv), summarize, var=var(mean))
#var.cmean0c <- ddply(cmeanvar0c, .(price, inv), summarize, var=var(mean))
var.cmean1c <- ddply(cmeanvar1c, .(price, inv), summarize, var=var(mean))
var.cmean3c.m <- ddply(cmeanvar3c.m, .(price, inv), summarize, var=var(mean))
var.cmean3c.m.EM <- ddply(cmeanvar3c.m.EM, .(price, inv), summarize, var=var(mean))






