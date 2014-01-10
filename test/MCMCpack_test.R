require("MCMCpack")
## Simulated data example with nonconstant choiceset
n <- 1000
y <- matrix(0, n, 4)
colnames(y) <- c("a", "b", "c", "d")
xa <- rnorm(n)
xb <- rnorm(n)
xc <- rnorm(n)
xd <- rnorm(n)
xchoice <- cbind(xa, xb, xc, xd)
z <- rnorm(n)
for (i in 1:n){
    ## randomly determine choiceset (c is always in choiceset)
    choiceset <- c(3, sample(c(1,2,4), 2, replace=FALSE))
    numer <- matrix(0, 4, 1)
    for (j in choiceset){
        if (j == 3){
            numer[j] <- exp(xchoice[i, j] )
        }
        else {
            numer[j] <- exp(xchoice[i, j] - z[i] )
        }
    }
    p <- numer / sum(numer)
    y[i,] <- rmultinom(1, 1, p)
    y[i,-choiceset] <- -999
}
post5 <- MCMCmnl(y~choicevar(xa, "x", "a") +
                     choicevar(xb, "x", "b") +
                     choicevar(xc, "x", "c") +
                     choicevar(xd, "x", "d") + z,
                 baseline="c", verbose=500,
                 mcmc=100000, thin=10, tune=0.01)
plot(post5)


summary(post5)
