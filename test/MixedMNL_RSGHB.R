####
# Mixed Multinomial Logit Model
#   Estimation using the RSGHB package
#   follows example 2 of http://cran.r-project.org/web/packages/RSGHB/vignettes/RSGHB_HowTo.pdf
#   underlying algorithm is based on Train's book (http://elsa.berkeley.edu/~train/software.html)
####

Sys.setenv(LANG = "en")
setwd("~/Dropbox/RCode/Choice_Gibbs.git/test")


require("RSGHB")



### True parameters
M <- 3 # number of alternatives (alternative 3 is dummy for no-purchase)
L <- 2 # number of covariates

#XMAT is the attributes of the alternatives; [Xij] is an M*L matrix, i=1...M, j=1...L.
#by col: [X11 X12; X21 X22; 0, 0]
XMAT <- matrix(c(0.4, 0.3, 
                 0.6, 0.1,
                 0, 0), 
               nrow=M, ncol=L, byrow=TRUE);

#true coefficient of beta~mvN(b,W)
b <- c(0.8, 0.2); # L-dimensional
W <- matrix(c(0.04, -0.01, -0.01, 0.01), nrow=L, ncol=L); # L*L matrix


### simulate data
N <- 10000 # number of data points
beta <- t(rmvnorm(N, mean=b, sigma=W)) # L*N matrix
#score of choice 1 and 2 (col) by users (row) is exp(X*beta)
score <- exp(XMAT %*% beta) # M*N matrix
#choice probabilities
choice.prob <- apply(score, 2, function(x) x/sum(x)) # M*N matrix
#simulate actual choice
choice.mat <- apply(choice.prob, 2, function(x) rmultinom(1, 1, x)) # M*N matrix
choice.vec <- t(1:M) %*%  choice.mat   # N-dimensional vector

rowMeans(choice.mat)



### use RSGHB package for estimation
# ------------------
# DATA PREPARATION
# ------------------

choicedata <- data.frame("ID"=seq(1:N))
choicedata$x11 <- rep(XMAT[1,1],N)
choicedata$x12 <- rep(XMAT[1,2],N)
choicedata$x21 <- rep(XMAT[2,1],N)
choicedata$x22 <- rep(XMAT[2,2],N)
choicedata$x31 <- rep(XMAT[3,1],N)
choicedata$x32 <- rep(XMAT[3,2],N)
choicedata$choice <- t(choice.vec)

x11 <- choicedata$x11
x12 <- choicedata$x12
x21 <- choicedata$x21
x22 <- choicedata$x22
x31 <- choicedata$x31
x32 <- choicedata$x32
choice1 <- (choicedata$choice==1)
choice2 <- (choicedata$choice==2)
choice3 <- (choicedata$choice==3)


# ------------------
# ESTIMATION CONTROL
# Setting control list for estimation
# ?doHB for more estimation options
# ------------------

modelname <- "MixedMNL"            # used for output

# Names for the Random Variables
gVarNamesNormal <- c("beta1","beta2")

# For each random variable, specify the distribution for its coefficient
# The options are:
# 1. normal
# 2. log-nomal
# 3. negative log-normal
# 4. normal with all values below zero massed at zero
# 5. Johnson SB with a specified min and max
# gDIST must have an entry for each value in gVarNamesNormal

gDIST <- c(1,1)

# STARTING VALUES
svN <- c(0,0)               # for the random coefficients
# The selection of the mean here is important when working with non-normal distributions

# ITERATION SETTINGS
gNCREP    <- 30000           # Number of iterations to use prior to convergence
gNEREP    <- 20000                   # Number of iterations to keep for averaging after convergence has been reached
gNSKIP    <- 10                 # Number of iterations to do in between retaining draws for averaging
gINFOSKIP <- 100              # How frequently to print info about the iteration process

# CONTROL LIST TO PASS TO doHB
control <- list(
    modelname=modelname,
    gVarNamesNormal=gVarNamesNormal,
    gDIST=gDIST,
    svN=svN,
    gNCREP=gNCREP,
    gNEREP=gNEREP,
    gNSKIP=gNSKIP,gINFOSKIP=gINFOSKIP,
    nodiagnostics=T
)


# ------------------
# likelihood
# USE:     Calculates the likelihood of choice | B
#          Returns likelihood values for each observation
# NOTES:   This is where the bulk of the computation resides so coding this efficiently
#              is essential to reducing run time.
# ------------------
likelihood <- function(fc,b)
{  
    b1  <- b[,1]
    b2  <- b[,2]
    
    v1 <- b1 * x11 + b2 * x12
    v2 <- b1 * x21 + b2 * x22
    v3 <- b1 * x31 + b2 * x32
    
    p  <- (exp(v1)*choice1 + exp(v2)*choice2 + exp(v3)*choice3) / (exp(v1) + exp(v2)+ exp(v3))
    
    return(p)
}



# Estimate the model
doHB(likelihood, choicedata, control)


# review the results
#The _A file contains the sample-level means of the underlying normal at each iteration
plotA("MixedMNL")
plotLog("MixedMNL")

#The _B file contains the average across iterations of the individual level draws for the underlying normals for the random parameters. 
#The _Bsd file provides the standard deviations of those individual draws

