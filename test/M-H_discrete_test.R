# implementation of discrete Metropolis-Hastings algorithm with Negative Binomial proposal
# tested for one-dimensional discrete distribution
discreteMH <-function (logpost, proposal, start, m, ...) 
{
    pb = length(start)
    Mpar = array(0, c(m, pb))
    b = matrix(t(start))
    lb = logpost(start, ...)

    size = proposal$size
    
    accept = 0
    for (i in 1:m) {
        #proposed next point is Negative Binomial using the current position as its mean, variance is controled by $size$
        bc <- rnbinom(pb, size=size, mu=b)
        
        lbc = logpost(t(bc), ...) 
        #since proposal is asymmetric, need to adjust the jumping probability accordingly        
        prob = exp(lbc + dnbinom(b, size=size, mu=bc, log=TRUE) - lb - dnbinom(bc, size=size, mu=b, log=TRUE) )
        
        if (is.na(prob) == FALSE) {
            if (runif(1) < prob) {
                lb = lbc
                b = bc
                accept = accept + 1
            }
        }
        Mpar[i, ] = b
    }
    accept = accept/m
    stuff = list(par = Mpar, accept = accept)
    return(stuff)
}



logpost <- function(x, lambda){
    dpois(x, lambda, log=TRUE)
}


#sample
proposal <- list(size=1)
start <- 1
m <- 10000

s <- discreteMH(logpost, proposal,start,m, lambda=15)



s$accept
plot(s$par, type="l")

s.truncated <- s$par[1000:10000]
mean(s.truncated)
hist(s.truncated, 50, freq=FALSE, main="histogram")

