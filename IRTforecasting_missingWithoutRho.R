## Model of informative missingness
cat('
model{
    for(i in 1:NM){
    ind[i] ~ dbern(prob[i])
    
    for(j in 1:6){
    ## This is in slope/intercept form, vs IRT form, to match the two-factor model
    ## of missingness:
    x[i,j] <- exp((1-(1-r[j])^2)*(lambda[ifpidx[i],1]*theta[uidx[i],1] - b0[ifpidx[i]]))
    }
    prob[i] <- x[i,fcast[i]]/sum(x[i,1:6])
    }
    
    for(i in 1:N){  ## Users
    theta[i,1:2] ~ dmnorm(zero[1:2], invcor)
    
    ## Model 0/1 vector indicating who responds to what IFP
    for(j in 1:M){
    d[i,j] ~ dbern(pd[i,j])
    probit(pd[i,j]) <- g0[j] + lambda[(M + j),1]*theta[i,1] + lambda[(M + j),2]*theta[i,2]
    }
    }
    
    zero[1] <- 0
    zero[2] <- 0
    invcor <- inverse(cormat)
    cormat[1,1] <- 1
    cormat[2,2] <- 1
    cormat[1,2] <- phi
    cormat[2,1] <- phi
    phi ~ dunif(-1,1)
    
    lambda[1,1] ~ dnorm(0, .04)T(0,)
    lambda[1,2] <- 0
    lambda[(M+1),1] ~ dnorm(0, .04)
    lambda[(M+1),2] ~ dnorm(0, .04)T(0,)
    
    b0[1] ~ dnorm(0, .04)
    g0[1] ~ dnorm(0, .04)
    for(j in 2:M){  ## IFPs
    ## forecast loadings
    lambda[j,1] ~ dnorm(0, 0.04)
    lambda[j,2] <- 0
    ## response loadings
    lambda[(M+j),1] ~ dnorm(0, .04)
    lambda[(M+j),2] ~ dnorm(0, .04)
    
    b0[j] ~ dnorm(0,0.04)
    g0[j] ~ dnorm(0,0.04)
    }
       
    }', file="missmod.jag")

load("fcasts_y2toy4_v3.rda")

dat$uidx <- as.numeric(dat$uidx)
dat$ifpidx <- as.numeric(dat$ifpidx)

## small dataset to make sure the model works
#smalldat <- subset(dat, ifpidx %in% 1:5)

fcast <- ifelse(dat$fcast<0.1, 0,
                ifelse(dat$fcast<0.3,0.2,
                       ifelse(dat$fcast<0.5,0.4,
                              ifelse(dat$fcast<0.7,0.6,
                                     ifelse(dat$fcast<0.9,0.8,1)))))

data <- list(NM = nrow(dat), N = length(unique(dat$uidx)), 
             M = length(unique(dat$ifpidx)), ifpidx = dat$ifpidx, 
             uidx = dat$uidx, fcast=as.numeric(factor(fcast)), r=seq(0,1,length.out=6))
data <- c(data, list(ind=rep(1, data$NM)))

## Define missingness matrix
d <- matrix(0, data$N, data$M)
for(i in 1:nrow(d)){
  for(j in 1:ncol(d)){
    if(any(data$ifpidx==j & data$uidx==i)) d[i,j] <- 1
  }
}

data <- c(data, list(d=d))

inits <- list(b0 = rep(0, data$M), theta = matrix(0, data$N, 2))

## smaller number of iterations, to see how things are working
library(runjags)
mfit2 <- run.jags("missmod.jag", data=data, inits=inits, n.chains=3, monitor=c("theta","b0","lambda","g0"), burnin=5000, sample=20000, method="parallel")
#mfit$summary$statistics[1:10,]

save(mfit2, file="IRTforecasting_missingwithoutRho.rda")

load("IRTforecasting_missingwithoutRho.rda")
gelman.diag(mfit2)
