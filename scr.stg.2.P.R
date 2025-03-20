scr.stg.2.P <- function( n, J, M, X, A, out.1, 
                         beta = c(-2,-.0001), beta.tune = c(.5,.00005),
                         mu.beta = c( -2, -.01 ),
                         s2.beta = c( 1, .1 ) ){
 
###
###  Libraries and Subroutines 
###
  
logit <- function(p){
  log(p)-log(1-p)
}

logit.inv <- function(x){
  exp(x)/(1+exp(x)) 
}

crossdist <- function(S1,S2) { 
    c1 <- complex(real=S1[,1], imaginary=S1[,2]) 
    c2 <- complex(real=S2[,1], imaginary=S2[,2]) 
    dist <- outer(c1, c2, function(z1, z2) Mod(z1-z2)) 
    dist 
} 

###
###  Setup Variables 
###

n.mcmc=out.1$n.mcmc
beta.save=matrix(0,2,n.mcmc)
psi.save=rep(0,n.mcmc)

beta.star=out.1$beta.star
psi.star=out.1$psi.star

###
###  Priors and Starting Values 
###

x.lims=range(A[,1])
y.lims=range(A[,2])


L <- nrow( X )

a=1
b=1

sig.beta=sqrt(s2.beta)
s.tune=1

###
###  Starting values for second stage 
###

n.sim=1000
S.sim=cbind(runif(n.sim,x.lims[1],x.lims[2]),runif(n.sim,y.lims[1],y.lims[2]))
D2.sim=crossdist(S.sim,X)^2  # n.sim x L matrix of distances b/w S.sim and X

idx.star=sample(1:n.mcmc,1)
beta=beta.star[,idx.star]
psi=psi.star[idx.star]

P.sim=logit.inv(beta[1]+beta[2]*D2.sim)  # n.sim x L matrix 
theta.full=psi*(1-apply((1-P.sim)^J,1,prod)) # n.sim x 1, product over L traps 
theta.mat=matrix(theta.full[sample(1:n.sim,M*n.sim,replace=TRUE)],M,n.sim)  # M x n.sim matrix
lambda.vec=apply(theta.mat,2,sum) # sum individ capture probs/intensities 

mh.2=log(mean(dpois(n,lambda.vec)))

###
###  Parallelized MC Integration to get Mh.1
###

library(foreach)
library(doParallel)

cores=detectCores()
print(cores)
cl <- makeCluster(cores[1], setup_strategy = "sequential")
registerDoParallel(cl)

mh.1.parallel=rep(0,n.mcmc)

mh.1.parallel <- foreach(i=1:n.mcmc) %dopar% {
  beta.tmp=beta.star[,i]
  psi.tmp=psi.star[i]

  P.sim.tmp=logit.inv(beta.tmp[1]+beta.tmp[2]*D2.sim)
  theta.tmp.full=psi.tmp*(1-apply((1-P.sim.tmp)^J,1,prod))
  theta.tmp.mat=matrix(theta.tmp.full[sample(1:n.sim,M*n.sim,replace=TRUE)],M,n.sim)
  lambda.tmp.vec=apply(theta.tmp.mat,2,sum)

  log(mean(dpois(n,lambda.tmp.vec)))
}

mh.1.save = unlist(mh.1.parallel)

###
###  Begin Second Stage MCMC Loop 
###

for(k in 1:n.mcmc){
  if(k%%1000==0) cat(k," ")

  ###
  ###  Update beta and psi 
  ###

  idx.star=sample(1:n.mcmc,1)
  mh = exp(mh.1.save[idx.star]-mh.2)
  
  if(mh>runif(1)){
    beta=beta.star[,idx.star]
    psi=psi.star[idx.star]
    mh.2=mh.1.save[idx.star]
  }

  ###
  ###  Save Samples 
  ###
 
  beta.save[,k]=beta
  psi.save[k]=psi

}
cat("\n")

###
###  Write Output 
###

list(psi.save=psi.save,beta.save=beta.save,n.mcmc=n.mcmc,beta.star=beta.star,psi.star=psi.star)

}
