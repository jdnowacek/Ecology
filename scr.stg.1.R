scr.stg.1 <- function( Y, J, M, X, A, n.mcmc, 
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

x.lims = range(A[,1])
y.lims = range(A[,2])

tiny <- 1e-12         # added 

psi = n/M
a=1
b=1

L <- nrow( X )

bmin <- matrix( c( -4, -.01 ), 1 )
bmax <- matrix( c( 0, 0 ), 1 )

sig.beta = sqrt(s2.beta)

n.sim = 1000
S.sim = cbind(runif(n.sim,x.lims[1],x.lims[2]),runif(n.sim,y.lims[1],y.lims[2]))
D2 = t(crossdist(S.sim,X)^2)   # L x n.sim matrix of distances b/w S.sim and X

P = logit.inv(beta[1]+beta[2]*D2)  

p.detect = 1-exp(J*apply(log(1-P),2,sum))
ldenom.int = n*log(mean(p.detect))

mh.2 = sum( apply( Y, 1, function(y,J,P){
  log( mean( exp( apply( dbinom( y, J, P, log = TRUE ), 2, sum) ) ) )
  }, J = J, P = P) ) - ldenom.int + sum( dnorm(beta, mu.beta, sig.beta, log=TRUE) )

###
###  Begin First-Stage MCMC Loop 
###

beta.star = matrix(0,2,n.mcmc)
psi.star  = rbeta(n.mcmc,a,b)

for(k in 1:n.mcmc){
  
  if( k %% 100 == 0 )cat(k," ")
  
  if( is.matrix( beta.tune ) ){
    beta.star.star <- as.vector( rmvnormRcpp( 1, beta, beta.tune ) )
  }else{
    beta.star.star <- rnorm( 2, beta, beta.tune )
  }
  
  if( beta[1] < bmin[1] | beta[2] < bmin[2] ){
    
    beta <- pmax( beta, bmin )
    
    print( beta )
    
    bt <- beta.tune/10
    if( is.matrix( bt ) )bt <- sqrt( diag( bt ) )
    
    beta.star.star <- .tnorm( 2, bmin, bmax, beta, bt )
  }
  

  P.star.star = logit.inv( beta.star.star[1] + beta.star.star[2]*D2 ) # + tiny
  P.star.star[ P.star.star < tiny ] <- tiny
  
  p.detect.star.star = 1 - exp( J*apply( log( 1 - P.star.star ), 2, sum) ) 
  ldenom.int.star.star = n*log( mean( p.detect.star.star ) )

  mh.1 = sum( apply(Y,1,function(y,J,P){ 
    log(mean(exp(apply(dbinom(y,J,P,log=TRUE),2,sum))))},
    J = J,P = P.star.star)) -
    ldenom.int.star.star + sum( dnorm( beta.star.star, mu.beta, sig.beta, log = TRUE) )

  mh = exp( mh.1 - mh.2 )
  if( mh > runif(1) ){
    beta = beta.star.star
    P = P.star.star
    mh.2 = mh.1
  }
  
  beta.star[,k]=beta

};cat("\n")

list(n.mcmc = n.mcmc, beta.star = beta.star, psi.star = psi.star, a = a, b = b)

}

