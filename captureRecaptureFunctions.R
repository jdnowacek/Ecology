

library(tictoc)
library(PoissonBinomial)

logit <- function(p){
  log(p/(1-p))
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

fitModel <- function(Y, J, M, X, A, n.mcmc, 
                     mu.beta, beta.tune, s2.beta, n.sim, 
                     mpath = "/Volumes/research/clark/clark.unix/smallMammals/" ){
  
  # mpath is path to hooten files
  
  source( paste( mpath, "scr.stg.1.R", sep = "" ) )  
  source( paste( mpath, "scr.stg.2.P.R", sep = "" ) )  
  
  logit.inv <- function(x){ exp(x)/(1+exp(x)) }
  
  n <- ncol( Y )
  L <- nrow( X )
  
  beta <- rnorm( 2, mu.beta, beta.tune )
  
  cores <- detectCores()
  print( cores )
  cl <- makeCluster( cores[1], setup_strategy = "sequential")
  registerDoParallel(cl)
  
  out.1    <- scr.stg.1( Y, J, M, X, A, n.mcmc, 
                         beta = beta, beta.tune = beta.tune,
                         mu.beta = mu.beta,
                         s2.beta = s2.beta )      ###  Fit SCR Model w/ Stage 1 Algorithm 
  
  accept <- length( which( diff( out.1$beta.star[1,] ) != 0 ) )/n.mcmc
  cat('\nAcceptance:\n' )
  print( accept )
  
  out.2.P <- scr.stg.2.P( n, J, M, X, A, out.1, 
                          beta = beta, beta.tune = beta.tune,
                          mu.beta = mu.beta, s2.beta = s2.beta) ###  Fit SCR Model w/ Stage 2 Algorithm (Poisson)
  
  ###  Sample N from Full-Conditional 
  S.sim  <- cbind( runif( n.sim, min(A[,1]), max(A[,1] ) ),
                   runif( n.sim, min(A[,2]), max(A[,2] ) ) )
  D2.sim <- crossdist( S.sim, X )^2  # n.sim x L matrix of distances b/w S.sim and X
  
  N.2.P.parallel <- rep( 0, n.mcmc)
  N.2.P.parallel <- foreach( kj = 1:n.mcmc ) %dopar% {
    P.2.P.sim    <- logit.inv( out.2.P$beta.save[1,kj] + out.2.P$beta.save[2,kj]*D2.sim)
    num.2.P.tmp  <- out.2.P$psi.save[kj]*apply( (1 - P.2.P.sim)^J, 1, prod)
    psi.2.P.bar  <- mean( num.2.P.tmp/(num.2.P.tmp + 1 - out.2.P$psi.save[kj]) )
    n + rbinom( 1, M - n, psi.2.P.bar)
  }
  N.out.2.P <- unlist(N.2.P.parallel)
  
  bvec <- as.matrix( .chain2tab( t( out.2.P$beta.save ) )[,1:4] )
  bb   <- colnames(bvec)
  bvec <- matrix( t(bvec), 1 )
  colnames( bvec ) <- bnames
  
  nvec <- .chain2tab( unlist(N.2.P.parallel) )
  colnames( nvec ) <- paste( 'N', colnames( nvec ), sep = '.' )
  
  df <- cbind( df, nvec[, 1:4], bvec )
  
  # detection
  P.2.P.mat <- matrix(0, n.mcmc, n.p)
  for(kj in 1:n.p){
    P.2.P.mat[,kj] <- logit.inv( out.2.P$beta.save[1,] + 
                                   out.2.P$beta.save[2,]*(d.vec[kj]^2) )
  }
  
  pvec <- .chain2tab( as.vector( P.2.P.mat ) )  # take over all traps
  colnames( pvec ) <- pnames
  df <- cbind( df, pvec[, 1:4] )
  
  psiVec <- .chain2tab( out.2.P$psi.save )
  colnames( psiVec ) <- qnames
  df <- cbind( df, psiVec[, 1:4] )
  
  ###
  ###  Posterior Power to Detect 
  ###
  
  S.X.sim  <- cbind( runif( n.sim, min(A[,1]), max(A[,1]) ), 
                     runif(n.sim, min(A[,2]), max(A[,2]) ) )
  D2.X.sim <- crossdist( S.sim, X)^2  # n.sim x L matrix of distances b/w S.sim and X
  
  ptd.2.P  <- 0
  for(kj in 1:n.mcmc){
    if(kj %% 10000 == 0){cat(k," ")}
    tmp.idx = sample(1:n.sim,1)
    p.2.P.vec = logit.inv(out.2.P$beta.save[1,kj] + out.2.P$beta.save[2,kj]*D2.sim[tmp.idx,])
    y.2.P.vec = rbinom( L, J, p.2.P.vec)
    ptd.2.P = ptd.2.P + ( sum(y.2.P.vec) > 0)/n.mcmc
  };cat("\n")
  
  df$powerDetect <- ptd.2.P
  
  # half intercept distance
  l0     <- logit.inv( out.2.P$beta.save[1,] )   # intercept (D = 0 )
  lscale <- sqrt( ( logit( l0/2 ) - out.2.P$beta.save[1,] )/out.2.P$beta.save[2,] )
  lvec <- .chain2tab( lscale )
  colnames( lvec ) <- paste( 'L', colnames( lvec ), sep = '.' )
  df <- cbind( df, lvec[,1:4] )
  
  stopCluster(cl)
  
  list( out.1 = out.1, out.2.P = out.2.P, N.out.2.P = N.out.2.P, fit = df )
}


fixKrugerRodentNames <- function( species ){
  
  require( stringr )
  ns <- str_count( species, ' ' )
  species[ ns == 0 ] <- paste( species[ ns == 0 ], 'UNKN' )
  
  species <- .replaceString( species, 'spp', 'UNKN' )
  species <- .replaceString( species, 'Aethemys', 'Aethomys' )
  species <- .replaceString( species, 'Aethmus', 'Aethomys' )
  species <- .replaceString( species, 'Lemniscus', 'Lemniscomys' )
  species <- .replaceString( species, 'Mastemus', 'Mastomys' )
  species
}



fixMammalNames <- function( specNames ){
  
  pcomp <- c('Peromyscus gossypinus/leucopus', 'Peromyscus gossypinus', 'Peromyscus leucopus',
             'Peromyscus leucopus/maniculatus', 'Peromyscus maniculatus', 
             'Peromyscus keeni/maniculatus', 'Peromyscus keeni', 
             'Peromyscus maniculatus/boylii', 'Peromyscus boylii' )
  
  specNames[ specNames == 'Spermophilus richardsonii' ] <- 'Urocitellus richardsonii'
  specNames[ specNames == 'Tamias amoenus' ] <- 'Neotamias amoenus'
  specNames[ specNames == 'Tamias minimus' ] <- 'Neotamias minimus'
  specNames[ specNames == 'Tamias townsendii' ] <- 'Neotamias townsendii'
  specNames[ specNames == 'Tamias speciosus' ] <- 'Neotamias speciosus'
  specNames[ specNames  %in% pcomp ] <- 'Peromyscus boylii-gossypinus-keeni-leucopus-maniculatus'
  specNames[ specNames  %in% c('Dipodomys ordii/microps', 'Dipodomys ordii', 
                               'Dipodomys microps') ] <- 'Dipodomys ordii-microps'
  specNames[ specNames  %in% c('Sorex cinereus/haydeni', 'Sorex cinereus', 'Sorex haydeni') ] <- 
    'Sorex cinereus-haydeni'
  specNames[ specNames == 'Ictidomys tridecemlineatus monticola' ] <- 'Ictidomys tridecemlineatus'
  specNames[ specNames == 'Sigmodon hispidus eremicus' ] <- 'Sigmodon hispidus'
  
  specNames
  
}


getNEONbodySize <- function( y ){
  
  
  y$scientificName <- fixMammalNames( y$scientificName )
  
  
  bout  <- .replaceString( as.character( y$collectDate ), ' GMT', '' )
  dcol  <- columnSplit( bout, '-' )
  dcol  <- matrix( sapply( dcol, as.numeric ), ncol = 3 )
  
  plotSpecYr   <- table( y$plotID, y$scientificName, dcol[, 1] )
  
  gramsMu <- tapply( y$weight, list( y$plotID, y$scientificName, dcol[,1] ),
                     mean, na.rm = T )
  gramsSd <- tapply( y$weight, list( y$plotID, y$scientificName, dcol[,1] ),
                     sd, na.rm = T )
  gramsN  <- tapply( y$weight*0 + 1, list( y$plotID, y$scientificName, dcol[,1] ),
                     sum, na.rm = T )
  
  specMu <- apply( gramsMu, 2, mean, na.rm = T )
  specSd <- sqrt( apply( gramsSd^2, 2, mean, na.rm = T ) )
  
  specGm <- cbind( specMu, specSd )
  colnames( specGm ) <- c('gMu','gSd')
  
  list( gramsMu = gramsMu, gramsSd = gramsSd, gramsN = gramsN, 
        specGm = specGm )
}


mapCapture <- function( n, A, X, Y){
  
  pl <- .getPlotLayout( n, WIDE = TRUE )$mfrow
  
  par( mar = c(1,1,2,1) )
  
  layout( matrix( 1:prod( pl ), pl[1], pl[2], byrow = TRUE) )

  
 # layout( matrix(1:15, 3, 5, byrow = TRUE) )
  
  for(i in 1:n){
    plot( A, type="n", asp = TRUE, main=i, xlab = "", ylab = "",
          xaxt="n", yaxt="n", bty = "n" )
    text( matrix( X[ Y[, i] == 0, ], ncol = 2 ), labels = Y[ Y[, i]==0, i ],
          cex = 1, col = rgb(0,0,0,.3) )
    text( matrix( X[ Y[ ,i ] > 0,], ncol = 2 ), labels = Y[ Y[, i] > 0, i ],
          cex = 1.5, col=1)
    polygon( A, lwd = 1.5, lty = 3, border = rgb(0,0,0,.5) )
  }
}


cor2Cov <- function(sigvec,cormat){ 
  
  #correlation matrix and variance vector to covariance
  
  n <- nrow(cormat)
  d <- length(sigvec)
  if(d == 1)sigvec <- rep( sigvec, n )
  
  s <- matrix(sigvec,n,n)
  cormat*sqrt(s*t(s))
}

kernelPlot <- function( P.2.P.mat, out.2.P  ){
  
  
  par( bty = 'n' )
  
  col.1.1=rgb(1,.65,0,alpha=.8)
  col.1.2=rgb(1,.65,0,alpha=.4)
  col.2.1=rgb(.4,.6,.8,alpha=.8)
  col.2.2=rgb(.4,.6,.8,alpha=.4)
  
  n.p = 100
  d.vec = seq(0, max(dist(X))/1.5, length = n.p)
  P.2.P.mat=matrix(0,n.mcmc,n.p)
  for(j in 1:n.p){
    lp <- out.2.P$beta.save[1,] + out.2.P$beta.save[2,]*(d.vec[j]^2)
    P.2.P.mat[,j] = logit.inv( lp )
  }
  P.2.P.l = apply(P.2.P.mat,2,quantile,.025)
  P.2.P.u = apply(P.2.P.mat,2,quantile,.975)
  P.2.P.m = apply(P.2.P.mat,2,mean)
  loHi <- cbind( P.2.P.l, P.2.P.u )
  
  plot( d.vec, P.2.P.m, ylim = c( 0, max( loHi ) ), type="n", xlab="Distance (m)", ylab="p")
  
  .shadeInterval( d.vec, loHi, add = T )
  lines(d.vec,P.2.P.m,col=col.1.1,lwd=2)

}

posteriorPlots <- function( out.2.Pb, N.out.2.Pb, n, M ){
  
  
  col.1.1=rgb(1,.65,0,alpha=.8)
  col.1.2=rgb(1,.65,0,alpha=.4)
  col.2.1=rgb(.4,.6,.8,alpha=.8)
  col.2.2=rgb(.4,.6,.8,alpha=.4)
  
  bw1 = bw.nrd( out.2.Pb$beta.save[1,] )*1.7
  bw2 = bw.nrd( out.2.Pb$beta.save[2,] )*1.7
  bw3 = bw.nrd( out.2.Pb$psi.save )*1.7
  
  tab.1 = table(c(N.out.2.Pb,seq(n,M,1)))/n.mcmc
  attributes(tab.1)$dimnames[[1]] = as.numeric(attributes(tab.1)$dimnames[[1]]) + .25
  
  par( bty = 'n', mar = c(4,4,1,1 ) )
  
  layout(matrix(1:4,2,2))
  
  plot( density( out.2.Pb$beta.save[1,], bw = bw1 ), lwd=2, col=col.1.1, main = '',
        ylab="density", xlab = bquote(beta[0]) )
  abline( v = quantile( out.2.Pb$beta.save[1,], c(.025, .5, .975 ) ), lty = 2 )
  
  plot( density(out.2.Pb$beta.save[2,], bw=bw2 ), lwd=2,col=col.1.1, main = '',
        ylab="density", xlab = bquote(beta[1]) )
  abline( v = quantile( out.2.Pb$beta.save[2,], c(.025, .5, .975 ) ), lty = 2 )
  
  plot( density(out.2.Pb$psi.save, bw=bw3 ), lwd =2 , col = col.1.1, main = '',
        ylab="density", xlab = bquote(psi), xlim = c(0,1) )
  abline( v = quantile( out.2.Pb$psi.save, c(.025, .5, .975 ) ), lty = 2 )
  
  plot( tab.1, col = col.1.1, type="h", xlim = c(0, M ),
       ylab = "probability", xlab = "N", lwd=2, xaxt="n", main = '' )
  
  axis( 1, round( seq( 0, M, 40 ) ) )
#  abline( v=n, lwd=1.5, lty=3, col=rgb(0,0,0,.7))
  
  yy <- par('usr')[4]*.2
  text( M, yy, 'M' )
  text( n, yy, 'n' )
}


.getPlotLayout <- function( np, WIDE = TRUE ){
  
  # np - no. plots
  
  if( np == 1 )return( list( mfrow = c( 1, 1 ), left = 1, bottom = c( 1, 2 ) ) )
  if( np == 2 ){
    if( WIDE )return( list( mfrow = c( 1, 2 ), left = 1, bottom = c( 1, 2 ) ) )
    return( list( mfrow = c( 2, 1 ), left = c( 1, 2 ), bottom = 2 ) )
  }
  
  if( np == 3 ){
    if( WIDE )return( list( mfrow = c( 1, 3 ), left = 1, bottom = c( 1:3 ) ) )
    return( list( mfrow = c( 3, 1 ), left = 1:3, bottom = 3 ) )
  }
  if( np <= 4 )return( list( mfrow = c( 2, 2 ), left = c( 1, 3 ), bottom = c( 3:4 ) ) )
  if( np <= 6 ){
    if( WIDE )return( list( mfrow = c( 2, 3 ), left = c( 1, 4 ), bottom = c( 4:6 ) ) )
    return( list( mfrow = c( 3, 2 ), left = c( 1, 3, 5 ), bottom = 5:6 ) )
  }
  if( np <= 9 )return( list( mfrow = c( 3, 3 ), left = c( 1, 4, 7 ), bottom = c( 7:9 ) ) )
  if( np <= 12 ){
    if( WIDE )return( list( mfrow = c( 3, 4 ), left = c( 1, 5, 9 ), bottom = c( 9:12 ) ) )
    return( list( mfrow = c( 4, 3 ), left = c( 1, 4, 7, 10 ), bottom = 10:12 ) )
  }
  if( np <= 16 )return( list( mfrow = c( 4, 4 ), left = c( 1, 5, 9, 13 ), 
                              bottom = c( 13:16 ) ) )
  if( np <= 20 ){
    if( WIDE )return( list( mfrow = c( 4, 5 ), left = c( 1, 6, 11, 15 ), 
                            bottom = c( 15:20 ) ) )
    return( list( mfrow = c( 5, 4 ), left = c( 1, 5, 9, 13 ), bottom = 17:20 ) )
  }
  if( np <= 25 )return( list( mfrow = c( 5, 5 ), left = c( 1, 6, 11, 16, 21 ), 
                              bottom = c( 20:25 ) ) )
  if( np <= 30 ){
    if( WIDE )return( list( mfrow = c( 5, 6 ), left = c( 1, 7, 13, 19, 25 ), 
                            bottom = c( 25:30 ) ) )
    return( list( mfrow = c( 6, 5 ), left = c( 1, 6, 11, 16, 21, 26 ), bottom = 26:30 ) )
  }
  if( np <= 36 ){
    return( list( mfrow = c( 6, 6 ), left = c( 1, 7, 13, 19, 25, 31 ), bottom = c( 31:36 ) ) )
  }
  return( list( mfrow = c( 7, 6 ), left = c( 1, 7, 13, 19, 25, 31, 37 ), bottom = c( 37:42 ) ) )
}


getNEONtee <- function( gramsMu, gramsSd, specGm ){
  
  # gramsMu, gramsSd are plot X species X year arrays of body size in g
  # returns TEE by plotYr and diet type
  
  plotData <- read.csv( 'plotData.csv', row.names = 1 )
  Jmat <- as.matrix (plotData[ , startsWith( colnames( plotData), 'X' ) ] ) #bouts
  colnames( Jmat ) <- .replaceString( colnames( Jmat ), 'X', '' )
  
  ocols    <- c('scientificName', 'commonName')
  dietCols <- c( 'det_invert', 'det_mamsBirds', 'det_reptAmph', 'det_fish',
                 'det_vertUnkn', 'det_scav', 'det_fruit', 'det_nect', 
                 'det_seed', 'det_plantother' )
  
  nspace <- stringr::str_count( rownames( specGm ), ' ' )
  w2     <- which( nspace > 1 )
  
  if( length( w2 ) > 0 ){
    print( specGm[w2,] )
    stop( 'check species name' )
  }
  
  genus <- columnSplit( rownames( specGm ), ' ' )[,1]
  genera <- unique( genus )
  family <- traits$family[ match( genus, traits$genus ) ]
  
  gwt <- tapply( traits$adult_mass_g, traits$genus, mean, na.rm = T )[genera]
  fwt <- tapply( traits$adult_mass_g, traits$family, mean, na.rm = T )[genera]
  
  dm  <- as.matrix( traits[, dietCols] )
  gdi <- tapply( as.vector(dm), list( rep( traits$genus, length( dietCols ) ),
                                           rep( dietCols, each = nrow(traits) ) ), 
                                           mean, na.rm = T )
  fdi <- tapply( as.vector(dm), list( rep( traits$family, length( dietCols ) ),
                                                        rep( dietCols, each = nrow(traits) ) ), 
                 mean, na.rm = T )
  
  traits <- read.csv( '../traitsByGroup/mammalTraits.csv' )
  diet   <- which( startsWith( colnames(traits), 'det_' ) )
  wna <- which( is.na( specGm[,1] ) )
  
  snew <- specGm[,1]
  mm     <- match( rownames( specGm ), traits$scientificName )
  
  diet <- as.matrix( traits[mm,  dietCols ] )
  diet <- diet[ , colSums( diet, na.rm = T ) > 0 ]
  rownames(diet) <- rownames( specGm )
  
  # genus
  wna <- unique( which( is.na( diet ), arr.ind = T )[,1] )
  mm  <- match( genus[ wna ], rownames(gdi) )
  gtab <- gdi[mm,]
  diet[ wna, ] <- gtab
  
  # family
  wna <- unique( which( is.na( diet ), arr.ind = T )[,1] )
  mm  <- match( family[ wna ], rownames(fdi) )
  gtab <- gdi[mm,]
  diet[ wna, ] <- gtab
  
  diet <-  sweep( diet, 1, rowSums( diet, na.rm = T ), '/')
  
  lf <- list.files( 'output/', recursive = T )
  
  specPlot <- columnSplit( lf, '/' )
  specPlot[,2] <- .replaceString(specPlot[,2], '.csv', '' )
  specPlot[,1] <- .replaceString(specPlot[,1], '_', ' ' )
  
  lf <- list.files( 'output/', recursive = T, full.names = T )
  
  Nmu <- NVr <- NN <- numeric(0)
  Wmu <- WVr <- numeric(0)
  
  teeByGroupMu <- teeByGroupVr <- teeByGroupN <- numeric(0)
  
  for( j in 1:length(lf) ){
    
    tmp <- read.csv( lf[j] )
    J <- Jmat[ specPlot[j,2], ]
    yr <- as.numeric( names( J ) )
    yr <- min( yr[ J > 0 ] ):max( yr[ J > 0 ] )
    J  <- J[ as.character(yr) ]
    
    if( !specPlot[j,1] %in% rownames( diet ) ){
      cat('\nspecies missing from diet matrix:\n')
      print( specPlot[j,1] )
      next
    }
    
    dfrac <- diet[ drop = F, specPlot[j,1], ]
    
    mmat <- matrix( NA, 1, length( yr ), 
                    dimnames = list( specPlot[j,2], yr ) )
    mmat[ 1, J > 0 ] <- 0
    mmat[ 1, as.character( tmp$year ) ] <- tmp$N.estimate 
    smat <- mmat
    smat[ 1, as.character( tmp$year ) ] <- tmp$N.SE 
    
    
    wmat <- gramsMu[ specPlot[j,2], specPlot[j,1], ]
    wmat[ !is.finite( wmat ) ] <- specGm[ specPlot[j, 1], 1 ]
    wmat <- wmat[ as.character( yr ) ]
    tMu  <- kg2TEE( wmat[ colnames(mmat) ]/1000 )*mmat
    tSd  <- sqrt( (kg2TEE( wmat[ colnames(mmat) ]/1000 )*smat)^2 )
    
    tfmu <- t(dfrac)%*%tMu
    tsmu <- t(dfrac)%*%tSd 
    wk   <- which( rowSums( tfmu, na.rm = T ) != 0 )
    
    if( length( wk ) == 0 ) next
    
    tfmu <- tfmu[ drop = F, wk, ]
    tsmu <- tsmu[ drop = F,wk, ]
    rownames( tfmu ) <- rownames( tsmu ) <- paste( specPlot[j, 2], rownames( tfmu ), sep = '-' )
    
    teeByGroupMu <- .appendMatrix( teeByGroupMu, tfmu )   # these are summed
    teeByGroupVr <- .appendMatrix( teeByGroupVr, tsmu^2 )
    nn           <- tfmu*0 + 1
    nn[ is.na(nn) ] <- 0
    
    teeByGroupN  <- .appendMatrix( teeByGroupN, nn ) # species that have contributed
    
    Nmu <- .appendMatrix( Nmu, tMu )                # sum for plot
    NVr <- .appendMatrix( NVr, tSd^2 )
    nn  <- mmat*0 + 1
    nn[ is.na(nn) ] <- 0
    NN  <- .appendMatrix( NN, mmat*0 + 1 )           # species that have contributed
  }
  
  list( TEEmu = Nmu, TEEsd = sqrt( NVr ), nspec = NN, 
        teeByGroupMu = teeByGroupMu, teeByGroupSd = sqrt( teeByGroupVr ), 
        teeByGroupNspec = teeByGroupN )
}





.tnormMVNmatrix <- function(avec, muvec, smat, 
                            lo=matrix(-1000,nrow(muvec),ncol(muvec)), 
                            hi=matrix(1000,nrow(muvec),ncol(muvec)),
                            whichSample = c(1:nrow(smat)) ){
  
  #lo, hi must be same dimensions as muvec, avec
  
  lo[lo < -1000] <- -1000
  hi[hi > 1000]  <- 1000
  
  if(max(whichSample) > length(muvec))
    stop('whichSample outside length(muvec)')
  
  r <- avec
  a <- trMVNmatrixRcpp(avec, muvec, smat, lo, hi, whichSample, 
                       idxALL = c(0:(nrow(smat)-1)) )  
  r[,whichSample] <- a[,whichSample]
  r
}

.tnorm <- function( n, lo, hi, mu, sig, tiny = 0 ){   
  
  #normal truncated lo and hi
  
  if( length( lo ) == 1 & length( mu ) > 1 )lo <- rep( lo, length( mu ) )
  if( length( hi ) == 1 & length( mu ) > 1 )hi <- rep( hi, length( mu ) )
  
  if( length( lo ) != length( mu ) ){
    print( length( lo ) )
    print( length( mu ) )
    stop( )
  }
  
  q1 <- pnorm( lo, mu, sig )
  q2 <- pnorm( hi, mu, sig ) 
  
  z <- runif( n, q1, q2 )
  z <- qnorm( z, mu, sig )
  
  z[ z == Inf]  <- lo[ z == Inf] + tiny
  z[ z == -Inf] <- hi[ z == -Inf] + tiny
  z[ z > hi] <- hi[ z > hi]
  z[ z < lo] <- lo[ z < lo]
  
  z
}



.replaceString <- function( xx, now = '_', new = ' ' ){  #replace now string in vector with new
  
  if( !is.character( xx ) )xx <- as.character( xx )
  
  ww <- grep( now[ 1], xx, fixed = TRUE )
  if( length( ww ) == 0 )return( xx )
  
  if( length( new ) == 1 ){
    for( k in ww ){
      s  <- unlist( strsplit( xx[ k], now, fixed = TRUE ) )
      ss <- s[ 1]
      if( length( s ) == 1 )ss <- paste( ss, new, sep = '' )
      if( length( s ) > 1 )for( kk in 2:length( s ) ) ss <- paste( ss, s[ kk], sep = new )
      xx[ k] <- ss
    }
  }else{
    # new is a vector
    s  <- unlist( strsplit( xx, now, fixed = TRUE ) )
    nx <- length( xx )
    nc <- length( s )/length( xx )
    
    ss <- matrix( s, ncol = nc, byrow = TRUE )
    nl <- nchar( ss[ 1, ] )
    
    if( nl[ 1] == 0 )ss <- paste( new, ss[, 2], sep = '' )
    if( nl[ 2] == 0 )ss <- paste( ss[, 1], new, sep = '' )
    
    xx <- ss
  }
  xx
}


columnSplit <- function( vec, sep = '_', ASFACTOR = F, ASNUMERIC = FALSE, 
                         LASTONLY = FALSE ){
  
  vec <- as.character( vec )
  nc  <- length( strsplit( vec[ 1], sep, fixed = TRUE )[[ 1]] )
  
  mat <- matrix( unlist( strsplit( vec, sep, fixed = TRUE ) ), ncol = nc, byrow = TRUE )
  if( LASTONLY & ncol( mat ) > 2 ){
    rnn <- mat[, 1]
    for( k in 2:( ncol( mat )-1 ) ){
      rnn <- columnPaste( rnn, mat[, k] )
    }
    mat <- cbind( rnn, mat[, ncol( mat )] )
  }
  if( ASNUMERIC ){
    mat <- matrix( as.numeric( mat ), ncol = nc )
  }
  if( ASFACTOR ){
    mat <- data.frame( mat )
  }
  if( LASTONLY )mat <- mat[, 2]
  mat
}


.chain2tab <- function( chain, sigfigs = 3 ){
  
  if( !is.matrix( chain ) )chain <- matrix( chain, ncol = 1 )
  
  mu <- colMeans( chain )    
  
  SE <- apply( chain, 2, sd, na.rm=T )
  CI <- apply( chain, 2, quantile, c( .025, .975 ), na.rm=T )
  splus <- rep( ' ', length = length( SE ) )
  splus[ CI[ 1, ] > 0 | CI[ 2, ] < 0] <- '*'
  
  tab <- cbind( mu, SE, t( CI ) )
  tab <- signif( tab, sigfigs )
  colnames( tab ) <- c( 'estimate', 'SE', 'CI_025', 'CI_975' )
  tab <- as.data.frame( tab )
  tab$sig95 <- splus
  attr( tab, 'note' ) <- '* indicates that zero is outside the 95% CI'
  
  point <- grep( '.', rownames( tab ), fixed = T )  # re-insert dashes in species names
  if( length( point ) > 0 ){
    under <- grep( '_', rownames( tab )[point] )
    if( length( under ) > 0 ){
      ss <- columnSplit( rownames(tab)[ point[under] ], '_' )
      sd <- .replaceString( ss[,1], '.', '-' )
      rownames( tab )[ point[under] ] <- columnPaste( sd, ss[,2], '_' )
    }else{
      rownames( tab )[point] <- .replaceString( rownames( tab )[point], '.', '-' )
    }
  }
  xint <- grep( ':', colnames(chain), fixed = T )
  if( length( xint ) > 0 ){
    ss <- columnSplit( rownames(tab)[xint], '_' )
    sd <- .replaceString( ss[,2], '.', ':' )
    rownames( tab )[xint] <- columnPaste( ss[,1], sd, '_' )
  }
  
  tab
}
