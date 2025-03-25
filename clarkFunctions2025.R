colRed2Blu <- c('#67001f','#b2182b','#d6604d','#f4a582','#fddbc7','white',
                '#d1e5f0','#92c5de','#4393c3','#2166ac','#053061')
colContrast <- c('#1f78b4','#33a02c','#e31a1c',
                 '#ff7f00','#6a3d9a','#b15928')


colorScaleLegend <- function( xleg, yleg = NULL, zlim, 
                              zlimLabs = zlim, box = T,
                              xlim = NULL, ylim = NULL, ncol = 20, units = '',
                              colorRamp, textColor = NULL, cex = 1.2 ){
  
  xchar <- max( strwidth( zlimLabs ) )
  
  xusr <- par('usr')[1:2]
  yusr <- par('usr')[3:4]
  
  if( is.null( xlim ) )xlim <- xusr
  if( is.null( ylim ) )ylim <- yusr
  
  dx <- diff( xlim )
  dy <- diff( ylim )
  
  tpos <- 2
  xtext <- xlim[2] - dx/150
  
  if( is.character( xleg ) ){
    pos <- c('bottomleft', 'bottomright','topleft','topright')
    if( !xleg %in% pos )
      stop( paste( 'legend position must be:', paste0( pos, collapse = ', ' ) ) )
    
    xwide <- dx/25
    yhigh <- dy/2
    
    pos <- xleg
    
    xleg  <- xlim[1] + c(-2*xwide, -xwide ) # bottomleft
    yleg  <- ylim[1] + c(0, yhigh )
    tpos  <- 4                              # right of scale
    xtext <- xleg[2] + dx/150
    
    if( pos == 'topleft' )yleg <- ylim[2] + c(-yhigh, -yhigh/5 )
    
    if( xusr[1] < xlim[1] ){
      xleg <- xleg - xwide
      xtext <- xtext - xwide
    }
    
    if( pos == 'bottomright' ){
      xleg <- xlim[2] + c(0, xwide )
      yleg <- ylim[1] + c(0, yhigh )
      tpos <- 2                     # left of scale
      xtext <- xlim[2] - dx/150
      if( xusr[2] > xlim[2] ){
        xleg <- xleg + xwide
        xtext <- xtext + xwide
      }
    }
    if( pos == 'topright' ){
      xleg <- xlim[2] + c(0, xwide )
      yleg <- ylim[2] + c( -yhigh, 0 )
      tpos <- 2                     # left of scale
      xtext <- xleg[1] - dx/150
      if( xusr[2] > xlim[2] ){
        xleg <- xleg + xwide
        xtext <- xtext + xwide
      }
    }
  }else{
    lpos <- 'left'
    xtext <- xleg[2] + xtext
    if( xleg[1] > (xlim[1] + dx/2) ){
      lpos <- 'right'
      xtext <- xleg[1] - dx/150 
    }
  }
  
  xbox <- xleg
  ybox <- yleg + dy/20*c(-1, 1)
  if( nchar( units ) > 0 )ybox[2] <- ybox[2] + dy/20
  if( tpos == 2 )xbox[1] <- xbox[1] - xchar - dx/10
  if( tpos == 4 )xbox[2] <- xbox[2] + xchar + dx/10
  
  lcol <- colorRampPalette(colorRamp)( ncol )
  zseq <- seq( yleg[1], yleg[2], length = ncol + 1 )
  
  par( xpd = T )
  
  if( box )rect( xbox[1], ybox[1], xbox[2], ybox[2], col = 'white', border = NA )
  for( i in 1:ncol ){
    yi <- zseq[i:(i+1)]
    rect(xleg[1], yi[1], xleg[2], yi[2], col = lcol[i], border = NA )
  }
  
  if( is.null( textColor ) )textColor <- colorRampPalette(colorRamp)( 2 )
  text( rep( xtext, 2 ), yleg, zlimLabs, pos = tpos, col = textColor, 
        cex = cex )
  
  if( nchar( units ) > 0 )text( mean( xleg ), yleg[2] + diff(yleg)/3, units,
                                cex = cex )
  par( xpd = F )
}
getMunicipality <- function(lonLat, collapse = T){
  
  require(revgeo)
  require(countrycode)
  
  ll <- columnPaste(lonLat[,1],lonLat[,2])
  ww <- which(!duplicated(ll))
  
  lc <- lonLat[ww,]
  
  cc  <- vector('character', length = length(ww))
  geo <- data.frame( city = cc, state = cc, country = cc,
                     stringsAsFactors = F)
  
  z <- 25
  m <- 1:z
  
  for(k in 1:length(ww)){
    
    mm  <- m[ m < length(ww) ]
    gk <- revgeo(longitude = lc[mm,1], latitude = lc[mm,2], output='frame')
    geo1 <- cleanMunicipality( gk$city )
    geo2 <- cleanMunicipality( gk$state )
    geo2 <- abb2state(geo2)
    geo3 <- as.character( gk$country )
    geo3 <- countrycode(geo3, 'country.name', 'iso3c')
    
    geo[mm,] <- cbind(geo1, geo2, geo3)
    
    if(max(m) >= length(ww))break
    m <- m + z
    
    print('pause')
    Sys.sleep(5)
    
  }
  
  mm <- match(ll, columnPaste( lc[,1], lc[,2]) )
  colnames(lonLat) <- c('lon','lat')
  
  detach( package:countrycode )
  
  data.frame( lonLat, geo[mm,], stringsAsFactors = F)
}

.shadeInterval <- function( xvalues, loHi, col = 'grey', PLOT = TRUE, add = TRUE, 
          xlab = ' ', ylab = ' ', xlim = NULL, ylim = NULL, 
          LOG = FALSE, trans = .5 ){
  tmp <- NULL
  
  #draw shaded interval
  
  loHi <- as.matrix( loHi )
  tmp  <- smooth.na( xvalues, loHi )
  
  xvalues <- tmp[, 1]
  loHi    <- tmp[, -1]
  
  xbound <- c( xvalues, rev( xvalues ) )
  ybound <- c( loHi[, 1], rev( loHi[, 2] ) )
  if( is.null( ylim ) )ylim <- range( as.numeric( loHi ) )
  if( is.null( xlim ) )xlim <- range( xvalues )
  
  if( !add ){
    if( !LOG )plot( NULL, xlim = xlim, ylim = ylim, 
                    xlab = xlab, ylab = ylab )
    if( LOG )suppressWarnings( plot( NULL, xlim = xlim, ylim = ylim, 
                                     xlab = xlab, ylab = ylab, log = 'y' ) )
  }
  
  if( PLOT )polygon( xbound, ybound, border = NA, col = .getColor( col, trans ) )
  
  invisible( cbind( xbound, ybound ) )
}

cleanMunicipality <- function(mname){
  
  require(stringi)
  
  mname <- as.character(mname)
  mname <- stringi::stri_trans_general(mname, "Latin-ASCII")
  mname <- .replaceString(mname, ' and ', '/')
  
  toRemove <- c('Community of ', 'District of ', 'Region of ')
  
  for(k in 1:length(toRemove)){
    mname <- .replaceString(mname, toRemove[k], '')
  }
  
  mname <- .replaceString(mname, ',', ', ')
  
  wspace <- which(startsWith(mname, ' '))
  if(length(wspace) > 0)mname[wspace] <- substr(mname[wspace], 2, 10000)
  
  mname
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
}

abb2state <- function(name, convert = T, strict = F, truncate = 100){
  
  require(stringi)
  
  data(state)
  
  # state data doesn't include DC
  state = list()
  state[['name']] = c(state.name,"District Of Columbia")
  state[['abb']] = c(state.abb,"DC")
  
  mm <- match(name, state$name)
  wf <- which(is.finite(mm))
  ss <- name
  ss[wf] <- state$abb[mm[wf]]
  
  ww <- which(nchar(ss) == 0 | is.na(ss) | nchar(ss) > truncate)
  
  # if(length(ww) > 0)ss[ww] <- missing
  
  if(length(ww) > 0){
    fix <- ss[ww]
    fix <- .replaceString(fix, '-', '')
    fix <- .replaceString(fix, ' ', '')
    ss[ww] <- substr(fix, 1, truncate)
  }
  stri_trans_general(ss, "Latin-ASCII")
}
countWhiteSpaces <- function(s){  # fix this
  inds <- gregexpr(" ", s)
  
  ll <- sapply( inds, length)
  
  wd   <- which(diff(inds) != 1)
  return( diff(c(0, wd, length(inds)) ) )
}
  


latlong2state <- function(pointsDF, level = 'state') {
  
  # level can be 'county'
  
  require(sp)
  require(maps)
  require(maptools)
  
  # Prepare SpatialPolygons object with one SpatialPolygon
  # per state (plus DC, minus HI & AK)
  states <- map(level, fill=TRUE, col="transparent", plot=FALSE)
  IDs <- sapply(strsplit(states$names, ":"), function(x) x[1])
  states_sp <- map2SpatialPolygons(states, IDs=IDs,
                                   proj4string=CRS("+proj=longlat +datum=WGS84"))
  
  # Convert pointsDF to a SpatialPoints object 
  pointsSP <- SpatialPoints(pointsDF, 
                            proj4string=CRS("+proj=longlat +datum=WGS84"))
  
  # Use 'over' to get _indices_ of the Polygons object containing each point 
  indices <- over(pointsSP, states_sp)
  
  # Return the state names of the Polygons object containing each point
  stateNames <- sapply(states_sp@polygons, function(x) x@ID)
  stateNames[indices]
}

.ls.objects <- function (pos = 1, pattern, order.by,
                         decreasing=FALSE, head=FALSE, n=5) {
  napply <- function(names, fn) sapply(names, function(x)
    fn(get(x, pos = pos)))
  names <- ls(pos = pos, pattern = pattern)
  obj.class <- napply(names, function(x) as.character(class(x))[1])
  obj.mode <- napply(names, mode)
  obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
  obj.size <- napply(names, object.size)
  if(!is.numeric(obj.size))return( numeric(0) )
  obj.size <- round(obj.size*1e-6, 2)
  
  obj.dim <- t(napply(names, function(x)
    as.numeric(dim(x))[1:2]))
  vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
  obj.dim[vec, 1] <- napply(names, length)[vec]
  
  out <- data.frame(obj.type, obj.size, obj.dim)
  names(out) <- c("Type", "Mb", "Rows", "Columns")
  
  if (!missing(order.by))
    out <- out[order(out[[order.by]], decreasing=decreasing), ]
  if (head)
    out <- head(out, n)
  out
}
# shorthand
lsos <- function(..., n=10) {
  
  # note: use pos = environment() to get non-global objects 
  .ls.objects(..., order.by="Mb", decreasing=TRUE, head=TRUE, n=n)
}

bigObjectsAll <- function(...,n = 20){
  
  envs <- search()
  for(k in 1:length(envs)){
    ktab <- .ls.objects(pos=k, order.by="Mb", decreasing=TRUE, head=TRUE, n=n)
    cat( paste('\n',envs[k], '\n') )
    print(ktab)
  }
}
        
        

bayesReg <- function(formula, data, ng = 5000, TOBIT=NULL){
  
  fc <- as.character(formula)
  
  yy <- unlist( strsplit( fc, '~' ) )
  yy <- yy[ nchar(yy) > 0]
  y  <- data[,yy[1]]
  
  ypos  <- which(y > 0)
  yzero <- which(y == 0)
  nzero <- length(yzero)
  npos  <- length(ypos)
  
  if(is.null(TOBIT)){
    TOBIT <- F
    if(nzero > 0)TOBIT <- T
  } 
  
  if(TOBIT)message('fitted as Tobit model')
  
  tmp <- model.frame(formula, data, na.action=NULL)
  x   <- model.matrix(formula, data=tmp)
  
  colnames(x)[1] <- 'intercept'
  
  xnames    <- colnames(x)
  snames    <- colnames(y)
  Q         <- ncol(x)
  n         <- nrow(x)
  predXcols <- 2:Q
  
  ymiss <- which(is.na(y))
  
  yy <- y
  xx <- x
  if(length(ymiss) > 0){
    yy <- yy[-ymiss]
    xx <- xx[-ymiss,]
  }
  
  XX  <- crossprod(xx)
  IXX <- solve(XX)
  bg  <- IXX%*%crossprod(xx,yy)
  py  <- x%*%bg
  w   <- y
  w[ymiss] <- mean(y,na.rm=T)
  if(TOBIT)w[w == 0] <- py[w == 0]
  wi  <- c(which(y == 0),ymiss)   
  XX  <- crossprod(x)
  
  
  priorB   <- matrix(0,Q,1)         #prior mean regression parameters
  priorIVB <- solve(diag(1000,Q))   #inverse prior covariance matrix
  s1       <- .1                    #variance prior values
  s2       <- .1
  
  bchains <- matrix(NA,ng,Q)
  schains <- rep(0,ng)  #store parameters
  colnames(bchains) <- xnames
  ychains <- matrix(0,ng,n)
  
  for(g in 1:ng){
    
    sg <- updateSigma(x, w, bg, s1, s2)
    bg <- updateBeta(x, w, sg, priorIVB, priorB, XX)
    mu <- x%*%bg
    py <- rnorm(n,mu,sqrt(sg))
    
    if(TOBIT){
      w[yzero] <- .tnorm(nzero,-Inf,0,mu[yzero],sqrt(sg))
      py[py <= 0] <- 0
    }
    
    w[ymiss] <- py[ymiss]
    
    bchains[g,] <- bg   #store estimates
    schains[g]  <- sg
    ychains[g,] <- py
  }
  
  beta <- signif( t( apply(bchains,2,quantile,c(.5,.025,.975)) ), 4)
  zero <- which(beta[,2] < 0 & beta[,3] > 0)
  notZero <- rep('*',Q)
  notZero[ zero ] <- ' '
  
  colnames(beta) <- c('median','0.025','0.975')
  
  py <- signif( t( apply( ychains, 2, quantile,c(.5,.025,.975) ) ), 3)
  py[,1] <- colMeans(ychains)
  py <- cbind(y,py)
  
  rmspe <- sqrt( mean((y - py[,2])^2, na.rm=T) )
  
  list(beta = beta, predict = py, sigma = median(schains), 
       rmspe = rmspe)
}

updateSigma <- function(x, y, beta, s1, s2){ # random value for residual variance
  n  <- length(y)
  u1 <- s1 + n/2
  u2 <- s2 + 0.5*crossprod(y - x%*%beta)
  1/rgamma(1,u1,u2)                          
}

updateBeta <- function(x, y, sigma, priorIVB, priorB, 
                       XX=NULL){  # random vector of coefficients
  
  if(is.null(XX))XX <- crossprod(x)
  
  V  <- solve( XX/sigma + priorIVB ) 
  v  <- crossprod(x,y)/sigma + priorIVB%*%priorB
  t( myrmvnorm(1,t(V%*%v),V) )                     
}

cols2rows <- function(data, dname, rname, cname, 
                      rows = NULL, columns = NULL){
  
  # column from data to fill in matrix
  # rname, cname - columns in data for row and column names
  # rname can be a subset of new rows, given as 'rows'
  # cname can be a subset of new columns, given as 'columns'
  
  id  <- data[,rname]
  ids <- sort(unique(id))
  jd  <- data[,cname]
  jds <- sort(unique(jd))
  
  if(!is.null(rows))   ids <- sort( unique( c(ids, rows) ) )
  if(!is.null(columns))jds <- sort( unique( c(jds, columns) ) )
  
  
  dmat <- matrix(NA, length(ids), length(jds))
  colnames(dmat) <- jds
  rownames(dmat) <- ids
  
  ij <- cbind( match(id, ids), match(jd, jds) )
  wf <- which(is.finite(ij[,1]) & is.finite(ij[,2]))
  
  # dmat[wf,][ ij[wm,] ] <- tdata$diamRC[wf[wm]]
  
  dmat[ ij[wf,] ] <- data[wf,dname]
  dmat
}

rows2cols <- function(dmat, sep='_'){
  
  # from i by j matrix to ij vector
  
  id <- rownames(dmat)
  jd <- colnames(dmat)
  
  ii <- rep(id, length(jd))
  jj <- rep(jd, each = nrow(dmat))
  
  ij <- cbind( match(ii, id), match(jj, jd) )
  
  rr <- columnPaste(ii, jj, sep)
  
  dvec <- dmat[ ij ]
  names(dvec) <- rr
  dvec
}

matchFinite <- function(c1, c2){
  
  # returns matches in mm and which are finite in wf
  
  c1 <- as.character(c1)
  c2 <- as.character(c2)
  
  mm <- match(c1, c2)
  wf <- which(is.finite(mm))
  wn <- which(!is.finite(mm))
  
  if(length(wn) > 0){
    wm <- paste0( unique(c1[wn]), collapse = ', ')
    print( paste( 'missing names:',wm ) )
  }
  list(mm = mm, wf = wf)
}


USDA2traitTable <- function(usdaCode = NULL, code4 = NULL,
                            ufile = '/Volumes/research/clark/clark.unix/USDAplants.csv',
                            tfile = '/Volumes/research/clark/clark.unix/traitTable.csv'){
  
  # update traitTable.csv with traitNew.csv that includes usdaCode
  
  # either usdaCode or code4 must be provided
  
  utab <- read.csv(ufile, stringsAsFactors = F)   # USDA file
  ttab <- read.csv(tfile, stringsAsFactors = F)   # trait file
  colnames(utab)[colnames(utab) == 'taxonID'] <- 'codeUSDA'
  
  if(!is.null(usdaCode)){
    
    mm <- match(usdaCode, utab$codeUSDA)
    
  }else{
    mm <- match(code4, utab$code4)
  }
  wf <- which(is.finite(mm))
    
  newSpecs <- data.frame( matrix(NA, length(wf), ncol(ttab)), stringsAsFactors = F)
  colnames(newSpecs) <- colnames(ttab)
  newSpecs[, colnames(utab)] <- utab[mm[wf],]
  
  c1 <- lowerFirstLetter( substr(newSpecs$genus, 1, 6) )
  c2 <- upperFirstLetter( substr(newSpecs$species, 1, 6) )
  newSpecs$code6 <- columnPaste(c1,c2,'')
  
  ttab <- rbind(newSpecs, ttab)
  ttab <- ttab[order(ttab$genus, ttab$species),]
  
  ofile <- tfile
  ofile <- .replaceString(ofile, 'traitTable', 'traitNew')
  write.csv(ttab, file = ofile, quote = F, row.names = F)  # updated file
  
  print(' check traitNew.csv, if ok, replace traitTable.csv')
  
}




aggData <- function(cnames, data, gather, FUN){
  
  #cnames - column names in data to operate on
  #gather   - list of indices, all of same length
  
  cnames <- cnames[cnames %in% colnames(data)]
  if(length(cnames) == 0)return( numeric(0) )
  
  FAC <- F
  df <- data[,cnames]
  
  if(FUN == 'factor'){
    FAC <- T
    tf <- fnames <- character(0)
    df <- numeric(0)
    
    for(j in 1:length(cnames)){
      kf <- as.character(data[,cnames[j]])
      tf <- cbind(tf, kf)
      fnames <- append(fnames,list(sort(unique(kf))))
      df <- cbind(df, match(kf,fnames[[j]]) )
    }
    colnames(df) <- cnames
    FUN <- 'max'
  }
  
  tmp <- aggregate(df, by=gather, FUN=FUN)
  #  ord <- do.call(order,list(tmp[,names(gather)]))
  colnames(tmp)[-c(1:length(gather))] <- cnames
  
  if(FAC){
    for(j in 1:length(cnames)){
      kf <- fnames[[j]][tmp[,cnames[j]]]
      tmp[,cnames[j]] <- as.factor(kf)
    }
  }
  if(length(tmp) == 0)tmp <- numeric(0)
  tmp
}

.myBy <- function(x, i, j, summat=matrix(0,max(i),max(j)), 
                  totmat=summat, fun='mean'){  
  
  nn <- length(x)
  if( nn != length(i) | nn != length(j) )
    stop('vectors unequal in byFunctionRcpp')
  if( nrow(summat) < max(i) | ncol(summat) < max(j) )
    stop('matrix too small')
  
  ww <- which(is.na(x))
  if(length(ww) > 0){
    x <- x[-ww]
    i <- i[-ww]
    j <- j[-ww]
  }
  
  frommat <- cbind(i,j,x)
  
  nr  <- nrow(frommat)
  
  maxmat <- summat*0 - Inf
  minmat <- summat*0 + Inf
  
  tmp <- byRcpp(nr, frommat, totmat, summat, minmat, maxmat)
  
  if(fun == 'sum')return(tmp$sum)
  if(fun == 'mean'){
    mu <- tmp$sum/tmp$total
    mu[is.na(mu)] <- 0
    return(mu)
  }
  if(fun == 'min'){
    return( tmp$min )
  }
  tmp$max
}

fillTraits <- function(specCodes){
  
  # construct a species by trait data.frame for a character vector of species codes specCodes
  
  source('/nfs/clark/clark.unix/GJAMfunctions/gjamHfunctions.R')
  
  #specCodes: six-code, four-code, or USDA format
  
  traits <- read.table('../traitTable.txt',header=T)
  
  code <- 'code6'                                     #GGGGGGSSSSSS
  if(max(nchar(specCodes)) == 8)code <- 'code4'       #GGGGSSSS
  if(min(nchar(specCodes)) == 4)code <- 'fiaCode'     #GGSS
  
  id <- c(code,"family","genus","species")
  ocols <- c('shade','drought','flood')               #ordinal variables
  
  cnames <- colnames(traits)
  
  imat <- traits[,id]
  
  tvals <- traits[,!cnames %in% id]
  
  facts <- which(sapply(tvals,is.factor))
  facts <- facts[!names(facts) %in% c('code4','fiaCode')]
  logic <- which(sapply(tvals,is.logical))
  numer <- which(sapply(tvals,is.numeric))
  
  tnum <- as.matrix( tvals[,numer] )
  
  # genus table
  genera <- sort(unique(traits[,'genus']))                               #numeric
  i <- match(traits[,'genus'],genera)
  ii <- rep(i, ncol(tnum))
  jj <- rep(1:ncol(tnum),each = nrow(tnum))
  tmp <- signif( .byGJAM(as.vector(as.matrix(tnum)),ii,jj,fun='mean'), 4)
  rownames(tmp) <- attr(genera, 'levels')
  colnames(tmp) <- colnames(tnum)
  tmp[,ocols] <- round(tmp[,ocols])
  genTable <- tmp
  
  gfac <- data.frame( matrix(0,nrow(genTable),length(facts) ) ) #factors
  k <- 1
  for(j in facts){
    tj <- table(tvals[,j], traits[,'genus'])
    jt <- rownames(tj)[apply(tj,2,which.max)]
    jt <- factor(jt)
    gfac[,k] <- jt
    k <- k + 1
  }
  colnames(gfac) <- names(facts)
  
  glog <- data.frame( matrix(0,nrow(genTable),length(logic) ) ) #logical
  k <- 1
  for(j in logic){
    tj <- table(tvals[,j], traits[,'genus'])
    jt <- rownames(tj)[apply(tj,2,which.max)]
    jt <- as.logical(jt)
    glog[,k] <- jt
    k <- k + 1
  }
  colnames(glog) <- names(logic)
  
  genTable <- cbind(gfac, glog, genTable)
  
  write.table(genTable, file='../traitGenera.txt', 
              row.names = T, quote=F, sep='\t')
  
  # family table
  family <- sort(unique(traits[,'family']))
  i <- match(traits[,'family'],family)
  ii <- rep(i, ncol(tnum))
  jj <- rep(1:ncol(tnum),each = nrow(tnum))
  tmp <- signif( .byGJAM(as.vector(as.matrix(tnum)),ii,jj,fun='mean'), 4)
  rownames(tmp) <- attr(family, 'levels')
  colnames(tmp) <- colnames(tnum)
  tmp[,ocols] <- round(tmp[,ocols])
  famTable <- tmp
  
  gfac <- data.frame( matrix(0,nrow(famTable),length(facts) ) ) #factors
  k <- 1
  for(j in facts){
    tj <- table(tvals[,j], traits[,'family'])
    jt <- rownames(tj)[apply(tj,2,which.max)]
    jt <- factor(jt)
    gfac[,k] <- jt
    k <- k + 1
  }
  colnames(gfac) <- names(facts)
  
  glog <- data.frame( matrix(0,nrow(famTable),length(logic) ) ) #logical
  k <- 1
  for(j in logic){
    tj <- table(tvals[,j], traits[,'family'])
    jt <- rownames(tj)[apply(tj,2,which.max)]
    jt <- as.logical(jt)
    glog[,k] <- jt
    k <- k + 1
  }
  colnames(glog) <- names(logic)
  
  famTable <- cbind(gfac, glog, famTable)
  
  write.table(famTable, file='../traitFamilies.txt', 
              row.names = T, quote=F, sep='\t')
  
  ################## generate table for specCodes
  
  wm <- match(specCodes,traits[,code])
  wf <- which(is.finite(wm))
  nt <- length(specCodes)
  
  gtab <- matrix(NA, nt, length(numer))
  rownames(gtab) <- specCodes
  colnames(gtab) <- names(numer)
  
  wm <- match(specCodes, traits[,code])
  wf <- which(is.finite(wf))
  gtab[wf,] <- tnum[wm[wf],names(numer)]
  
  wna <- which(is.na(gtab),arr.ind=T)                     #genus
  if(length(wna) > 0){
    wm  <- match(specCodes[wna[,1]], traits[,code])
    wf  <- which(is.finite(wm))
    gen <- as.character( traits[wm[wf],'genus'] ) 
    gd <-  genTable[gen, names(numer)]
    gtab[wna[wf,]] <- gd[ cbind(1:length(wf),wna[wf,2]) ]
  }
  
  wna <- which(is.na(gtab),arr.ind=T)                     #family
  if(length(wna) > 0){
    wm  <- match(specCodes[wna[,1]], traits[,code])
    wf  <- which(is.finite(wm))
    fam <- as.character( traits[wm[wf],'family'] ) 
    gd <-  famTable[fam, names(numer)]
    gtab[wna[wf,]] <- gd[ cbind(1:length(wf),wna[wf,2]) ]
  }
  
  gfac <- data.frame( matrix(0,nt,length(facts) ) ) #factors
  k <- 1
  
  for(j in facts){
    tj <- factor(rep(NA,nt) , levels = attr(tvals[[j]],'levels') )
    tj[wf] <- tvals[wm[wf],j]
    
    wna  <- which(is.na(tj))
    wtab <- match( specCodes[wna] , traits[,code] )
    if( !all(is.na(wtab)) ){             # try genus
      wgen <- wtab[ !is.na(wtab) ]
      gen  <- traits[wgen,'genus']
      gval <- genTable[gen,colnames(tvals)[j]]
      tj[wna[!is.na(wtab)]] <- gval
    }
    
    wna  <- which(is.na(tj))
    wtab <- match( specCodes[wna] , traits[,code] )
    if( !all(is.na(wtab)) ){             # try genus
      wgen <- wtab[ !is.na(wtab) ]
      gen  <- traits[wgen,'family']
      gval <- famTable[gen,colnames(tvals)[j]]
      tj[wna[!is.na(wtab)]] <- gval
    }
    gfac[,k] <- tj
    k <- k + 1
  }
  colnames(gfac) <- names(facts)
  
  glog <- data.frame( matrix(0,nt,length(logic) ) ) #factors
  k <- 1
  
  for(j in logic){
    
    tj <- rep(F,nt)
    wm <- match(specCodes,traits[,code])
    wf <- which(is.finite(wm))
    tj[wf] <- traits[wm[wf],names(logic)[k]]
    
    wna  <- which(is.na(tj))
    wtab <- match( specCodes[wna] , traits[,code] )
    if( !all(is.na(wtab)) ){             # try genus
      wgen <- wtab[ !is.na(wtab) ]
      gen  <- traits[wgen,'genus']
      gval <- genTable[gen,names(logic)[k]]
      tj[wna[!is.na(wtab)]] <- gval
    }
    
    wna  <- which(is.na(tj))
    wtab <- match( specCodes[wna] , traits[,code] )
    if( !all(is.na(wtab)) ){             # try genus
      wgen <- wtab[ !is.na(wtab) ]
      gen  <- traits[wgen,'family']
      gval <- famTable[gen,names(logic)[k]]
      tj[wna[!is.na(wtab)]] <- gval
    }
    glog[,k] <- tj
    k <- k + 1
  }
  colnames(glog) <- names(logic)
  
  rownames(gtab) <- NULL
  
  cbind(specCodes, gfac, glog, gtab)
  
}



extractRasterSeries <- function(sites, dataFolder){
  
  require(rgdal)
  
  tifFiles <- dir(dataFolder, pattern = '*.tif$', full.names = T)
  tifFilesshort <- dir(dataFolder, pattern = '*.tif$', full.names = F)
  if(length(tifFiles)==0) stop('No .tif file was found!')
  
  nt <- length(tifFiles)
  ns <- nrow(sites)  
  
  extMat <- matrix(NA, nrow = ns, ncol = nt)
  colnames(extMat) <- 1:nt
  rownames(extMat) <- sites[,'name']
  
  for(i in 1:nt){
    r <- raster(tifFiles[i])
    extMat[,i] <- extract(r, sites[,c('lon','lat')])
  }
  
  extMat
}



replaceNonAscii <- function(filename){
  l   <- 'latin1'
  a   <- 'ASCII'
  x   <- readLines(filename)
  nar <- grep("notASCII", iconv(x, l, a, sub="notASCII"))
  if(length(nar) > 0) x[nar] <- iconv(x[nar], l, a, sub="-")
  list(notAsciRows = nar, x = x)
}


covarGrid <- function(lonlat, path='../dataFiles/'){
  
  require(RANN)
  file <- paste(path,'UScovarGrid.rdata',sep='')
  load(file)
  
  ii <- nn2(gridded.data[,c('lon','lat')], lonlat,k=1)[[1]]
  xx <- gridded.data[ii,c('soil','temp','prec','therm','def','nlcd')]
  colnames(xx)[colnames(xx) == 'temp'] <- 'winterTemp'  # two temp variables
  xx
}

.insertRows <- function(xnow,xinsert,ibefore){
  
  # length(ibefore) = nrow(xinsert)
  
  nnow <- nrow(xnow)
  inew <- sort( c(1:nnow,ibefore) )
  nnew <- length(inew)
  xnew <- as.data.frame( matrix(NA,nnew,ncol(xnow)) )
  
  wi   <- which(duplicated(inew)) - 1
  wm   <- c(1:nnew)[-wi]
  
  for(j in 1:ncol(xnow)){
    xj <- xnow[,j]
    if( is.factor(xj) ){
      levs <- attr(xj,'levels')
      xm <- factor(xnew[[j]], levels=levs)
      xm[wm] <- xnow[,j]
      xm[wi] <- xinsert[,j]
      xnew[[j]] <- xm
    } else {
      xnew[wm,j] <- xnow[,j]
      xnew[wi,j] <- xinsert[,j]
    }
  }
  colnames(xnew) <- colnames(xnow)
  xnew
}

.combineFacLevels <- function(xfactor,fname=NULL, aname = 'reference', 
                              vminF=1){
  tmp <- as.character(xfactor)
  tmp[tmp %in% fname] <- aname
  tab <- table(tmp)
  wm  <- names(tab)[tab < vminF]
  tmp[tmp %in% wm] <- aname
  as.factor(tmp)
}

.factor2Numeric <- function(xfactor) as.numeric(as.character(xfactor))

.getColor <- function(col,trans){
  
  # trans - transparency fraction [0, 1]
  
  tmp <- col2rgb(col)
  rgb(tmp[1,], tmp[2,], tmp[3,], maxColorValue = 255, 
      alpha = 255*trans, names = paste(col,trans,sep='_'))
}

getColorNames <- function(cols){
  
  # takes color codes from colorRampPalette
  
  xcol <- t( col2rgb(cols) )
  acol <- t(  col2rgb(colors()) )
  
  xd <- 0
  
  for(k in 1:3){
    xt <- xcol[,k]
    xs <- acol[,k]
    xd <- xd + outer(xt,xs,function(xt,xs) (xt - xs)^2)
  }
  
  xd <- t(sqrt(xd))
  dr <- apply(xd,2,which.min)
  colors()[dr]
}


.insetPlot <- function(func, x1, y1, xwide, ywide, 
                       bgCol = 'white'){
  
  pnow <- par('plt')
  unow <- par('usr')
  on.exit( par(plt = pnow, usr = unow) )
  dplt <- c(pnow[2] - pnow[1], pnow[4] - pnow[3])
  dusr <- c(unow[2] - unow[1], unow[4] - unow[3])
  
  posx <- pnow[1] + (x1 - unow[1])*dplt[1]/dusr[1]
  posy <- pnow[3] + (y1 - unow[3])*dplt[2]/dusr[2]
  
  px2 <- pnow[1] + (x1 + xwide - unow[1])*dplt[1]/dusr[1]
  py2 <- pnow[3] + (y1 + ywide - unow[3])*dplt[2]/dusr[2] 
  
  par(plt = c(posx, px2, posy, py2), new = TRUE)
  plot.new()
  polygon(c(-0.1, 1.1, 1.1, -0.1), c(-0.1, -0.1, 1.1, 1.1), 
          border = NA, col = bgCol)
  par(plt = c(posx, px2, posy, py2), new = TRUE)
  
  eval(func)
}

getEDF <- function(events){
  
  # events are values at which events occur
  
  ord <- order(events)
  xp  <- events[ord]/max(events)
  yp  <- c(1:length(events))/length(events)
  cbind(xp,yp)
}

.combineFacLevels <- function(xfactor,fname=NULL, aname = 'reference', 
                              vminF=1){
  
  tmp <- as.character(xfactor)
  tmp[tmp %in% fname] <- aname
  tab <- table(tmp)
  wm  <- names(tab)[tab < vminF]
  tmp[tmp %in% wm] <- aname
  as.factor(tmp)
}

.factor2Numeric <- function(xfactor){
  
  as.numeric(as.character(xfactor))
}

.getColor <- function(col,trans){
  
  # trans - transparency fraction [0, 1]
  
  tmp <- col2rgb(col)
  rgb(tmp[1,], tmp[2,], tmp[3,], maxColorValue = 255, 
      alpha = 255*trans, names = paste(col,trans,sep='_'))
}

.bins4data <- function(obs, nPerBin=NULL, breaks=NULL, nbin=NULL, LOG=F){
  
  if( is.null(breaks) ){
    
    if(is.null(nbin))nbin <- 20
    
    br   <- range(obs[is.finite(obs)],na.rm=T)
    bins <- seq(br[1],br[2],length=nbin)
    if(LOG){
      yy <- obs[obs > 0]
      oo <- min( yy,na.rm=T )
      
      nper <- length(yy)/nbin
      nbb <- nper/length(yy)
      nbb <- seq(0,1,length=100)
      if(max(nbb) < 1)nbb <- c(nbb,1)
      bins <- 10^quantile(log10(yy),nbb,na.rm=T)
      bins <- sort(unique(bins))
      br[1] <- .5*oo
      
    }
    if( !is.null(nPerBin) ){
      nbb <- nPerBin/length(obs)
      nbb <- seq(0,1,by=nbb)
      if(max(nbb) < 1)nbb <- c(nbb,1)
      bins <- quantile(obs,nbb,na.rm=T)
      bins <- sort(unique(bins))
      
      db <- diff(bins)
      wb <- which( db/diff(range(quantile(obs,c(.1,.9)))) < .02)
      wb <- wb[wb != 1]
      if(length(wb) > 0)bins <- bins[-wb]
      
      nbin <- length(bins)
    }
  } else {
    bins <- breaks
    nbin <- length(bins)
  }
  
  list(breaks = breaks, bins = bins, nbin = nbin)
}

.plotObsPred <- function(obs,yMean,ySE=NULL,nbin=NULL,nPerBin=NULL,breaks=NULL,
                         LOG=F,xlimit=NULL,ylimit=NULL,xlabel='Observed',
                         ylabel='Predicted',
                         ptcol=NULL, boxPerc = .6826895, whiskerPerc = .95,
                         fill=NULL, add=F, box.col='black', wide=NULL, POINTS=T,
                         MEDIAN=T){
  
  aa <- (1 - boxPerc)/2
  boxQuant <- c(aa, 1 - aa )
  aa <- (1 - whiskerPerc)/2
  whiskerQuant <- c(aa, 1 - aa )
  
  if(is.null(ptcol)){
    ptcol <- 'black'
    ptcol <- 'grey'
    if(!is.null(nbin))ptcol <- 'grey'
  }
  if(length(ptcol) == 1)ptcol <- rep(ptcol,length(obs))
  
  if(is.null(xlimit))xlimit <- quantile(obs[is.finite(obs)],c(.01,.99),na.rm=T)
  if(is.null(ylimit))ylimit <- range(yMean[is.finite(yMean)],na.rm=T)
  
  xxx <- obs
  yyy <- yMean
  
  if(LOG){
    if(is.null(xlimit))xlimit <- range( obs[obs > 0],na.rm=T )
    if(is.null(ylimit))ylimit <- range( yMean[yMean > 0],na.rm=T )
    if(xlimit[1] <= 0)xlimit[1] <- .001
  }
  
  if(!POINTS){
    xxx <- xlimit[1]
    yyy <- ylimit[1]
  }
  
  if(!add){
    if(is.null(ylimit)){
      if(!LOG)plot(xxx,yyy,col=ptcol,cex=.3,xlab=xlabel,ylab=ylabel)
      if(LOG) plot(xxx,yyy,col=ptcol,cex=.3,xlab=xlabel,ylab=ylabel,log='xy')
    }
    if(!is.null(ylimit)){
      if(!LOG)plot(xxx,yyy,col=ptcol,cex=.3,xlab=xlabel,ylab=ylabel,
                   xlim=xlimit,ylim=ylimit)
      if(LOG) plot(xxx,yyy,col=ptcol,cex=.3,xlab=xlabel,ylab=ylabel,
                   xlim=xlimit,log='xy',ylim=ylimit)
    }
  }
  
  if(!is.null(ySE)){
    ylo <- yMean - 1.96*ySE
    yhi <- yMean + 1.96*ySE
    for(i in 1:length(obs))lines(c(obs[i],obs[i]),c(ylo[i],yhi[i]),
                                 col='grey',lwd=2)
  }
  
  tmp    <- .bins4data(obs,nPerBin=nPerBin,breaks=breaks,LOG=LOG)
  breaks <- tmp$breaks
  bins   <- tmp$bins
  nbin   <- tmp$nbin
  
  if(is.null(wide))wide <- diff(bins)/2.1
  if(length(wide) == 1)wide <- rep(wide,nbin)
  minmax <- par('usr')[1:2]
  dff    <- diff(minmax)
  if(!LOG)wide[wide > dff/5] <- dff/5
  
  maxx <- 0
  last <- F
  
  for(k in 1:(nbin-1)){
    
    mb <- bins[k+1]
    if(mb >= xlimit[2]){
      last <- T
      mb   <- Inf
    }
    ok <- which(obs >= bins[k] & obs < mb)
    if(length(ok) == 0)next
    qk <- which(is.finite(yMean) & obs >= bins[k] & obs <= mb)
    q  <- quantile(yMean[qk],c(.5,whiskerQuant[1],boxQuant[1],
                               boxQuant[2],whiskerQuant[2]),na.rm=T)
    if(LOG)q[q <= 0] <- ylimit[1]
    
    ym <- q[1]
    xx <- mean(bins[k:(k+1)])      # bounded by bins
    if(!LOG){
      #    if(!is.null(nPerBin) & MEDIAN)xx <- median(obs[ok],na.rm=T)
      if(MEDIAN)xx <- median(obs[ok],na.rm=T)
    } else {
      xx <-  sqrt( prod(bins[k:(k+1)]) )
    }
    points(xx,q[1],pch=3,col=box.col)
    yy    <- q[c(2,5)]
    yy[1] <- max(c(yy[1],ylimit[1]),na.rm=T) + .0000001
    yy[2] <- max(yy)
    
    yy1 <- q[3]
    yy1 <- max(yy1,ylimit[1],na.rm=T) + .00000001
    yy2 <- max(yy1,q[4])
    
    minx <- xx - .4*(xx - bins[k])
    maxx <- xx + .4*(mb - xx)
    
    figRange <- par('usr')
    totalwide <- (maxx - minx)/diff(figRange[1:2])
    
    if(is.null(nPerBin)){
      
      if(maxx >= xlimit[2])maxx <- xlimit[2]
      
      if(LOG){
        
        dx   <- log10(bins[k+1]) - log10(xx)
        maxx <- 10^(log10(xx) + .2*dx)
        if(k == 1){
          dx <- -log10(xlimit[1]) + log10(xx)
        } else {
          dx <- -log10(bins[k-1]) + log10(xx)
        }
        minx <- 10^(log10(xx) - .2*dx)
        if(minx < xlimit[1])minx <- xlimit[1]
        totalwide <- (log10(maxx) - log10(minx))/diff(figRange[1:2])
      }
      
      
      rect(minx,yy1,maxx,yy2,col=fill,border=box.col)
      lines(c(minx,maxx),c(ym,ym),lwd=2,col=box.col)
      xx <- mean(c(minx,maxx))
    }
    
    if(!is.null(nPerBin)){
      qo <- quantile(obs[ok],c(.3,.7,.25,.75),na.rm=T)
      if(qo[1] == qo[2] | !MEDIAN)qo <- c(xx-.2*wide[k],
                                          xx+.2*wide[k],xx-.3*wide[k],
                                          xx+.3*wide[k])
      rect(qo[1],yy1,qo[2],yy2,col=fill,border=box.col)
      lines(c(qo[3],qo[4]),c(ym,ym),lwd=2,col=box.col)
      lines(rep(mean(qo[1:2]),2),yy,lwd=2,col=box.col)
    }
    lines(c(xx,xx),yy,lwd=2,col=box.col)
    if(last)break
  }
  
  if(POINTS)points(xxx,yyy,col=ptcol)
  
  invisible( bins )
}
.byIndex <-function(xx,INDICES,FUN,coerce=F,...){  
  
  #INDICES is list, each same length as  x
  
  #  fun <- match.fun(FUN)
  
  nl <- length(INDICES)
  
  tmp  <-  unlist(by( as.vector(xx),INDICES,FUN,...) ) 
  nd   <- dim(tmp)
  tmp  <- array(tmp,dim=nd, dimnames=dimnames(tmp))
  
  tmp[is.na(tmp)] <- 0
  
  if(!coerce)return(tmp)
  
  dname <- dimnames(tmp)
  mk    <- rep(0,length(nd))
  
  for(k in 1:length(nd))mk[k] <- max(as.numeric(dimnames(tmp)[[k]]))
  
  if(length(nd) == 2){
    newr <- c(1:mk[1])
    newc <- c(1:mk[2])
    tnew <- array(0,dim=mk)
    rownames(tnew) <- newr
    colnames(tnew) <- newc
    tnew[rownames(tmp),colnames(tmp)] <- tmp
  }
  
  tmp
}

#################################gjam functions

.directIndirectCoeffs <- function( snames, xvector, chains,
                                   MEAN=T, factorList = NULL,
                                   keepNames, omitY, standardX, 
                                   sdScaleY, nsim = 50,
                                   otherpar = NULL, REDUCT = F, ng, burnin){
  
  # if MEAN, then use means, otherwise median
  # indirect do not change with x, can choose not to calculate
  #          - a list of vectors, one for each multilevel factor, 
  #            where hostNames appear in colnames of bchain
  #indirFrom - effect from all others
  #indirTo   - effect on all others
  
  if(is.matrix(xvector))stop('xvector must be a row vector with variable names')
  xnames <- names(xvector)
  
  N      <- otherpar$N
  r      <- otherpar$r
  bchain <- chains$bgibbs
  schain <- chains$sgibbs
  sigErrGibbs <- kchain <- NULL
  if(REDUCT){
    kchain      <- chains$kgibbs
    sigErrGibbs <- chains$sigErrGibbs
  }
  
  ns <- nsim
  simIndex <- sample(burnin:ng,ns,replace=T)
  
  if(sdScaleY){
    tmp <- .expandSigmaChains(snames, simIndex = simIndex, schain, 
                              sigErrGibbs, kchain, 
                              bchain, otherpar, REDUCT, CHAINS = T)
    bchain <- tmp$chainList$bchainCor    # standardized
    kchain <- kchain[simIndex,]
    schain <- schain[simIndex,]          # not standardized
    sigErrGibbs <- sigErrGibbs[simIndex]
  } else {
    bchain <- bchain[simIndex,]
    schain <- schain[simIndex,]
  }
  
  if(length(factorList) > 0){
    factorNames <- factorList
    for(j in 1:length(factorList)){
      tmp <- matrix( unlist(strsplit(factorList[[j]],names(factorList)[j])),
                     ncol=2,byrow=T)[,2]
      tmp[nchar(tmp) == 0] <- paste(names(factorList)[j],c(1:length(tmp)),sep='')
      factorNames[[j]] <- tmp
    }
  }
  
  S <- S1 <- length(snames)
  sindex <- c(1:S)
  knames <- snames
  
  nc <- nrow(bchain)
  
  gs <- 1:nrow(bchain)
  
  if(length(omitY) > 0){
    wob <- grep(omitY,colnames(bchain))
    bchain[,wob] <- 0
    sindex <- sindex[!snames %in% omitY]
    knames <- snames[sindex]
    S1     <- length(knames)
  }
  
  nspec <- length(snames)
  
  ww   <- grep(':',xnames)
  main <- xnames
  if(length(ww) > 0)main <- xnames[-ww]
  main <- main[main != 'intercept']
  int  <- unique( unlist( strsplit(xnames[ww],':') ) ) 
  
  mainEffect <- matrix(NA,nspec,length(main))
  colnames(mainEffect) <- main
  rownames(mainEffect) <- snames
  intEffect  <- dirEffect <- indEffectTo <- mainEffect
  mainSd <- dirSd <- intSd <- indSdTo <- mainEffect 
  
  maxg <- length(main)*length(sindex)*length(gs)
  pbar <- txtProgressBar(min=1,max=maxg,style=1)
  ig   <- 0
  
  for(j in 1:length(main)){
    
    ttt <- .interactionsFromGibbs(mainx=main[j],bchain=bchain,
                                  specs=snames,xmnames=names(xvector),xx=xvector,
                                  omitY = omitY)
    maine   <- ttt$main
    inter   <- ttt$inter  #
    indirTo <- maine*0
    direct  <- maine + inter
    
    if(MEAN){
      dmain  <- colMeans(maine)
      inte   <- colMeans(inter)
      dir    <- colMeans(direct)
    } else {
      dmain  <- apply(maine,2,median)
      inte   <- apply(inter,2,median)
      dir    <- apply(direct,2,median)
    }
    
    mainEffect[sindex,j] <- dmain
    intEffect[sindex,j]  <- inte
    dirEffect[sindex,j]  <- dir
    
    mainSd[sindex,j] <- apply(maine,2,sd)
    intSd[sindex,j]  <- apply(inter,2,sd)
    dirSd[sindex,j]  <- apply(direct,2,sd)
    
    for(g in gs){
      
      if(REDUCT){
        Z  <- matrix(schain[g,],N,r)
        ss <- .expandSigma(sigErrGibbs[g], S, Z = Z, kchain[g,], 
                           REDUCT = T)[sindex,sindex]
        if(sdScaleY)cc <- .cov2Cor(ss)
      } else {
        ss <- .expandSigma(schain[g,], S = S, REDUCT = F)[sindex,sindex]
        if(sdScaleY)cc <- .cov2Cor(ss)
      }
      
      for(s in sindex){
        
        if(REDUCT){
          si <- .invWoodburryArma(sigErrGibbs[g], Z[kchain[g,sindex[-s]],])
          if(sdScaleY){
            dc <- diag(sqrt(diag(ss)))[-s,-s]
            ci <- dc%*%si%*%dc
          }
        } else {
          si <- solve(ss[-s,-s])
          if(sdScaleY)ci <- solve(cc[-s,-s])
        }
        
        if(!sdScaleY){
          sonk <- ss[drop=F,s,-s]
          e2   <- sonk%*%si%*%direct[g,-s]
        } else {
          sonk <- cc[drop=F,s,-s]
          e2   <- sonk%*%ci%*%direct[g,-s]
        }
        indirTo[g,s] <- e2
        
        ig <- ig + 1
        setTxtProgressBar(pbar,ig)
        
      } ##############
    }
    
    if(MEAN){
      indirectTo   <- colMeans(indirTo[gs,])
    } else {
      indirectTo   <- apply(indirTo[gs,],2,median)
    }
    indEffectTo[sindex,j]   <- indirectTo
    indSdTo[sindex,j]       <- apply(indirTo[gs,],2,sd)
  } ######################################
  
  if(!is.null(keepNames)){
    wk <- which(rownames(mainEffect) %in% keepNames)
    mainEffect <- mainEffect[wk,]
    intEffect <- intEffect[wk,]
    dirEffect <- dirEffect[wk,]
    indEffectTo   <- indEffectTo[wk,]
    mainSd    <- mainSd[wk,]
    dirSd     <- dirSd[wk,]
    indSdTo   <- indSdTo[wk,]
  }
  
  if(!is.null(standardX)){
    wx <- match(colnames(mainEffect),names(standardX))
    sx <- matrix(standardX[wx],nrow(mainEffect),length(wx),byrow=T)
    mainEffect <- mainEffect*sx
    intEffect  <- intEffect*sx
    dirEffect  <- dirEffect*sx
    indEffectTo <- indEffectTo*sx
  }
  
  list(mainEffect = mainEffect, intEffect = intEffect, dirEffect = dirEffect,
       indEffectTo = indEffectTo, mainSd = mainSd, dirSd = dirSd,
       intSd = intSd, indSdTo = indSdTo)
}

.interactionsFromGibbs <- function(mainx,bchain,specs,xmnames=names(xx),
                                   xx=colMeans(xx),omitY=NULL){
  
  # returns main effects and interactions for variable named main
  # xx are values of covariates to condition on
  # mainx is the name of a main effect
  
  if(length(omitY) > 0){
    wob <- grep(omitY,colnames(bchain))
    bchain[,wob] <- 0
    specs <- specs[!specs %in% omitY]
  }
  
  ww   <- grep(':',xmnames)
  int  <- unique( unlist( strsplit(xmnames[ww],':') ) ) 
  int  <- int[int != mainx]
  
  xj <- paste(mainx,specs,sep='_')
  wj <- which(colnames(bchain) %in%  xj)
  if(length(wj) == 0){
    xj <- paste(specs,mainx,sep='_')
    wj <- which(colnames(bchain) %in%  xj)
  }
  
  maine <- bchain[,xj]
  inter  <- maine*0
  
  m1 <- paste(mainx,':',sep='')
  m2 <- paste(':',mainx,sep='')
  i1 <- grep( m1,xmnames )
  i2 <- grep( m2,xmnames )
  
  if( length(i1) > 0 ){
    
    ww <- match(unlist( strsplit(xmnames[i1],m1) ),xmnames)
    ox <- xmnames[ww[is.finite(ww)]]
    for(kk in 1:length(i1)){
      xi <- paste(xmnames[i1[kk]],specs,sep='_')
      wi <- which(colnames(bchain) %in%  xi)
      if(length(wi) == 0){
        xi <- paste(specs,xmnames[i1[kk]],sep='_')
        wi <- which(colnames(bchain) %in%  xi)
      }
      inter <- inter + bchain[,xi]*xx[ox[kk]]
    }
  }
  if( length(i2) > 0 ){
    
    ww <- match(unlist( strsplit(xmnames[i2],m2) ),xmnames)
    ox <- xmnames[ww[is.finite(ww)]]
    for(kk in 1:length(i2)){
      xi    <- paste(xmnames[i2[kk]],specs,sep='_')
      wi <- which(colnames(bchain) %in%  xi)
      if(length(wi) == 0){
        xi    <- paste(specs,xmnames[i2[kk]],sep='_')
        wi <- which(colnames(bchain) %in%  xi)
      }
      inter <- inter + bchain[,xi]*xx[ox[kk]]
    }
  }
  list(main = maine, inter = inter)
}


.stackedBoxPlot <- function( stackList, stackSd=character(0),
                             ylim=NULL,sortBy = NULL, barnames=NULL,
                             col=rep(NULL,length(stackList)),
                             border=rep(NA,length(stackList)),
                             decreasing=T, nsd=1.96, cex=1,
                             legend=NULL, scaleLegend=.1){
  
  # sortBy - if length 1 indicates which variable in stackList to sort by
  #        - if a vector it is the order to plot
  # nds    - no. standard deviations for whiskers
  
  nn  <- length(stackList)
  ord <- c(1:length(stackList[[1]]))
  nx  <- length(ord)
  
  xx <- 0:(nx-1)
  
  if(is.null(ylim)){
    
    ymax <- ymin <- 0
    
    for(j in 1:nn){
      ymax <- ymax + max( c(0,stackList[[j]]),na.rm=T )
      ymin <- ymin + min( c(0,stackList[[j]]),na.rm=T )
    }
    
    ylim <- c(ymin,ymax)
    
    yscale <- diff(ylim,na.rm=T)*.4
    ylim[1] <- ylim[1] - yscale
    ylim[2] <- ylim[2] + yscale
  }
  
  if(!is.null(sortBy)){
    
    if(length(sortBy) > 1){
      ord <- sortBy
    } else {
      ord <- order( stackList[[sortBy]], decreasing = decreasing)
    }
    if(!is.null(barnames))barnames <- barnames[ord]
  }
  
  dy   <- diff(ylim)
  xlim <- c(0,1.2*length(ord))
  
  
  add <- F
  
  offset <- offsetPos <- offsetNeg <- rep(0,length(stackList[[1]]))
  
  if(is.null(col))col <- c(1:nn)
  
  for(j in 1:nn){
    
    xj <- stackList[[j]][ord]
    names(xj) <- NULL
    
    wp <- which(xj > 0)     # increase
    wn <- which(xj < 0)     # decrease
    
    offset[wp] <- offsetPos[wp]
    offset[wn] <- offsetNeg[wn]
    
    hj <- xj 
    
    barplot(height= hj,offset=offset,xlim=xlim,ylim=ylim,
            col=col[j],border=border[j],add=add)
    
    ww <- grep(names(stackList)[j],names(stackSd))
    if(length(ww) > 0){
      xz <- xx + .5
      xz <- xz*1.2
      
      tall <-  nsd*stackSd[[ww]]
      y1   <-  hj + offset + tall
      y2   <-  hj + offset - tall
      
      for(i in 1:length(ord)){
        lines(xz[c(i,i)],c(y1[i],y2[i]),lwd=6,col='white')
        lines(c(xz[i]-.1,xz[i]+.1),y1[c(i,i)],lwd=6,col='white')
        lines(c(xz[i]-.1,xz[i]+.1),y2[c(i,i)],lwd=6,col='white')
        
        lines(xz[c(i,i)],c(y1[i],y2[i]),lwd=2,col=col[j])
        lines(c(xz[i]-.1,xz[i]+.1),y1[c(i,i)],lwd=2,col=col[j])
        lines(c(xz[i]-.1,xz[i]+.1),y2[c(i,i)],lwd=2,col=col[j])
      }
    }
    
    if(j == 1)add <- T
    
    offsetPos[wp] <- offsetPos[wp] + hj[wp]
    offsetNeg[wn] <- offsetNeg[wn] + hj[wn]
    
    if(j == nn & !is.null(barnames)){
      
      xall <- par('usr')[1:2]
      xtic <- c(1:nx)*(diff(xall) - 1)/nx - .8
      
      yy <- offsetPos + .2*dy
      pos <- yy*0 + 1
      wl <- which(abs(offsetNeg) < abs(offsetPos))
      yy[wl] <- offsetNeg[wl] - .2*dy
      pos[wl] <- 4
      text(xtic,yy,barnames,srt=90.,pos=pos,cex=cex)
    }
  } 
  
  if(!is.null(legend)){
    
    dy <- diff(ylim)*scaleLegend
    dx <- 1.2
    x1 <- length(ord)*.02 + 1
    y1 <- ylim[1]
    pos <- 4
    if(legend == 'topright'){
      x1  <- length(ord)
      y1  <- ylim[2]
      dy  <- -dy
      dx <- -dx
      pos <- 2
    }
    if(legend == 'topleft'){
      y1  <- ylim[2]
      dy  <- -dy
    }
    if(legend == 'bottomright'){
      x1  <- length(ord)
      dx <- -dx
      pos <- 2
    }
    for(j in 1:length(stackList)){
      y2 <- y1 + dy
      rect(x1,y1,x1 + 1,y2,col=col[j],border=border[j])
      y1 <- y2
      colj <- col[j]
      if(colj == 'white')colj <- border[j]
      text(x1 + dx,y1 - dy/2,names(stackList)[[j]],col=colj,pos=pos,cex=cex)
    }
  }
  
  invisible( ord )
}  
############################## end gjam functions


regressSegment <- function(xx,yy,nn=nrow(xx),whichx=c(1:nn),whichy = c(1:nn),
                           whichPredx=c(1:nn),whichPredy=c(1:nn),fixBeta = NULL,
                           whichFactor=NULL){
  
  # whichFit  - samples to include in fitting
  # whichPred - samples to include in prediction
  # fixBeta   - any beta values that should be fixed, one row per value, 
  #             rownames(fixBeta) match colnames(xx)
  # whichFactor - names of factors
  
  if(!is.matrix(yy))yy <- matrix(yy,ncol=1)
  
  FIX <- F
  if(!is.null(fixBeta))FIX <- T
  
  fixb <- 0
  ym   <- yy[whichy,]
  Q    <- ncol(xx)
  mb   <- c(1:Q)
  bj   <- matrix(mb*0,ncol=1)
  if(FIX){
    wb   <- match(rownames(fixBeta),colnames(xx))
    fixb <- xx[whichx,wb]%*%fixBeta
    ym   <- ym - fixb
    mb   <- mb[-wb]
  }
  
  xm <- xx[whichx,]
  rx <- qr(xm[,mb])$rank
  if(rx < length(mb)){
    xc  <- cor(xm[,mb])[,2]
    wna <- which( is.na(xc) )
    wna <- wna[-1]
    if(length(wna) > 0)mb  <- mb[-wna]
    
    if(qr(xm[,mb])$rank < length(mb) & !is.null(whichFactor)){
      wf <- match(whichFactor,colnames(xm[,mb]))
      if(length(wf) > 0)mb  <- mb[-wf]
    }
  }
  
  bj[mb,]  <- solve(crossprod(xm[,mb]))%*%crossprod(xm[,mb],ym)
  py       <- xx[whichPredx,]%*%bj 
  
  if(FIX){
    fixb  <- xx[whichPredx,wb]%*%fixBeta
    py    <- py + fixb
  }
  yy[whichPredy,] <- py
  yy
}


removeRows <- function(ww,...){  # ... are names of matrices
	xx <- list(...)
	nx <- length(xx)
	out <- numeric(0)
	
	for(k in 1:nx)out <- append(out, list(get(xx[[k]])[-ww,]))
	names(out) <- xx
	out
}

strsplit2Mat <- function(xx,value){
  
  # value occurs the same number of times in each element of
  #   character vector xx
  
  npart <- length( unlist( strsplit(xx[1],value) ) )
  matrix( unlist( strsplit(xx,value) ),ncol=npart,byrow=T)
}

outFile <- function(outfolder=character(0),file){
  
  if(length(outfolder) > 0){
    W <- file.exists(outfolder)
    if(!W)dir.create(outfolder)
  }
  ov <- paste(outfolder,file,sep='/')
  ww <- grep('//',ov)
  if(length(ww) > 0)ov <- .replaceString(ov,'//','/')
  ov
}

minMatrix <- function(mat,minNumber=5,ncol=NULL,OTHER=F){  
  
  # minNumber - minimum no. of non-zero values in a column of mat
  # ncol      - number of columns to retain, those with highest values
  # if(OTHER) sum of not selected is 'other' column
  
  nc     <- ncol(mat)
  mnames <- colnames(mat)
  
  io <- mat
  io[io > 0] <- 1
  
  csum <- colSums(io)
  ww <- which(csum > minNumber)
  
  ww <-  ww[ order(csum[ww],decreasing=T) ]
  
  if(!is.null(ncol))ww <-ww[1:ncol]
  
  out <- mat[,ww]
  mnames <- mnames[ww]
  
  if(length(ww) < nc & OTHER){
    other <- rowSums(mat[,-ww],na.rm=T)
    out <- cbind(out,other)
    mnames <- c(mnames,'other')
  }
  
  if(!is.matrix(out)){
    out <- matrix(out,ncol=1)
    colnames(out) <- mnames
  }
  csum <- csum[ww]
  
  list(ww = ww, mat = out, csum = csum)
}

plotLines <- function(x,ymat,vline=NULL,ylim=range(ymat),
                      lwd=rep(2,ncol(ymat)),lty=rep(1,ncol(ymat))){
  # ymat  - each column is plotted against x
  # vline - x locations for vertical lines
  
  plot(x,ymat[,1],type='l',ylim=ylim)
  for(j in 1:ncol(ymat))lines(x,ymat[,j],lwd=lwd[j],lty=lty[j])
  
  if(!is.null(vline)){
    for(j in 1:length(vline))abline(v=vline[j],lty=lty[j])
  }
}



imageMap <- function( lonLat, z, colRamp = NULL, 
                      xlim = range(lonLat[,1]), ylim = range(lonLat[,2]),
                      zlim = range( z, na.rm = T), ncol = 20, xaxt = 's', yaxt = 's', 
                      xlab = '', ylab = '', add = F ){
  
  # lonLat has rownames Elon_Nlat, e.g., 'E30.5_N-22.067797'
  # z is a climate vector with names(z) == rownames(lonLat)
  
  if( is.null( colRamp ) )
    colRamp <- c('#8c510a','#bf812d','#dfc27d','#80cdc1','#35978f','#01665e')
  
  cfun  <- colorRampPalette( colRamp )
  cfn   <- cfun(ncol)
  
  ww <- which( lonLat[,1] >= xlim[1] & lonLat[,1] <= xlim[2] &
                 lonLat[,2] >= ylim[1] & lonLat[,2] <= ylim[2] )
  lonLat <- lonLat[ww,]
  z <- z[ww]
  
  x <- sort( unique( lonLat[,1] ) )
  y <- sort( unique( lonLat[,2] ) )
  names(x) <- paste( 'E', x, sep = '' )
  names(y) <- paste( 'N', y, sep = '' )
  zmat <- matrix( NA, length(x), length(y), dimnames = list( names(x), names(y) ) )
  zmat[ columnSplit( names(z), '_' ) ] <- z
  
  image( x, y, zmat, xlim = xlim, ylim = ylim, zlim = zlim, xlab = xlab, ylab = ylab,
         asp = 1, col = colRamp, xaxt = xaxt, yaxt = yaxt )
}

polyMap <- function( polyList, type, lwd = 1, xlim = NULL, ylim = NULL,
                     border = terrain.colors(length(unique(type)), alpha = .5 ), 
                     col = NULL, xaxt = 's', yaxt = 's', 
                     xlab = '', ylab = '', add = F ){
  
  # polyList is a list of polygons with columns x, y
  # type is used for the color legend
  # border and col have length equal to the number of unique values in type
  
  subs <- unique( type )
  if( length( border ) == 1 )border <- rep( border, length( subs ) )
  
  names( border ) <- subs
  cc <- fill <- NA
  if( !is.null( col ) ){
    fill <- col
    names( fill ) <- subs
  }
  
  if( !add ){
    if( is.null( xlim ) ){
      sapply( polyList, range )
      rr <- sapply( polyList, function(x) apply( x, 2, range ) )
      xlim <- range( rr[1:2,] )
      ylim <- range( rr[3:4,] )
    }
    plot( NA, xlim = xlim, ylim = ylim, asp = 1, bty = 'n',
          xlab = xlab, ylab = ylab, xaxt = xaxt, yaxt = yaxt )
  }
  
  for( k in 1:length( polyList ) ){
    kbord <- border[ type[k] ]
    if( !is.null( col ) )cc <- fill[ type[k] ]
    
    polygon( polyList[[k]][,1], polyList[[k]][,2], border = kbord, col = cc,
             lwd = lwd )
  }
}



rampColor <- function( values, 
                       ramp = c('#d53e4f','#fc8d59','#e6f598','#99d594','#3288bd'),
                       nlevs = 10 ){
  # assign colors to the magnitude of values
  
  s  <- seq( min( values, na.rm = T ), max( values, na.rm = T ), length = nlevs )
  si <- findInterval( values, s, all.inside = T )
  cfun <- colorRampPalette( ramp )
  cfun( nlevs )[ si ]
}

colorSegment <- function(xx, yy, colors, nc = 20, lwd=1, lty=1,
                         clear=T){
  
  # add line to plot, colored by magnitude
  # colors - colors to interpolate, see colorRampPalette
  #          interpolate for range(y)
  # x, y   - vectors to plot
  # clear  - clear background of line
  
  yr    <- range(yy)
  sc    <- seq(yr[1],yr[2],length=nc)
  yc    <- c(1:nc)
  y1    <- 1:(nc-1)
  y2    <- y1 + 1
  sc    <- round( yr[1] + yc/nc*diff(yr) ,0)
  
  cf   <- colorRampPalette(colors)(nc)
  ycol <- findInterval(yy,sc,all.inside=T)
  
  if(clear) lines( new$density, new$mids, lwd=lwd*3, col='white')
  segments(xx[y1], yy[y1], xx[y2], yy[y2], col=cf[ycol],lwd=lwd,lty=lty)
}


makeAdvectionMat <- function(v,delta,n){
  
  # v     - velocity
  # delta - discretization
  # n     - dimension
  
  A <- diag(n)
  A[row(A) == col(A) - 1] <- -v*delta    # subdiagonals
  A[row(A) == col(A) + 1] <-  v*delta
  
  A[n,n]   <- 1 - 2*delta*v              #lower boundary
  A[n,n-1] <- 2*delta*v  
  A
}

gaussianRatio <- function(q1,q2,mu,variance){ 
  #returns log ratio of 2 gaussians (e.g., metropolis)
  
  -.5*( (q1^2 - q2^2 -2*mu*(-q1 + q2))/variance)
}

binomialRatio <- function(nn,yy,t1,t2){
  
  yy*(log(t1) - log(t2)) - (nn - yy)*(log(1 - t1) - log(1 - t2))
}

poissonRatio <- function(y,lambda1,lambda2){
  
  y*(log(lambda1) - log(lambda2)) - lambda1 + lambda2
  
}
  
comp_outindex <- function(index_in, total) return((1:total)[-index_in]-1)

bUpdateMVN_Rcpp <- function(xx,yy,bb = NULL,lo = NULL, hi = NULL,sigma,times=1){
  
  nc <- 1
  if(is.matrix(yy))nc <- ncol(yy)
  if(!is.null(bb))dm <- dimnames(bb)
  if(is.null(bb)) dm <- list(colnames(xx),colnames(yy))
  
  testv <- try(chol(crossprod(xx)),T)
  if(inherits(testv,'try-error')){
    message('X not full rank')
    return( bb ) 
  }
  
  cx   <- chol2inv(testv)
  mu   <- matrix(cx %*% crossprod(xx,yy),1)
  vv   <- kronecker(sigma,cx)
  
  if(is.null(lo)){
    b <- matrix( rmvnorm(1,mu,vv),ncol=nc)
    dimnames(b) <- dm
    return(b)
  }
  
  lo <- as.vector(lo)
  hi <- as.vector(hi)

  testvv <- try(chol(vv),T)
  if(inherits(testvv,'try-error')){
    message('kronecker singular')
    return( bb ) 
  }

  tmp <- tnorm.mvtRcpp(as.vector(bb),mu,vv,lo,hi,times=times)
  bb <- matrix(tmp[times,],ncol=nc)
  dimnames(bb) <- dm
  bb
}

makeXY <- function(xdata,ydata,xnames,ynames,factors=NULL,standard=T,complete=F,
                   minLevels=10,combineFactors=NULL){
  
  #xnames columns in xdata for design matrix; should include 'intercept', unless
  #       all factor levels included; if the latter, each is a full intercept
  #factors are names in x taken to be factors
  #interactions have this form: xname1:xname2--identified by ':'
  #ynames columns in ydata for response (univarite will have one value)
  # if 'complete' only complete obs, otherwise NAs are in data;
  #factors cannot be NA
  #minLevels - observations with too few factor levels removed
  
  
  y <- as.matrix( ydata[,ynames] )
  ww <- rowSums(y)
  
  y <- y[ww > 0,]
  xdata <- xdata[ww > 0,]
  
  n <- nrow(xdata)
  x <- numeric(0)
  
  if(!is.null(factors)){
    
    for(j in 1:length(factors)){
      
      xj <- as.character( xdata[,factors[j]] )
      xj[is.na(xj)] <- 'bogus'

      jtab <- table(xj)
      wmin  <- which(jtab < minLevels | names(jtab) == 'bogus')
      wkeep <- which(jtab >= minLevels & names(jtab) != 'bogus')
      
      fnames <- names(jtab)[wkeep]
      other  <- names(jtab)[wmin]
      
      xj[xj %in% other] <- 'other'
      fnames <- c(fnames,'other')
      
      nf     <- length(fnames)
      f      <- matrix(0,n,nf)
      xx <- match(xj,c(fnames))
     
      
      f[ cbind(c(1:n),xx) ] <- 1
      cj <- paste(factors[j],fnames,sep='_')
      colnames(f) <- cj
   #   if('intercept' %in% xnames & j == 1){
        f <- as.matrix(f[,-1],n)
   #     colnames(f) <- cj[-1]
   #   }
      
      x <- cbind(x,f)
    }
    
    otherCols <- grep('other',colnames(x))
    if(!is.null(combineFactors)){
      otherCols <- otherCols[ !otherCols %in% grep(combineFactors,colnames(x)) ]
    }
    otherRows <- unique( which(x[,otherCols] == 1,arr.ind=T)[,1] )
    
    x <- x[-otherRows,]
    x <- x[,-otherCols]
    y <- y[-otherRows,]
    xdata <- xdata[-otherRows,]
    n <- nrow(x)
  }

  
  vnames <- xnames[xnames %in% colnames(xdata) & !xnames %in% factors]
  p      <- length(vnames)
  
  if(p == 1){
    xmean <-  mean(xdata[,vnames],na.rm=T)
    xsd   <-  sd(xdata[,vnames],na.rm=T)
    names(xmean) <- names(xsd) <- vnames
  }
  if(p > 1){
    xmean <- colMeans(xdata[,vnames],na.rm=T)
    xsd   <- apply(xdata[,vnames],2,sd,na.rm=T)
  }
  
  if(p > 0){
    xmat    <- as.matrix(xdata[,vnames])
    
    xcenter <- matrix( xmat - matrix(xmean[vnames],n,p,byrow=T), n,p)
    xstand  <- matrix( xcenter/matrix(xsd[vnames],n,p,byrow=T), n, p)
    colnames(xcenter) <- colnames(xstand) <- vnames
    
    if(!standard)x <- cbind(x,xmat)
    if(standard) x <- cbind(x,xstand)
    
    ww <- grep(':',xnames)
    
    if(length(ww) > 0){
      
      for(j in ww){
        inames <- unlist( strsplit(xnames[j],':') )
        if(!inames[1] %in% factors){
          xi1 <- x[,inames[1]]
          name1 <- inames[1]
        }
        if(inames[1] %in% factors){
          xi1 <- x[,grep(inames[1],colnames(x))]
          name1 <- colnames(x)[grep(inames[1],colnames(x))]
        }
        if(!inames[2] %in% factors){
          xi2 <- x[,inames[2]]
          name2 <- inames[2]
        }
        if(inames[2] %in% factors){
          xi2 <- x[,grep(inames[2],colnames(x))]
          name2 <- colnames(x)[grep(inames[2],colnames(x))]
        }
        xi <- matrix(xi1*xi2,n)
        colnames(xi) <- paste( name1,name2,sep=':')
        x <- cbind(x,xi)
      }
    }
  }
  
  n <- nrow(x)
  
  if('intercept' %in% xnames){
    intercept <- matrix(1,n,1)
    x <- cbind(intercept,x)
    colnames(x)[1] <- 'intercept'
  }
  
  x <- as.matrix(x)
  wmiss <- which(!is.finite(x),arr.ind=T)
  
  xx <- x
  xx[is.na(xx)] <- 0
  rank <- qr(xx)$rank
  
  if(rank < ncol(x))warning('x not full rank')
  
  
  list(x = x, y = y, wmiss = wmiss)
}

makeAR1 <- function(nd,rho,sigma){   # construct AR(1) correlation and covariance matrices of dimension nd
  
  x    <- diag(nd)
  cmat <- rho^abs(row(x) - col(x))
  vmat <- sigma*cmat
  
  list(cormat = cmat, varmat = vmat)
}

invertAR1old <- function(nn,phi,sig){  #dimension, correlation, variance
  
  x1 <- 1/(1 - phi^2)/sig
  x  <- diag( c(x1, rep( (1 + phi^2)/(1 - phi^2)/sig, nn-2), x1) )
  x[ row(x) == col(x) - 1 ] <- x[ row(x) == col(x) + 1 ] <-  -phi/(1 - phi^2)/sig 
  x
}

invertAR1 <- function(nd,phi,sig){  #dimension, correlation, variance
  
  x  <- diag( c(1, rep( (1 + phi^2), nd-2), 1) )
  x[ row(x) == col(x) - 1 ] <- x[ row(x) == col(x) + 1 ] <-  -phi 
  x/sig
}

xy2area <- function(xy){  # convex hull and polygon area

  require(splancs)
  
  xy <-  xy[is.finite(xy[,1]) & is.finite(xy[,2]),] 
  hull  <- as.matrix(xy[chull(xy),])
  area  <- areapl(hull)
  
  list(hull = hull, area = area)
}

pasteVector <- function(vec,sep=''){     #paste a vector to create one variable

  x <- vec[1]
  if(length(vec) > 1)for(k in 2:length(vec))x <- paste(x,vec[k],sep=sep)
  x
}

aggregateSequence <- function(index,dat,action='mean',minObs=NULL){    
  # index - column index for aggregation
  # dat is n by t
  
  nc <- length(index)
  if(!is.matrix(dat)){
    dat <- matrix(dat,ncol=nc)
  }
  nn <- nrow(dat)
	
	allDates <- sort(unique(index))
  nm <- length(allDates)
  
	newDat   <- matrix(NA,nn,nm)
	rownames(newDat) <- rownames(dat)
  colnames(newDat) <- allDates
	
	jk <- 0
	for(j in allDates){
		jk <- jk + 1
		wj <- which(index == j)
		mj <- dat[,wj]
    
    if(!is.null(minObs))if(length( which(is.finite(mj)) ) < minObs)next

    if(action == 'mean'){
		  if(is.matrix(mj)) newDat[,jk] <- rowMeans(mj,na.rm=T)
		  if(!is.matrix(mj))newDat[,jk] <- mean(mj,na.rm=T)
    }
    if(action == 'sum'){
		  if(is.matrix(mj)) newDat[,jk] <- rowSums(mj,na.rm=T)
		  if(!is.matrix(mj))newDat[,jk] <- sum(mj,na.rm=T)
    }
	}
	list(dates = allDates, data = newDat)
}

gapFillSequence <- function(y,baseline,fitwidth,x=NULL,tooLo=-Inf,tooHi=Inf){
	
	#y        - sequence with gaps
	#baseline - reference sequence
	#fitwidth - number of preceding values for calibration
	#x        - additional predictors
	
	nf <- length(y)
	if(length(x) > 0 & !is.matrix(x))x <- matrix(x,ncol=1)
	
	xmat <- matrix(NA,nrow(x),ncol=fitwidth)
		
   ki <- c(fitwidth:nf)
   kj <- ki
		
   for(k in 1:fitwidth){
       xmat[ki,k] <- baseline[kj]
		 kj <- kj - 1
	}
	if(length(x) > 0)xmat <- cbind(xmat,x)
	
	wr   <- which(is.finite(rowSums(xmat)) & is.finite(y))
	wr   <- wr[y[wr] > tooLo]              #bad values
	wr   <- wr[y[wr] < tooHi]
	xx   <- xmat[wr,]
	yy   <- y[wr]

	xcol <- c(1:ncol(xx))
	if(sum(xx[,1]) == sum(xx[,2])){
		xcol <- c(1:ncol(xx))[-1]
	}
	b  <- solve(crossprod(xx[,xcol]))%*%crossprod(xx[,xcol],yy)
		
	predy <- xmat[,xcol]%*%b
		
	wf    <- which(!is.finite(y))
	y[wf] <- predy[wf]
	y
}


getDateSequence <- function(year, dayInterval){    #dates from Julian dates
  
  jd <- seq(dayInterval, 365, by=dayInterval)
  da <- as.Date(jd - 1, origin = paste(year,"-01-1",sep=""))
  list(jd = jd,date = da)
}


daysSinceDate <- function(month0,day0,yr0,mo,da,yr,nineteen=F){  #Julian dates from dates
	
	#vectors for mo, da, yr since initial date 0

	dd <- clark_mdy.date(mo,da,yr,nineteen)
	dd - clark_mdy.date(month0, day0, yr0,nineteen) + 1
	
}


clark_mdy.date <- function (month, day, year, 
                            nineteen = TRUE, fillday = FALSE, 
                            fillmonth = FALSE){
                            	
    #days since 1/1/1960

    temp <- any( (month != trunc(month)) | (day != trunc(day)) | 
        (year != trunc(year)))
    if (!is.na(temp) && temp) {
        warning("Non integer input values were truncated in mdy.date")
        month <- trunc(month)
        day   <- trunc(day)
        year  <- trunc(year)
    }
    if (nineteen) 
      year <- ifelse(year < 100, year + 1900, year)
     temp  <- numeric(length(month + day + year))
     month <- month + temp
     day   <- day + temp
     year  <- year + temp
    if (fillmonth) {
        temp <- is.na(month)
        month[temp] <- 7
        day[temp]   <- 1
    }
    if (fillday) 
        day[is.na(day)] <- 15
        month[month < 1 | month > 12] <- NA
        day[day < 1] <- NA
        year[year == 0] <- NA
    year   <- ifelse(year < 0, year + 1, year)
    tyear  <- ifelse(month > 2, year, year - 1)
    tmon   <- ifelse(month > 2, month + 1, month + 13)
    julian <- trunc(365.25 * tyear) + trunc(30.6001 * tmon) + 
        day - 715940
    temp <- trunc(0.01 * tyear)
    save <- ifelse(julian >= -137774, 
                   julian + 2 + trunc(0.25*temp) - temp, 
                   julian)
    year   <- ifelse(month == 12, year + 1, year)
    month  <- ifelse(month == 12, 1, month + 1)
    day    <- 1
    tyear  <- ifelse(month > 2, year, year - 1)
    tmon   <- ifelse(month > 2, month + 1, month + 13)
    julian <- trunc(365.25 * tyear) + trunc(30.6001 * tmon) + day - 715940
    temp   <- trunc(0.01 * tyear)
    save2  <- ifelse(julian >= -137774, julian + 2 + trunc(0.25*temp) - temp, julian)
    temp   <- as.integer(ifelse(save2 > save, save, NA))
    temp
}

yearMonthVecCens <- function(yrvec){   
	
	#returns days to end of each month for years in yrvec
	
	allYears <- sort(unique(yrvec))
	nyr <- length(allYears)
	
	mNames   <- c('jan','feb','mar','apr','may','jun','jul','aug',
	              'sep','oct','nov','dec')
	mDays <- rep(31,12)
	mDays[c(4,6,9,11)] <- 30
	mDays[2] <- 28
	names(mDays) <- mNames
	endMonth <- cumsum(rep(mDays,times=nyr))
	endYr    <- cumsum(rep(sum(mDays),nyr))
	list(endMonth = endMonth,endYr = endYr)
}

bigObjects <- function(nn=10){  #find big objects

  #find large objects
  zz <- sapply(ls(pos = 1), function(x)
                 object.size( get(x, envir = globalenv())) )
  zz <- as.matrix(rev(sort(zz))[1:nn])
  colnames(zz) <- 'Mb'
  print( round(zz*1e-6, 2), units = 'Mb')
}


logit2Prob <- function(z){   #multivar logit to fractions
  
  # z has 1 fewer columns than probabilities
	
  if(!is.matrix(z))z <- matrix(z,nrow=1)
  
  S    <- ncol(z)
  ze   <- exp(z[,1:(S-1)])
  zs   <- rowSums( ze )
  zS   <- 1/(1 + zs)
  zz   <- ze*zS
  
  cbind(zz,zS)

}

prob2Logit <- function(y){     #fractions to multivar logit      
  
  S <- ncol(y)
  
  pS <- log(y) - log(y[,S])
  pS[,-S]
  
}

################################
stateSpaceMVNlogit <- function(){
	#z is the observed proportions of plot/time by species
	
#	krons <- kronecker(diag(1,2),sg)

   omat <- matrix(obsError,n,ns,byrow=T)
   accept <- 0
	
	for(t in 1:nt){
				
		ti <- t + tindex
		propy <- matrix(rnorm(n*ns,yg[ti,],.3),n,ns)
		
		pnow <- rowSums(dnorm(yg[ti,],yobs[ti,],sqrt(obsError),log=T))
		pnew <- rowSums(dnorm(propy,yobs[ti,],sqrt(obsError),log=T))
		
		if(t > 1){
			mu <- xtime[ti-1,]%*%bg + yg[ti-1,]%*%ag
			pnow <- pnow + diag(-(yg[ti,] - mu)%*%sinv%*%t(yg[ti,] - mu)/2)
			pnew <- pnew + diag(-(propy - mu)%*%sinv%*%t(propy - mu)/2)
		}	
		if(t < nt){
			mu1 <- xtime[ti,]%*%bg + yg[ti,]%*%ag
			mu2 <- xtime[ti,]%*%bg + propy%*%ag
			pnow <- pnow + diag(-(yg[ti+1,] - mu1)%*%sinv%*%t(yg[ti+1,] - mu1)/2)
			pnew <- pnew + diag(-(yg[ti+1,] - mu2)%*%sinv%*%t(yg[ti+1,] - mu2)/2)
		}	
		a <- exp(pnew - pnow)
		za <- runif(n,0,1)
		wp <- which(za < a)
		if(length(wp) > 0){
			yg[ti[wp],] <- propy[wp,]
			accept <- accept + length(wp)
		}
   }
   list(yg = yg, accept = accept/nrow(yg))
}

obsErrorMVNlogit <- function(){
	
	u1 <- s1 + .5*n*ns*nt
	u2 <- s2 + .5*sum( (yg - yobs)^2 )
	1/rgamma(1,u1,u2)
}		

vec2Mat <- function(ir,ic,x,nr=max(ir),nc=max(ic)){   #vector x to matrix x[ir,ic]

  xmat <- matrix(NA,nr,nc)
  xmat[cbind(ir,ic)] <- x
  xmat

}

mat2Vec <- function(x,GETINDEX=F,lastTime=NULL){  #matrix x to vector, after lastTime removed

  xx  <- t(x)
  nr <- nrow(xx)
  nc <- ncol(xx)

  rowCol <- index <- numeric(0)

  if(is.null(colnames(xx)))colnames(xx) <- 1:nc
  if(is.null(rownames(xx)))rownames(xx) <- 1:nr

  
  ir <- as.vector( matrix(c(1:nr),nr,nc) )
  ic <- as.vector( matrix(c(1:nc),nr,nc,byrow=T) )
  y <- xx[cbind(ir,ic)]

  if(GETINDEX){
    index  <- expand.grid(col = c(1:nr),row = c(1:nc) )[,c(2,1)]

    if(!is.null(lastTime)){

      wl <- which(is.finite(lastTime))
      if(length(wl) > 0){
         wk <- numeric(0)
         for(k in wl)wk <- c(wk, which(index[,1] == k & index[,2] > lastTime[k]) )

         if(length(wk) > 0){
            y <- y[-wk]
            index <- index[-wk,]
         }
       }
    }
  }

  list(x = y, index = index)

}

interp <- function(y,INCREASING=F,minVal=-Inf,maxVal=Inf,defaultValue=NULL,
                   tinySlope=NULL){  #interpolate vector x

  if(is.null(defaultValue))defaultValue <- NA

  tiny <- .00001
  if(!is.null(tinySlope))tiny <- tinySlope

  y[y < minVal] <- minVal
  y[y > maxVal] <- maxVal

  n  <- length(y)
  wi <- which(is.finite(y))

  if(length(wi) == 0)return(rep(defaultValue,n))
  if(length(wi) == 1)ss <- tiny

  xx  <- c(1:n)
  z  <- y

  if(wi[1] != 1) wi <- c(1,wi)
  if(max(wi) < n)wi <- c(wi,n)

  ss <- diff(z[wi])/diff(xx[wi])

  ss[is.na(ss)] <- 0

  if(length(ss) > 1){
    if(length(ss) > 2)ss[1] <- ss[2]
    ss[length(ss)] <- ss[length(ss)-1]
  }
  if(INCREASING)ss[ss < tiny] <- tiny

  if(is.na(y[1]))  z[1] <- z[wi[2]] - xx[wi[2]]*ss[1]
  if(z[1] < minVal)z[1] <- minVal
  if(z[1] > maxVal)z[1] <- maxVal

  for(k in 2:length(wi)){

     ki <- c(wi[k-1]:wi[k])
     yk <- z[wi[k-1]] + (xx[ki] - xx[wi[k-1]])*ss[k-1]
     yk[yk < minVal] <- minVal
     yk[yk > maxVal] <- maxVal
     z[ki] <- yk
  }
  z
}
  
interpRows <- function(x,startIndex=rep(1,nrow(x)),endIndex=rep(ncol(x),nrow(x)),
                       INCREASING=F,minVal=-Inf,maxVal=Inf,
                       defaultValue=NULL,tinySlope=.001){  
  #interpolate rows of x subject to increasing

  nn  <- nrow(x)
  p  <- ncol(x)
  xx <- c(1:p)

  if(length(minVal) == 1)minVal <- rep(minVal,nn)
  if(length(maxVal) == 1)maxVal <- rep(maxVal,nn)

  ni   <- rep(NA,nn)
  flag <- numeric(0)

  z <- x

  for(i in 1:nn){
    if(startIndex[i] == endIndex[i]){
      z[i,-startIndex[i]] <- NA
      next
    }
    z[i,startIndex[i]:endIndex[i]] <- interp(x[i,startIndex[i]:endIndex[i]],
                                             INCREASING,minVal[i],maxVal[i],
                                             defaultValue,tinySlope)
  }
  
  z
}

byFunction <- function(x,i,j,mat=matrix(0,max(i),max(j)),FUN=sum){  
  
  # for 2-D, ... can be i,j
  # mat is max(i) by max(j)
  
  FUN <- match.fun(FUN)
  
  g <- interaction(i,j)
  split(x, g) <- lapply( split(x,g), FUN, na.rm=T )
  mat[cbind(i,j)] <- x
  mat
  
}

byIndex <- function(xx,INDICES,FUN,coerce=F,...){  
  
#INDICES is list, each same length as  x
  
#  fun <- match.fun(FUN)
  
  if(!is.list(INDICES))INDICES <- list(INDICES)
  
  nl <- length(INDICES)
  
  tmp  <-  unlist(by( as.vector(xx),INDICES,FUN,...) ) 
  nd   <- dim(tmp)
  nm   <- names(tmp)
  tmp  <- array(tmp,dim=nd, dimnames=dimnames(tmp))
  tmp2 <- matrix(tmp,ncol=1)
  rownames(tmp2) <- nm
  
  tmp2[is.na(tmp2)] <- 0
  
  if(!coerce)return(tmp2)
  
  if(nl != 2)return(tmp)
  
  #following only works for 2 dimensions:
  
  drange <- numeric(0)
  for(j in 1:nl)drange <- rbind(drange,range(INDICES[[j]]))
  
  dname <- dimnames(tmp)
  mk    <- rep(0,length(nd))
  tnew  <- array(0,dim=drange[,2])
  id    <- numeric(0)
  
  for(k in 1:length(nl)){
    cc <- drange[k,1]:drange[k,2]
    id <- append(id,list(match(dimnames(tmp)[[k]],cc)))
  }
  
  tnew[id[[1]],id[[2]]] <- tmp
  rownames(tnew) <- 1:drange[1,2]
  colnames(tnew) <- 1:drange[2,2]
    
  tnew
}

sumByIndex <- function(x,pvec,imax){  

  #x - vector to sum, vector of indices,max index
  #    if x is a matrix, then one value of pvec per row of x
  #returns sum by index, inclusive from 1:imax
  
  if(is.matrix(x))x <- rowSums(x,na.rm=T)

     svec <- as.matrix( by(x,pvec,sum,na.rm=T) )

     psum <- rep(0,imax)
     psum[ match(rownames(svec),c(1:imax)) ] <- svec
     names(psum) <- c(1:imax)
     psum
}

sumBy2Index <- function(x,pvec1,pvec2,
                        imax1=max(pvec1),imax2=max(pvec2)){   #

  smat <- matrix(NA,imax1,imax2)

  svec <- by( x, list(i = pvec1, j =  pvec2), sum, na.rm=T) 
  svec <- as( unlist(svec),'matrix' )

  ir <- match(rownames(svec),c(1:imax1))
  ic <- match(colnames(svec),c(1:imax2))
  ii <- as.matrix( expand.grid(ir,ic) )
  smat[ii] <- as.vector(svec)
  smat
}

inData <- function(filename, xnames = NULL, ynames = NULL, tname = NULL, 
                   iname = NULL, na.rm = F, INTERCEPT = F){  
                   	
  #read in data file, return design matrix x, response y
  #xnames, ynames, tname, iname are column headings in filename
  #time indicator t
  #individual indicator i

  data <- read.table(filename,header=T)
  
  if(is.atomic(xnames)){
    x    <- data[,xnames]
    if(!is.matrix(x))x <- as.matrix(x)
    if(INTERCEPT){
      intercept <- rep(1,nrow(data))
  	   x <- cbind(intercept,x)
  	   colnames(x) <- c('intercept',xnames)
    }
  }
  if(is.atomic(ynames)){
    y <- data[,ynames]
    if(!is.matrix(y))y <- as.matrix(y) 
    y  <- matrix(y,nrow(data),length(ynames))
    colnames(y) <- ynames
  }
  
  wf <- c(1:nrow(data))
  
  if(na.rm){
  	 wf <- which(is.finite(rowSums(x)) & is.finite(rowSums(y)))
  	 x  <- x[wf,]
  	 y  <- y[wf,]
  }
  
  z  <- list(x = x, y = y)
  
  if(is.atomic(tname))z$t <- data[wf,tname]
  if(is.atomic(iname))z$i <- data[wf,iname]
  z
}

treebytime <- function(iindex,tindex,x,nr=max(iindex),nc=max(tindex)){ #matrix for x with individuals by time
	
  y <- matrix(NA,nr,nc)
  y[cbind(iindex,tindex)] <- x
  y
}

pRate <- function(par){

  g <- par[1]
  r <- par[2]
  c <- par[3]
  s <- par[4]
  
  f <- g*(1 - exp(-r*x))
  -sum(dnorm(y,c + f,s*f,log=T))
}

byHour <- function(q){  #for capacitance data (Ward et al. 2013)
	
	hrSeq <- sort(unique(b[,'hour']))
	dySeq <- sort(unique(b[,'JD']))
	nh    <- length(hrSeq)
	nd    <- length(dySeq)
	
	qmat <- matrix(NA,nd,nh)
	colnames(qmat) <- hrSeq
	rownames(qmat) <- dySeq
	
	di <- match(b[,'JD'],dySeq)
	hi <- match(b[,'hour'],hrSeq)
	qmat[cbind(di,hi)] <- q
	list(q = qmat, days = dySeq, hours = hrSeq)
	
}

############################### map species

mapSpecies <- function(x,y,z,mapx=range(x),mapy=range(y),
                       scale=0,add=F,sym='circles',colVec=rep(1,length(x)),fill=F){
	
   fillCol <- NA
   if(fill)fillCol <- colVec
   if(scale > 0)mapSetup(mapx,mapy,scale)
   if(sym == 'circles')symbols(x,y,circles=z/10,inches=F,xlim=mapx,ylim=mapy,
                               fg=colVec,bg=fillCol,lwd=2,add=add)
   if(sym == 'squares')symbols(x,y,squares=z/10,inches=F,xlim=mapx,ylim=mapy,
                               fg=colVec,bg=fillCol,lwd=2,add=add)
}

myECDF <- function(x){
  n  <- length(x)
  xx <- sort(unique(x))                 
  fx <- cumsum(tabulate(match(x,xx))/n) 
  list(x = xx, fx = fx)                 
}

samplePlots <- function(mapx,mapy,wide,nplot,mapscale=1,PLOTIT = T){

  yt      <- seq(mapy[1],mapy[2],by=wide)      #y grid locations
  yl      <- length(yt)
  xt      <- seq(mapx[1],mapx[2],by=wide)      #x grid locations
  xl      <- length(xt)
  mapgrid <- cbind(rep(xt,each=yl),rep(yt,xl))       #x and y locations
  loc     <- mapgrid[sample(yl*xl,nplot,replace=F),] #samples from grid

  sl <- .5*wide                             #plot edges
  xbound <- cbind((loc[,1]-sl),(loc[,1]+sl))
  ybound <- cbind((loc[,2]-sl),(loc[,2]+sl))
  xindex <- c(1,2,2,1,1)
  yindex <- c(1,1,2,2,1)
  
  specname  <- sort(unique(treedata[,'species']))
  nspec     <- length(specname)

 # plot.data  <- numeric(0)                   #list of observations
  tableSpec <- matrix(0,nspec,nplot)
  rownames(tableSpec) <- specname

  if(PLOTIT)plot(-1000,0,xlim=mapx,ylim=mapy)

  for(i in 1:nplot){

  # extract trees on sample plot i and store them in table.spec
    xt <- !is.na(cut(treedata[,'x'],breaks=xbound[i,],exclude=NA))
    yt <- !is.na(cut(treedata[,'y'],breaks=ybound[i,],exclude=NA))
    
    tmp <- treedata[xt & yt,]
    if(nrow(tmp) > 0){
      tableSpec[,i] <- table(tmp[,'species'])
      if(PLOTIT)symbols(tmp[,'x'],tmp[,'y'],circles=tmp[,'dbh']/20,inches=F,add=T)
    }

  # draw a box around each plot
    if(PLOTIT){
      xvec <- xbound[i,xindex]
      yvec <- ybound[i,yindex]
      lines(xvec,yvec)
    }
  }

  tableSpec
}

appendData <- function(oldfile,newfile,oldDates,newDates){  #append new plot data
    	
    rold <- rownames(oldfile)
    rnew <- rownames(newfile)

    	if(!is.matrix(newfile)){
    		newfile <- matrix(newfile,length(newfile),1)
    		colnames(newfile) <- newDates
    	}

    	if(length(oldfile) == 0)return(newfile)
    	
    	if(!is.matrix(oldfile)){
    		oldfile <- matrix(oldfile,length(oldfile),1)
    		colnames(oldfile) <- oldDates
    	}
    	
      allDates <- unique(c(oldDates,newDates))
    	
      newMat  <- matrix(NA,nrow(newfile),length(allDates))
      colnames(newMat) <- allDates

      wwc     <- match(colnames(newfile), allDates)
      newMat[,wwc] <- newfile
      newMat  <- matrix(newMat,nrow(newfile),length(allDates))
      colnames(newMat) <- allDates
      rownames(newMat) <- rnew

      oldMat <- matrix(NA,nrow(oldfile),length(allDates))
      colnames(oldMat) <- allDates
      wwc    <- match(colnames(oldfile), allDates)
      oldMat[,wwc] <- oldfile
      oldMat <- matrix(oldMat,nrow(oldfile),length(allDates))
      colnames(oldMat) <- allDates
      rownames(oldMat) <- rold
      
      rbind(oldMat,newMat)
}

missingCol <- function(x,colName,action='warn',value=NA){   

# check for missing columns in vector colName
# action = 'warn','stop','add'

  nx <- nrow(x)
  
  newCols <- character(0)
  
  checkx <- which( duplicated(colnames(x)) )
  checkc <- which( duplicated(colName ) )
  
  if(length(checkc) > 0){
    warning('duplicated columns in species list')
  }

  for(j in 1:length(colName)){

    if(colName[j] %in% colnames(x))next

    if(action == 'warn')warning( paste('missing',colName[j]) )
    if(action == 'stop')   stop( paste('missing',colName[j]) )

    if(action == 'add'){

      if(length(value) == 1) new <- rep(value,nx)
      if(length(value) == nx)new <- value
      if(!length(value) %in% c(1,nx))stop('value must be 1 or nrow(x)')
     
      x   <- cbind(x,new)
      colnames(x)[ncol(x)] <- colName[j]
      newCols <- c(newCols,colName[j])
    }
  }
#  x <- x[,colName]
  invisible(x)
}

hist2Add <- function(xx,xlim=NULL,ylim=NULL,title='',labels=F,cex=1){
  
  # add histogram to existing plot using add.scatter, library ade4
  
  if(is.null(xlim))xlim <- quantile(xx,c(.001,.999),na.rm=T)
  
  xx[xx > xlim[2]] <- xlim[2]
  opar <- par('mar','yaxt','plt')
  on.exit(par(opar))
  par(mar = rep(.1,4),yaxt='n',plt=par('plt'))
  hist(xx,xlab='',ylab='',main=' ',col='brown',probability=T,nclass=60,
       xlim=xlim,ylim=ylim,labels=labels,cex=cex)
}

KLdivergenceMVN <- function(mu0=0,mu1=0,s0,s1){ 
  
  m <- nrow(s0)
  if(length(mu0) == 1)mu0 <- rep(mu0,m)
  if(length(mu1) == 1)mu1 <- rep(mu1,m)
  
  si <- solve(s1)
  .5*( sum(diag(si%*%s1)) + t(mu1 - mu0)%*%si%*%(mu1 - mu0) - m - + log(det(s1)/det(s0)) )
}
  
KLdivergence <- function(probs,q,NORM=T){  # discrete K-L divergence
  
  # probs, q are normalized
  # if matrices, rows are distributions
  # if !NORM, then must be normalized
  
  if(!is.matrix(probs)){
    probs <- matrix(probs,1)
    q     <- matrix(q,1)
  }
  if(!NORM){
    p1 <- rowSums(probs)
    q1 <- rowSums(q)
    w0 <- which(p1 == 0)
    probs <- probs/matrix(p1,nrow(probs),ncol(probs))
    q     <- q/matrix(q1,nrow(probs),ncol(probs))
  }
  
  rat <- probs/q
  prr <- log(rat)
  prr[rat == 0 | !is.finite(prr)] <- 0
  
  list(KL = rowSums(prr*probs), zero = w0)
}

appendMatrix <- function(m1,m2,fill=NA,SORT=F,asNumbers=F){  
  
  # matches matrices by column names
  # asNumbers: if column heads are numbers and SORT, then sort numerically
  
  if( is.list(m2) )m2 <- as.matrix(m2)

   if( length(m1) == 0 ){
     if( is.matrix(m2) ){
       m3 <- m2
     } else {
       m3 <- as.matrix(m2,nrow=1)
     }
     if( !is.null(names(m2)) )colnames(m3) <- names(m2)
     return(m3)
   }
   if(length(m2) == 0){
     if(!is.matrix(m1))m1 <- matrix(m1,nrow=1)
     return(m1)
   }
   if( is.vector(m1) | (length(m1) > 0 & !is.matrix(m1)) ){
     nn <- names(m1)
     if(is.null(nn))warning('cannot append matrix without names')
     m1 <- matrix(m1,1)
     colnames(m1) <- nn
   }  
   if( is.vector(m2) | (length(m2) > 0 & !is.matrix(m2)) ){
     nn <- names(m2)
     if(is.null(nn))warning('cannot append matrix without names')
     m2 <- matrix(m2,1)
     colnames(m2) <- nn
   }

   c1 <- colnames(m1)
   c2 <- colnames(m2)
   r1 <- rownames(m1)
   r2 <- rownames(m2)
   n1 <- nrow(m1)
   n2 <- nrow(m2)

   allc <-  unique( c(c1,c2) ) 
   if(SORT & !asNumbers)allc <- sort(allc)
   if(SORT & asNumbers){
     ac <- as.numeric(allc)
     allc <- as.character( sort(ac) )
   }

   nr <- n1 + n2
   nc <- length(allc)

   if(is.null(r1))r1 <- paste('r',c(1:n1),sep='-')
   if(is.null(r2))r2 <- paste('r',c((n1+1):nr),sep='-')
   new <- c(r1,r2)

   mat1 <- match(c1,allc)
   mat2 <- match(c2,allc)

   out <- matrix(fill,nr,nc)
   colnames(out) <- allc
   rownames(out) <- new

   out[1:n1,mat1] <- m1
   out[(n1+1):nr,mat2] <- m2
   out
}
   
row2Mat <- function(vec){

  if(is.matrix(vec))return(vec)
  vn  <- names(vec)
  vec <- matrix(vec,1)
  colnames(vec) <- vn
  vec
}
############################################
col2Mat <- function(vec,namecol=NULL){

  if(is.matrix(vec))return(vec)
  vn  <- names(vec)
  vec <- matrix(vec,ncol=1)
  rownames(vec) <- vn
  colnames(vec) <- namecol
  vec
}
##############################################
rowBind <- function(matnow,row2add,rowName){

  if(length(matnow) == 0){
    matnow <- row2add
    if(!is.matrix(row2add)){
        matnow <- matrix(matnow,nrow=1)
        if(!is.null(names(row2add)))colnames(matnow) <- names(row2add)
    }
    rownames(matnow) <- rowName
    return(matnow)
  }

  if(!is.matrix(row2add))row2add <- t(as.matrix(row2add))

  matnow <- rbind(matnow,row2add)
  lastrows <- c( (nrow(matnow) - nrow(row2add) + 1):nrow(matnow) )
  rownames(matnow)[lastrows] <- rowName
  matnow
}

variancePrior <- function(mu,wt,maxFactor=100){

  s1 <- wt
  s2 <- mu*(s1 - 1)
  lo <- 1e-10
  hi <- maxFactor*mu

  list(mu = mu, s1 = s1, s2 = s2, lo = lo, hi = hi)
}

#############################################

values2grid <- function(x, y, z, nx=NULL, ny=NULL, dx=NULL, dy=NULL,
                        ksearch = 5, MATFORMAT=T){
  
  xl <- range(x)
  yl <- range(y)
  
  xs <- seq(xl[1],xl[2],length=nx)
  ys <- seq(yl[1],yl[2],length=ny)
  
  grid <- as.matrix(expand.grid(xs,ys))
  
  tmp <- RANN::nn2(cbind(x,y),grid,k=ksearch)
  nn  <- tmp[[1]]
  wt  <- tmp[[2]]
  mn  <- min(wt[wt > 0])
  wt  <- 1/(wt + mn)
  
  zz  <- matrix(z[nn],nrow(nn),ncol(nn))
  zz  <- rowSums(zz*wt)/rowSums(wt)
  
  if(!MATFORMAT)return(  cbind(grid,zz) )
  
  zmat <- matrix(NA, nx, ny)
  ix  <- match(grid[,1],xs)
  iy  <- match(grid[,2],ys)
  zmat[ cbind(ix,iy) ] <- zz
  
  rownames(zmat) <- xs
  colnames(zmat) <- ys
  
  list( x = xs, y = ys, z = zmat )
}



values2contour <- function(xx, yy, z, nx=100,ny=100,lty=1,labcex=.7,
                           xlim = NULL, ylim = NULL,
                           col='black',lwd=1,zlevs=NULL,add=T,fill=F,
                           drawlabels=F){    
  
  # contours where x,y is not a uniform grid, requires 'spatial' library
  
  require(MBA)
  
  xyzmat <- cbind(xx,yy,z)
  
  wna <- apply(xyzmat,1,sum)
  wna <- which(is.na(wna))
  if(length(wna) > 0)xyzmat <- xyzmat[-wna,]
  
  colnames(xyzmat) <- c('x','y','z')
  
  if(!is.null(xlim)){
    wx <- which(xyzmat[,'x'] >= xlim[1] & xyzmat[,'x'] <= xlim[2])
    xyzmat <- xyzmat[wx,]
  }
  if(!is.null(ylim)){
    wx <- which(xyzmat[,'y'] >= ylim[1] & xyzmat[,'y'] <= ylim[2])
    xyzmat <- xyzmat[wx,]
  }
  
  maxz <- max(xyzmat[,'z'])
  if(!is.null(zlevs)){
    zlevs <- zlevs[zlevs <= maxz]
    colM <- colorRampPalette(c('mintcream','azure','green','darkgreen','darkgreen'))
    nz     <- length(zlevs) - 1
    cols <- colM(nz)
  }
  
  surf  <- mba.surf(xyz=xyzmat,no.X=nx,no.Y=ny,h=7,sp=F,extend=F)$xyz.est
  
  # surf  <- myMBA.surf(xyz=xyzmat,no.X=nx,no.Y=ny,h=7,sp=F,extend=F)$xyz.est
  
  if(is.null(zlevs)){
    zlevs <- signif(seq(min(z), max(z), length=3),1)
  }
  contour(surf, levels=zlevs,lwd=lwd,lty=lty,col=col,add=add,labcex=labcex,
          drawlabels=drawlabels)
  
  if(fill){
    zl <- zlevs
    if(length(zl) == 1)stop('fill.contour() needs at least 2 contour lines')
    .filled.contour(surf$x,surf$y,surf$z,levels=zl,col=col)
  }
  invisible(surf)
}

values2contourOld <- function(xx,yy,z,nx=100,ny=100,lty=1,labcex=.7,
                           col='black',lwd=1,zlevs=NULL,add=T,fill=F,
                           drawlabels=F){    

  # contours where x,y is not a uniform grid, requires 'spatial' library

  require(MBA)

  xyzmat <- cbind(xx,yy,z)
  
  wna <- apply(xyzmat,1,sum)
  wna <- which(is.na(wna))
  if(length(wna) > 0)xyzmat <- xyzmat[-wna,]
  
  colnames(xyzmat) <- c('x','y','z')
  
  print(range(xyzmat))
  
  surf  <- mba.surf(xyz=xyzmat,no.X=nx,no.Y=ny,h=7,sp=F,extend=F)$xyz.est
  
 # surf  <- myMBA.surf(xyz=xyzmat,no.X=nx,no.Y=ny,h=7,sp=F,extend=F)$xyz.est
  
  if(is.null(zlevs)){
    zlevs <- signif(seq(min(z), max(z), length=3),1)
  }
  contour(surf, levels=zlevs,lwd=lwd,lty=lty,col=col,add=add,labcex=labcex,
          drawlabels=drawlabels)
  
  if(fill){
    zl <- zlevs
    if(length(zl) == 1)stop('fill.contour() needs at least 2 contour lines')
    .filled.contour(surf$x,surf$y,surf$z,levels=zl,col=col)
  }
  invisible(surf)
}


nearestNeighbors <- function(xy1,xy2,k=1){ #find k nearest neighbots
  
  require(RANN)
  nn2(xy1,xy2,k)$nn.idx
}

#####################################################

crosscorByRow <- function(xmat,ymat,lag=ncol(xmat),BOOTSTRAP=F,PLOT=F){  
  #cross correlation for each row of xmat[i,] vs ymat[i,]

  xmat <- row2Mat(xmat)
  ymat <- row2Mat(ymat)

  nn <- nrow(xmat)
  xx <- c(-lag:lag)
  nc <- length(xx)
  yy <- matrix(NA,nn,nc)
  ii <- numeric(0)

  ciMean <- numeric(0)

  for(i in 1:nn){
    di <- xmat[i,]
    fi <- ymat[i,]
    wi <- which(is.finite(fi) & is.finite(di) & fi > 0)
    if(length(wi) < lag/2)next
    if(var(fi[wi]) == 0)next

  #  xxx <- xx[wi]
    xxx <- wi - mean(wi)
    ddd <- di[wi]
    fff <- fi[wi]
    di  <- lm(ddd ~ xxx)$residuals
    fi  <- lm(fff ~ xxx)$residuals

    cross <- ccf(di,fi,type='correlation',plot=F)
    cx <- cross$lag
    cy <- cross$acf

    cy <- cy[cx %in% xx]
    cx <- cx[cx %in% xx]
    yy[i,match(cx,xx)] <- cy
    ii <- c(ii,i)
  }

  nk <- length(ii)   #sample size for good series
  ci <- matrix(NA,3,ncol(yy))
  rownames(ci) <- c('50%','2.5%','97.5%')

  if(nk == 1)ci[1,] <- yy[ii,]
  if(nk > 1 & nk < 10)ci[1,] <- apply(yy,2,mean,na.rm=T)
  if(nk > 10){
    ci <- apply(yy,2,quantile,c(.5,.025,.975),na.rm=T)
 #   ci[1,] <- apply(yy,2,mean,na.rm=T)
    colnames(ci) <- xx
  }
  ciMean <- ci

  if(BOOTSTRAP & nk > 10){

    nboot <- 2000
    mu <- matrix(NA,nboot,nc)
    for(g in 1:nboot){
      isamp <- sample(ii,nk,replace=T)
      mu[g,] <- apply(yy[isamp,],2,mean,na.rm=T)
    }

    ciMean <- apply(mu,2,quantile,c(.5,.025,.975),na.rm=T)
    colnames(ciMean) <- xx
  }


  if(PLOT){
    par(bty='n')
    plot(xx,ci[1,],type='l',lwd=2,ylim=c(-.6,.6),ylab='Correlation',xlab='Lag',col=2)
    abline(h=0,lwd=2,col='grey')
    abline(v=0,lwd=2,col='grey')

    for(j in 1:3)lines(xx,ci[j,],lty=2)
    for(j in 1:2)lines(xx,ciMean[j,],lty=2,col=2,lwd=2)

    text(xx[1],.5,paste('n = ',nk),pos=4)
  }

  list(lag = xx, ci = ci,  ciMean = ciMean, n = nk)
}

PDForPS <- function(PDF,fname){

  if(PDF)dev.copy2pdf(file=paste(fname,'.pdf',sep=''))
  if(!PDF)dev.print(device=postscript,file=paste(fname,'.ps',sep=''),width=6,horizontal=F)
}

####################################################



grid2map <- function( xy, z, nx = 20, ny = 20, 
                      colRamp = rampWhite2Green,
                      xlim = NULL, ylim = NULL, zlim = NULL,
                      xaxt = 's', yaxt = 's', units = '', 
                      MAP = T, LEGEND = T, 
                      xleg = 'bottomright', 
                      yleg = c( ylim[1] - .2, ylim[1] + 1 ) ){
  
  # xy comes from plot file or expand.grid( x, y )
  # z is x by y matrix
  
  if( is.null( xlim ) ){
    xlim = range( xy[,1] ) + 1.2*c(-1, 1)
    ylim = range( xy[,2] ) + 1.2*c(-1, 1)
  }
  if( is.null( zlim ) )zlim <- range( z, na.rm = T ) + 1.2*c(-1, 1)
  
  wx <- which( xy[,1] >= xlim[1] & xy[,1] <= xlim[2] &
                 xy[,2] >= ylim[1] & xy[,2] <= ylim[2] )
  xy <- xy[wx,]
  z  <- z[wx]
  
  tgrid <- points2grid( x = xy[,1], y = xy[,2], z = z, nx = nx, ny = ny )
  z <- tgrid$z
  z[ !is.finite(z) ] <- NA
  x <- tgrid$x
  y <- tgrid$y
  
  surf2map( x, y, z, colRamp = colRamp,
            xlim = xlim, ylim = ylim, zlim = zlim,
            xaxt = xaxt, yaxt = yaxt, 
            MAP = MAP, LEGEND = LEGEND, 
            units = units, xleg = xleg, yleg = yleg )
}


points2grid <- function( x, y, z, FUN = 'mean', nx = 50, ny = 50, 
                         xlim = NULL, ylim = NULL, LOG = F ){
  
  # FUN can be 'mean' or 'sum'
  
  KEEPX <- KEEPY <- T
  
  if( is.null(xlim) ){
    xlim <- quantile(x, c(.01, .99), na.rm=T )
    KEEPX <- F
  }
  if( is.null(ylim) ){
    ylim <- quantile(y, c(.01, .99), na.rm=T )
    KEEPY <- F
  }
  
  xseq <- seq( xlim[1], xlim[2], length = nx )
  yseq <- seq( ylim[1], ylim[2], length = ny )
  
  if( !KEEPX ){
    xnew <- xseq
    nnx <- nx
    d1 <- 10
    
    while( nnx == nx ){
      nnx <- length( unique( round(xnew, d1 ) ) )
      d1  <- d1 - 1
    }
    dx <- min( diff( round( xseq, d1+2 ) ) )
    xseq <- sort( round( seq( xlim[1], xlim[2], by = dx ), d1+2 ) )
    nx   <- length(xseq)
  }
  
  if( !KEEPY ){
    ynew <- yseq
    nny <- ny
    d2 <- 10
    while( nny == ny ){
      nny <- length( unique( round(ynew, d2 ) ) )
      d2  <- d2 - 1
    }
    dy <- min( diff( round( yseq, d2+2 ) ) )
    yseq <- sort( round( seq( ylim[1], ylim[2], by = dy ), d2+2 ) )
    ny   <- length(yseq)
  }
  
  ix <- findInterval( x, xseq, all.inside = T )
  iy <- findInterval( y, yseq, all.inside = T )
  
  if( FUN == 'mean' ){
    z0 <- tapply( z, list(x = ix, y = iy), mean, na.rm=T)
  }else{
    z0 <- tapply( z, list(x = ix, y = iy), sum, na.rm=T)
  }
  
  wx <- sort( unique( which( !c(1:nx) %in% as.numeric( rownames(z0) ) ) ) )
  wy <- sort( unique( which( !c(1:ny) %in% as.numeric( colnames(z0) ) ) ) )
  
  if( length(wx) > 0 ){
    xx <- matrix( NA, length(wx), ncol(z0) )
    rownames(xx) <- wx
    z0 <- rbind( z0, xx )
    z0 <- z0[order( as.numeric( rownames(z0)) ),]
  }
  if( length(wy) > 0 ){
    xx <- matrix( NA, nrow(z0), length(wy) )
    colnames(xx) <- wy
    z0 <- cbind( z0, xx )
    z0 <- z0[, order( as.numeric( colnames(z0)) )]
  }
  
  rownames(z0) <- xseq[ as.numeric(rownames(z0)) ]
  colnames(z0) <- yseq[ as.numeric(colnames(z0)) ]
  
  kn <- 4
  
  wna <- which( is.na( z0 ), arr.ind = T )
  if( length(wna) > 0 ){
    wf <- which( is.finite( z0 ), arr.ind = T )
    tmp <- RANN::nn2( wf, wna, k = kn, radius = 2 )
    mm <- tmp[[1]]
    dd <- tmp[[2]] 
    dd[ dd > 3 ] <- NA
    wm <- which( rowSums( (dd*0 + 1), na.rm = T ) > 1 )
    
    fix <- wna[wm,]
    dd  <- dd[wm,]
    mm  <- mm[wm,]
    
    ss <- matrix( z0[ wf[mm,] ], ncol = kn )
    ss[ is.na(dd) ] <- 0
    
    z0[ fix ] <- rowSums( ss/dd^2, na.rm = T )/rowSums( dd^2, na.rm = T )
  }
  
  xseq <- xseq[-1] - .5*diff(xseq)[1] # midpoints
  yseq <- yseq[-1] - .5*diff(yseq)[1]
  
  
  if( LOG )z0 <- log10( z0 )
  
  z0[ !is.finite(z0) ] <- NA
  
  xseq <- as.numeric(rownames(z0))
  yseq <- as.numeric(colnames(z0))
  
  list( x = xseq, y = yseq, z = z0 )
}


predictTrajectory <- function( x, ts, afreq, time, nh = length( ts ), 
                               periods = NULL, verbose = F ) {
  
  # x     - vector returned by stats::fft
  # nh - number of terms used for prediction
  
  N   <- length(ts)
  
  if( nh < N ){
    ord <- order( x )[ (nh+1): N ]
    x[ ord ] <- 0
    if( verbose )print( ts[ x != 0 ]*afreq*time )
  }
  if( !is.null( periods ) ){
    
    ti <- ts*afreq*time
    
    dd <- outer( ti, periods, function( x1, x2 ) ( x1 - x2 )^2 )
    ip <- apply( dd, 2, which.min )
    x[ -ip ] <- 0
  }
  
  i   <- complex(real = 0, imaginary = 1)
  x.n <- rep(0,N)           # create vector to keep the trajectory
  ks  <- 0:(length(x) - 1)
  
  for(n in 0:(N-1)) {       
    x.n[n+1] <- sum(x * exp(i*2*pi*ks*n/N))/N
  }
  #  x.n * afreq 
  x.n
}

plotSpectrum <- function( x, afreq, xlim = c( 0, length(x)) ) {
  
  # x     - vector returned by stats::fft
  # afreq - acquisition frequency, i.e., time interval between observations
  
  plot.data     <- cbind( 0:(length(x) - 1), Mod(x) )
  plot.data[,1] <- plot.data[,1]/afreq/2/pi
  plot.data[ 2:length( x), 2] <- 2*plot.data[2:length(x), 2] 
  
  plot(plot.data, t = "h", lwd = 4, bty = 'n',
       xlab = "Frequency (time)", ylab = "Strength", 
       xlim = xlim, ylim = c(0, max( Mod(plot.data[,2]))))
}


surf2map <- function( x, y, z, colRamp = c('#a6611a','#dfc27d','#80cdc1','#018571'),
                      xlim = NULL, ylim = NULL, zlim = NULL,
                      xaxt = 's', yaxt = 's', 
                      MAP = T, LEGEND = T,
                      units = '', xleg = 'bottomright', 
                      yleg = c( ylim[1] - .2, ylim[1] + 1 ) ){
  
  if( is.null( xlim ) ){
    xlim = range( x ) + .1*c(-1, 1)
    ylim = range( y ) + .1*c(-1, 1)
  }
  wx <- which( x >= xlim[1] & x <= xlim[2] )
  x  <- x[ wx ]
  wy <- which( y >= ylim[1] & y <= ylim[2] )
  y  <- y[ wy ]
  z  <- z[wx, wy]
  
  if( is.null( zlim ) )zlim <- range( z, na.rm = T ) + 1.2*c(-1, 1)
  
  z[ z > zlim[2] ] <- zlim[2]
  
  image( x, y, z, xlim = xlim, ylim = ylim, col = colRamp, bty = 'n', 
         zlim = zlim, xlab = '', ylab = '', asp = 1, xaxt = xaxt, yaxt = xaxt,
         add = F )
  if( MAP ){
    maps::map('world', xlim = xlim, ylim = ylim, add = T, lwd = 2.5, col = 'white' )
    maps::map('world', xlim = xlim, ylim = ylim, add = T, lwd = 2, col = 'grey' )
  }
  
  zlim <- round( zlim, -2 )
  
  if( LEGEND )
    colorScaleLegend( xleg = xleg, yleg = yleg, 
                      xlim = xlim, ylim = ylim, 
                      zlim = zlim,units = units, 
                      colorRamp = colRamp ) # in mastifFunctions.R
}

makeGrid <- function( xlim, ylim, by = 1, byy = by ){
  
  # if byy provided then can have different increment for x and y
  
  require( sf )
  require( spData ) ## For `world`, an sf MULTIPOLYGON object
 # require( spDataLarge )
  require( terra )
  
  if( by < 1 ){
    digits <- nchar( by ) - 1
  }else{
    digits <- nchar( by )
  }
  
  # supplement grid
  xlim[1] <- floor( xlim[1] ) - 1
  xlim[2] <- ceiling( xlim[2] ) + 1
  ylim[1] <- floor( ylim[1] ) - 1
  ylim[2] <- ceiling( ylim[2] ) + 1
  xgrid   <- round( seq( xlim[1], xlim[2], by = by ), digits )
  ygrid   <- round( seq( ylim[1], ylim[2], by = byy ), digits )
  xygrid  <- expand.grid( xgrid, ygrid ) 
  colnames( xygrid ) <- c('lon','lat')
  rownames(xygrid)   <- paste( 'E', xygrid[,1], 'N', xygrid[,2], sep = '' )
  
  ## Create an sf POINTS object
  pts  <- st_as_sf(xygrid, coords=1:2, crs=4326)
  stpt <- st_intersects(pts, world )
  
  wi   <- try( which( !is.na( as.numeric(stpt) ) ), silent = T )
  
  if( inherits( wi, 'try-error' ) ){
    
    countries <- geodata::world(path = tempdir())
    
    # make a polygon map delimiting the entire extent of the Earth:
    earth <- terra::vect(terra::ext(), crs = "EPSG:4326")
    # erase the countries (land parts) to get just the marine polygon:
    marine <- terra::erase(earth, countries)
    
    stpt <- st_intersects(pts, sf::st_as_sf(marine) )
    
    #    sdf <- as.data.frame(stpt)
    
    wi   <- which( !is.na( as.numeric(stpt) ) )
  }
  
  xygrid[wi,]
}


points2grid <- function( x, y, z = NULL, FUN = 'mean', nx = 50, ny = 50, 
                         nk = 10, dwt = 1, dmax = 5,
                         xlim = NULL, ylim = NULL, LOG = F ){
  
  # FUN can be 'mean' or 'sum'
  # nk is the number of neighbors
  # dwk is the exponential weight given distance as exp( -(dwt*distance)^2 )
  # dmax is maximum distance
  
  require( RANN )
  
  KEEPX <- KEEPY <- T
  tiny  <- 1e-10
  tinc  <- c( tiny, 1 - tiny )
  
  if( is.null( z ) )z <- rep( 1, length(x) )
  
  if( is.null(xlim) ){
    xlim <- quantile(x, tinc, na.rm=T )
    #   KEEPX <- F
  }
  if( is.null(ylim) ){
    ylim <- quantile(y, tinc, na.rm=T )
    #   KEEPY <- F
  }
  
  xseq <- seq( xlim[1], xlim[2], length = nx )
  yseq <- seq( ylim[1], ylim[2], length = ny )
  
  by  <- signif( diff( xseq )[1], 3 )
  byy <- signif( diff( yseq)[1], 3 )
  xygrid <- makeGrid( xlim, ylim, by = by, byy = byy )
  
  xy <- cbind( x, y )
  xycols <- c( 'lon', 'lat' )
  colnames( xy ) <- xycols
  tmp <- nn2( xy[ , xycols ], xygrid[ ,xycols], k = nk )
  
  wsum <- ssum <- rep( 0, nrow( xygrid ) )
  
  for( k in 1:nk ){
    
    dk <- tmp[[2]][,k]
    wk <- exp( -(dwt*dk)^2 )
    wk[ dk > dmax ] <- 0
    zk <- z[ tmp[[1]][,k] ]
    wsum <- wsum + wk
    ssum <- ssum + wk*zk
  }
  z <- ssum/wsum
  xygrid <- cbind( xygrid, z )
  
  xc <- sort( unique( xygrid$lon ) )
  yc <- sort( unique( xygrid$lat ) )
  
  zmat <- matrix(NA, length( xc ), length( yc ), dimnames = list( xc, yc ) )
  
  mr <- match( as.character( xygrid$lon ), rownames( zmat ) )
  mc <- match( as.character( xygrid$lat ), colnames( zmat ) )
  zmat[ cbind( mr, mc ) ] <- xygrid$z
  
  list( x = xc, y = yc, z = zmat, xyz = xygrid )
}


points2grid2 <- function(xx,yy,grid,wt=1){
  
  # grid is 2 cols (xgrid,ygrid)
  # wt are weigths applied to each point
  
  require(RANN)
  
  tmp <- nn2(grid,cbind(xx,yy),k=1)[[1]]
  out <- grid[tmp,]
  tt  <- unique(tmp)
  fraction <- length(tt)/nrow(grid)
  
  xv  <- sort(unique(grid[,1]))
  yv  <- sort(unique(grid[,2]))
  
  i  <- match(out[,1],xv)
  j  <- match(out[,2],yv)
  
  mat <- .byIndex(i*0 + wt,INDICES=list('x' = i, 'y' = j),sum,coerce=T)
  ix  <- as.numeric(rownames(mat))
  iy  <- as.numeric(colnames(mat))
  
  rownames(mat) <- signif(xv[ix],3)
  colnames(mat) <- signif(yv[iy],3)

  list(mat = mat, outValues = out, gridIndex = tmp, gridFraction = fraction)
}
  

points2contour <- function(x,y,q=NULL,xlabel=NULL,ylabel=NULL,main=NULL,
                           xp=NULL,yp=NULL,levs=NULL,add=F,maxz=Inf,col='black',fill=F){  

  #creates contours for (x,y) at density q
  
  if(is.null(q))q <- 10
  if(is.null(xp))xp <- range(x,na.rm=T)
  if(is.null(yp))yp <- range(y,na.rm=T)

  xr    <- range(x,na.rm=T)
  yr    <- range(y,na.rm=T)
  xd <- (xr[2] - xr[1])/q
  yd <- (yr[2] - yr[1])/q

  xgrid <- seq(xr[1],xr[2],by=xd)
  ygrid <- seq(yr[1],yr[2],by=yd)

  xf <- cut(x,xgrid)
  yf <- cut(y,ygrid)

  z <- table(xf,yf)
  z[z > maxz] <- maxz

  xmids <- (xgrid - xd/2)[-1]
  ymids <- (ygrid - yd/2)[-1]
  
  d     <- max(c(diff(xr),diff(yr)),na.rm=T)
  if(is.null(q))q <- d/20
  
  if(is.null(levs))levs <- signif(seq(min(z),max(z),length.out=4),1)
  
  lwdd <- 2
  if(length(levs) > 1)lwdd <- seq(2,length(levs),by=1)
  
  if(!add)image(xmids,ymids,z,xlab=xlabel,ylab=ylabel,xlim=c(xp[1],xp[2]),ylim=c(yp[1],yp[2]))
  contour(xmids,ymids,z,add=T,levels=levs,lwd=lwdd,col=col)
  if(fill).filled.contour(xmids,ymids,z,levels=levs,col=col)
  if(!is.null(main))title(main)
  
  invisible( list(x = xmids, y = ymids, z = z) )
}

####################################################

histf <- function(vec,minv,maxv)hist(vec,breaks=seq(minv,maxv,by=.02))


####################################################

myrmultinom <- function(size,p,ASVECTOR=F){  

  # n multinomial r.v. for a n by ncol(p) matrix of probs
  # each row of p is a probability vector
  # size is one integer or a length-n vector of integers
  # if ASVECTOR = T all size == 1, returns a vector of columns, otherwise a matrix

  p <- row2Mat(p)

  n     <- nrow(p)
  J     <- ncol(p)

  if(length(size) == 1)size <- rep(size,n)

  jord  <- sample(J,J)    #randomize order
  
  rs <- rowSums(p)
  ws <- which(rs != 1)
  if(length(ws) > 0){
    p[ws,] <- p[ws,]/rs[ws]
  }

  p <- row2Mat(p[,jord])

  sizej <- size
  sumj  <- rep(0,n)
  dpj   <- rep(1,n)
  pj    <- p
  wj    <- c(1:n)
  
  if(ASVECTOR){        #  only if all size == 1
    
    yy <- size*0
    
    for(j in 1:(J-1)){
      a     <- round(pj[wj,1],10)
      tmp  <- rbinom(length(wj),sizej[wj],a)
      yy[wj[tmp == 1]] <- j
      sumj[wj]  <- sumj[wj] + tmp
      sizej <- size - sumj                       # no. remaining to choose
      dpj   <- dpj - p[,j]                       # Pr for remainder
      pj    <- matrix(p[,c((j+1):J)]/dpj,nrow(p))
      wj    <- which(sumj < size,arr.ind=T) 
    }
    
    yy[yy == 0] <- J

    return(yy)
  }
  
  yy  <- matrix(0,n,J)

  for(j in 1:(J-1)){
    a     <- round(pj[wj,1],10)
    yy[wj,j] <- rbinom(length(wj),sizej[wj],a)
    sumj  <- sumj + yy[,j]
    sizej <- size - sumj                       # no. remaining to choose
    dpj   <- dpj - p[,j]                       # Pr for remainder
    pj    <- matrix(p[,c((j+1):J)]/dpj,nrow(p))
    wj    <- which(sumj < size,arr.ind=T) 
  }

  if(n == 1)yy[,J] <- size - sum(yy)
  if(n > 1) yy[,J] <- size - rowSums(yy)

  yy[,jord] <- yy
  yy
}

####################################################

truncpars <- function(x,lo,hi){        

  #JS Clark
  #fit truncated multivariate normal to posterior x, known lo and hi

  if(is.vector(x))x <- matrix(x,length(x),1)  
  px <- ncol(x)          #dimension of original x

  ww <- which(!is.finite(x),arr.ind=T)
  if(length(ww) > 0){
    ww <- unique(ww[,2])
    wk <- which(!c(1:ncol(x)) %in% ww,arr.ind=T)
    x <- x[,wk]
  }
  muvec <- apply(x,2,mean,na.rm=T)
  sig   <- cov(x,use="complete.obs")
  pk <- ncol(x)
  nn <- nrow(x)

  ngg   <- 1000
  nkeep <- ngg - 300
  nk    <- 0
  mug   <- rep(0,pk)
  cvg   <- rep(0,pk^2)

  for(g in 1:ngg){
    tmp   <- truncmvtnorm(x,muvec,sig,lo,hi)
    muvec <- tmp$mu
    sig   <- tmp$sig

    if(g > nkeep){
      print(muvec)
      nk       <- nk + 1
      mug      <- mug + muvec
      cvg      <- cvg + as.vector(sig)
    }
  }
  mvec <- mug/nk
  covmat <- matrix(cvg/nk,pk,pk)

  if(length(ww) > 0){
    mnew     <- rep(NA,px)
    mnew[wk] <- mvec
    covnew   <- matrix(NA,px,px)
    covnew[wk,wk] <- covmat
    mvec   <- mnew
    covmat <- covnew
  }

  list(mu = mvec, cm = covmat)
}

####################################################

trunclogis <- function(n,lo,hi,bpars,xvars,wp){ 

  #truncated logistic
  # bpars - parameter vector for logit
  # xvars - variables corresponding to bpars, e.g., c(1,x1,x2)
  # wp    - which variable and parameter in the logit model 
  #         (not 1, because bpars[1] is intercept)
  #         (xvars[wp] is not used)

  xlo     <- xvars
  xlo[wp] <- lo
  xhi     <- xvars
  xhi[wp] <- hi
  sl      <- sum(xlo*bpars)
  sh      <- sum(xhi*bpars)
  lflo    <- exp(sl)/(1 + exp(sl))
  lfhi    <- exp(sh)/(1 + exp(sh))

  z <- runif(n,lflo,lfhi)
  (log(z/(1 - z)) - sum(xvars[-wp]*bpars[-wp]))/bpars[wp]
}

####################################################

truncmvtnorm <- function(x,muvec,sig,lo,hi){

  #sample from a truncated normal

  if(is.vector(x))x <- matrix(x,length(x),1)
  if(length(sig) == 1)sig <- matrix(sig,1,1)
  y  <- x*0
  pk <- ncol(x)
  n  <- nrow(x)

  for(j in 1:pk){

    if(j == 1){
      t <- muvec[1]
      w <- sqrt(sig[1,1])
    }
    if(j > 1){

       svec <- c(1:(j-1))
       vmat <- sig[j,svec] %*% solve(sig[svec,svec])
       ymu  <- y[,svec] - matrix(rep(muvec[svec],n),n,(j-1),byrow=T)
       t    <- t(muvec[j] +  vmat %*% t(ymu))
       w    <- as.numeric(sqrt( sig[j,j] - vmat %*% c(sig[svec,j]) ))
    }

    up <- pnorm(x[,j],t,w) - pnorm(lo[j],t,w)
    do <- pnorm(hi[j],t,w) - pnorm(lo[j],t,w)

    add <- w*qnorm(up/do)
    add[!is.finite(add)] <- 0
    y[,j] <- t + add

  }

  muy   <- colMeans(y)
  muvec <- myrmvnorm(1,muy,sig/n)

 #the covariance matrix

  mumat <- matrix(muvec,n,pk,byrow=T)
  sy    <- crossprod(y - mumat) 
   ss   <- solve(sy)
   df   <- n 
  alpha <- myrmvnorm(df,rep(0,pk),ss)
  th    <- solve(crossprod(alpha))

  list(mu = muvec, sig = th)

} 

####################################################

logit <- function(x){log(x/(1-x))}  #logit

####################################################

invlogit <- function(x, log = FALSE){  #inverse logit

  if(log)return(-log(1 + exp(-x)))
  1/(1 + exp(-x))

}

##############################
invMat <- function(SS,NEARPD=F){  #matrix inversion, if NEARPD find closest PD matrix
    
    testv <- try(chol(SS),T)

    if( inherits(testv,'try-error') ){
       message('near pos definite used in invMat')
       if(NEARPD){
         SS     <- as.matrix( nearPD(SS)$mat )
         testv <- try(chol(SS),T)
       }
    }

    chol2inv(testv)
}

mydmvnorm <- function(xx,mu,sigma,sinv=NULL,log=FALSE){
  
  xx <- xx - mu
  if(!is.matrix(xx))xx <- matrix(xx,1)
  
  if(!is.null(sinv)){
    distval <- diag( xx%*%sinv%*%t(xx) )
    ev      <- eigen(sinv, only.values = TRUE)$values
    logdet  <- -sum(log(ev))
  }
  
  if(is.null(sinv)){
    testv <- try(chol(sigma),T)
    if(inherits(testv,'try-error')){
      tiny <- min(abs(xx))/100 + 1e-5
      sigma <- sigma + diag(diag(sigma + tiny))
      testv <- try(chol(sigma),T)
    }
    covm <- chol2inv(testv)
    distval <- rowSums((xx %*% covm) * xx)
    ev <- eigen(sigma, only.values = TRUE)$values 
    if(min(ev) < 0)ev <- nearPD(sigma)$eigenvalues
    logdet   <- sum(log( ev ))
  }
  
  logret <- -(ncol(xx) * log(2 * pi) + logdet + distval)/2
  if(log)return(logret)
  exp(logret)
}
  
################################################
mydmvnormOld <- function(xx,mu,sigma,log=FALSE){

  #mv normal density

    if (is.vector(xx))xx <- matrix(xx, ncol = length(xx))
    if (is.vector(mu))mu <- matrix(mu, ncol = length(xx))

 #   zz   <- sweep(xx, 2, mu)
    zz <- xx - mu
    
    ss <- colSums(zz^2)%*%sigma

    testv <- try(chol(sigma),T)
    if(inherits(testv,'try-error')){
       message('error in mydmvnorm')
       return(mu)
    }

    cov <- chol2inv(testv)
    distval <- rowSums((zz %*% cov) * zz)
    names(distval) <- rownames(zz)

  #  distval <- mahalanobis(zz, mu, sigma)
    logdet   <- sum(log( eigen(sigma, symmetric = TRUE, only.values = TRUE)$values ))
    if(is.na(logdet))return(logdet)
    logretval <- -(ncol(zz) * log(2 * pi) + logdet + distval)/2
    if(log)return(logretval)
    exp(logretval)
}


bUpdateMVN <- function(xx,yy,bb = NULL,lo = NULL, hi = NULL,sigma){

  #  yy ~ MVN(xx%*%bb,sigma)
  #   xx - design matrix (n X p)
  #   yy - response matrix (n X q)
  #   bb - parameter matrix to update (p X q)
  #   lo - lower truncation (p X q)
  #   hi - upper truncation (p X q)
  #   if truncation (lo and/or hi provided) bb must also be provided

  testv <- try(chol(crossprod(xx)),T)
  if(inherits(testv,'try-error')){
    message('X not full rank')
    return( bb ) 
  }
 
  cx   <- chol2inv(testv)
  mu   <- matrix(cx %*% crossprod(xx,yy),1)
  vv   <- kronecker(sigma,cx)
  
  if(is.null(lo)){
    b <- matrix( rmvnorm(1,mu,vv),ncol=2)
    dimnames(b) <- dimnames(bb)
    return(b)
  }

  lo <- as.vector(lo)
  hi <- as.vector(hi)

  b <- matrix( tnorm.mvt(as.vector(bb),mu,vv,lo,hi,times=3), ncol=2)
  dimnames(b) <- dimnames(bb)
  b
}


##############################################
rbvnormFromCols <- function(muMat,sigMat,lo=NULL,hi=NULL){  
	
	#one random vector per row for bivariate normal, based on cholesky
	#sigMat has 4 cols (s11,s12,s12,s22)
	
	n   <- nrow(muMat)
	rho <- sigMat[,2]/sqrt(sigMat[,1]*sigMat[,4])
	
        if(is.null(lo)){
	  z1 <- rnorm(n)
	  z2 <- rnorm(n)
        }
        if(!is.null(lo)){
          l1 <- (lo - muMat[,1])/sqrt(sigMat[,1])
          h1 <- (hi - muMat[,1])/sqrt(sigMat[,1])
          z1 <- tnorm(n,l1,h1,0,1)

          l2 <- ((lo - muMat[,2])/sqrt(sigMat[,4]) - rho*z1)/sqrt(1 - rho^2)
          h2 <- ((hi - muMat[,2])/sqrt(sigMat[,4]) - rho*z1)/sqrt(1 - rho^2)
          z2 <- tnorm(n,l2,h2,0,1)
        }
	
	cbind(muMat[,1] + sqrt(sigMat[,1])*z1,
	      muMat[,2] + sqrt(sigMat[,4])*(rho*z1 + sqrt(1 - rho^2)*z2) )
}

dbvnormFromCols <- function(y,muMat,sigMat){	
  
	sinv <- invertcol2(sigMat)
	z    <- y - muMat
	
	q <- y*0
	q[,1] <- z[,1]*sinv[,1] + z[,2]*sinv[,2]
	q[,2] <- z[,1]*sinv[,3] + z[,2]*sinv[,4]
	q     <- q[,1]*z[,1] + q[,2]*z[,2]
	
	logdet <- log(sigMat[,1]*sigMat[,4] - 2*sigMat[,2])
	
	-(log(2 * pi) + logdet + q)/2
	
}
	


pmvnormApprox <- function(q,mu=rep(0,nrow(sigma)),sigma ,times=100){
  
  #'joint' - pr that each vector < q and
  #'marginal' - pr each each element marginally, i.e., pnorm(q,mu,sqrt(diag(sigma)) )

if(!is.matrix(mu))mu <- matrix(mu,1)
if(length(q) == 1)q <- rep(q,nrow(sigma))

nr <- nrow(mu)
nc <- nrow(sigma)

avec <-  mu
z    <- mu*0

joint <- rep(0,nr)

for(i in 1:times){
  
  zi <- mu*0
  
  for(k in 1:nc){
    
    #   r <- tnorm.mvtRcpp(mu, mu, sigma, lo= -Inf, hi=Inf, whichSample=k, times=1) 
    
    tmp <- conditionalMVNRcpp(avec,mu,sigma,k)
    muk <- tmp$mu
    sgk <- tmp$vr
    
    if(length(muk) == 0)next
    
    avec[,k] <- rnorm(nr,muk,sqrt(sgk))
    wj  <- which(avec[,k] < q[k])
    zi[wj,k] <- zi[wj,k] + 1
    z[wj,k]  <- z[wj,k] + 1
  }
  joint[rowSums(z) == nc] <- joint[rowSums(z) == nc] + 1
}

 list(joint = joint/times, marginal = z/times)
}


####################################################

myrmvnormOld <- function (nn, mu, sigma){
  
  
  testv <- try(svd(sigma),T)
  
  if( inherits(testv,'try-error') ){
    
    ev <- eigen(sigma, symmetric = TRUE)
    testv <- t(ev$vectors %*% (t(ev$vectors) * sqrt(ev$values)))
  } else {
    testv <- t(testv$v %*% (t(testv$u) * sqrt(testv$d)))
  }
  
    retval <- matrix(rnorm(nn * ncol(sigma)), nn) %*% testv
    retval + mu
}
################################

myrmvnorm <- function (nn, mu, sigma){
  
  # nn - no. samples from one mu vector or nrow(mu) for matrix
  
  if(!is.matrix(mu)) mu <- matrix(mu,1)
  if(length(mu) == 1)mu <- matrix(mu,1,nrow(sigma))
  if(ncol(mu) == 1)  mu <- t(mu)
  
  m <- ncol(sigma)
  
  if(ncol(mu) != m)stop('dimension mismatch mu, sigma')
  
  if(nn > 1 & nrow(mu) == 1)mu <- matrix(mu,nn,m,byrow=T)
  
  if(nn != nrow(mu))stop('sample size does not match mu')
  
  testv <- try(svd(sigma),T)
  
  if( inherits(testv,'try-error') ){
    ev <- eigen(sigma, symmetric = TRUE)
    testv <- t(ev$vectors %*% (t(ev$vectors) * sqrt(ev$values)))
  } else {
    testv <- t(testv$v %*% (t(testv$u) * sqrt(testv$d)))
  }
  
  p <- matrix(rnorm(nn * m), nn) %*% testv
  p + mu
}

####################################################

mypmvnorm <- function (lower, upper, mu,sigma){

   corr <- cov2cor(sigma)
   lower <- (lower - mu)/sqrt(diag(sigma))
   upper <- (upper - mu)/sqrt(diag(sigma))
   mean <- rep(0, length(lower))
   RET  <- mvt(lower = lower, upper = upper, df = 0,
                corr = corr, delta = mu, maxpts = 25000, abseps = 0.001,
                releps = 0)
   return(RET$value)
}
####################################################

myqmvnorm <- function (p, interval = c(-10, 10), tail = c("lower.tail", "upper.tail",
    "both.tails"), mean = 0, corr = NULL, sigma = NULL, maxpts = 25000,
    abseps = 0.001, releps = 0, ...)
{
    if (length(p) != 1 || (p <= 0 || p >= 1))
        stop(sQuote("p"), " is not a double between zero and one")
    tail <- match.arg(tail)
    dim <- length(mean)
    if (is.matrix(corr))
        dim <- nrow(corr)
    if (is.matrix(sigma))
        dim <- nrow(sigma)
    lower <- rep(0, dim)
    upper <- rep(0, dim)
    args <- checkmvArgs(lower, upper, mean, corr, sigma)
    dim <- length(args$mean)
    pfct <- function(q) {
        switch(tail, both.tails = {
            low <- rep(-abs(q), dim)
            upp <- rep(abs(q), dim)
        }, upper.tail = {
            low <- rep(q, dim)
            upp <- rep(Inf, dim)
        }, lower.tail = {
            low <- rep(-Inf, dim)
            upp <- rep(q, dim)
        }, )
        pmvnorm(lower = low, upper = upp, mean = args$mean, corr = args$corr,
            sigma = args$sigma, abseps = abseps, maxpts = maxpts,
            releps = releps) - p
    }
    if (tail == "both.tails") {
        interval[1] <- 0
        interval <- abs(interval)
    }
    qroot <- uniroot(pfct, interval = interval, ...)
    names(qroot)[1:2] <- c("quantile", "f.quantile")
    qroot
}


####################################################

fit.tnorm <- function(parvec){  

  #fit parameters for truncated normal
  #vector parvec includes mean, diagonal, & offdiagonals for cov matrix

  mu <- parvec[mi]
  sigmat <- matrix(0,length(mu),length(mu))
  sigmat <- diag(parvec[p2])
  sigmat[si] <- parvec[p3]
  sigmat[cbind(si[,2],si[,1])] <- parvec[p3]

  c1     <- mydmvnorm(x,mu,sigmat,log=T)
  c2     <- log(mypmvnorm(lo,hi,mu,sigmat))
  -sum(c1 - c2)

}

patch2patch <- function(from,to,q,dx,sigma){

   n2   <- length(to)
   n1   <- length(from)
   r    <- matrix(q[to],n2,n1)*exp(-(dx/sigma)^2)

   if(n1 > 1) p <- r/matrix(colSums(r),n2,n1,byrow=T)
   if(n1 == 1)p <- r/sum(r)
   p
}

move <- function(from,to,q,dx,sigma){
   pr  <- patch2patch(from,to,q,dx,sigma) #prob moving to other patches
   wk  <- myrmultinom(1,t(pr))            #movement vector
   which(wk == 1,arr.ind=T)[,2]           #extract new patch location
} 


####################################################

smooth.na <- function(x,y){   

  #remove missing values
  #x is the index
  #y is a matrix with rows indexed by x
  
  if(!is.matrix(y))y <- matrix(y,ncol=1)

    wy <- which(!is.finite(y),arr.ind =T)
    if(length(wy) == 0)return(cbind(x,y))
    wy <- unique(wy[,1])
      ynew <- y[-wy,]
      xnew <- x[-wy]

    return(cbind(xnew,ynew))
}
####################################################

smooth.ma <- function(y,wt){   

  #moving average filter with weights (w0,w1,...), assumed symmetric

  if(length(wt) > length(y))wt <- wt[1:length(y)]
  nw <- length(wt)
  ny <- length(y)
  w <- c(rev(wt),wt[2:nw])
  ymat <- matrix(NA,ny,length(w))
  kb <- nw
  ke <- ny
  ky <- ny - kb + 1

  kb <- c(nw:-(nw-2)); kb[kb < 1] <- 1
  ke <- c((ny+nw-1): (ny-nw+1)); ke[ke > ny] <- ny
  yb <- rev(kb)
  ye <- rev(ke)

  for(kj in 1:(2*nw-1)){
    ymat[kb[kj]:ke[kj],kj] <- y[yb[kj]:ye[kj]]
  }

  wmat <- matrix(rep(w,ny),ny,length(w),byrow=T)
  wmat[which(is.na(ymat))] <- NA
  wmat <- wmat/rowSums(wmat,na.rm=T)
  newy <- rowSums(ymat*wmat,na.rm=T)
  newy

}

eucDist <- function(x1, x2){
  
  #euclidean distance between two matrices with same number of columns
  #x1 is n X p
  #x2 is m X p
  #distance is n X m
  
  p <- ncol(x1)
  d <- matrix(0,nrow(x1),nrow(x2))
  
  for(j in 1:p){
    d <- d + outer(x1[,p],x2[,p],function(x1,x2) (x1 - x2)^2)
  }
  sqrt(d)
}


####################################################

distmat <- function(x1,y1,x2,y2){
    xd <- outer(x1,x2,function(x1,x2) (x1 - x2)^2)
    yd <- outer(y1,y2,function(y1,y2) (y1 - y2)^2)
    t(sqrt(xd + yd)) 
}
####################################################

tnorm <- function(n,lo,hi,mu,sig){   

  #normal truncated lo and hi

  if(length(lo) == 1 & length(mu) > 1)lo <- rep(lo,length(mu))
  if(length(hi) == 1 & length(mu) > 1)hi <- rep(hi,length(mu))

  q1 <- pnorm(lo,mu,sig)
  q2 <- pnorm(hi,mu,sig) 

  z <- runif(n,q1,q2)
  z <- qnorm(z,mu,sig)
  z[z == Inf]  <- lo[z == Inf]
  z[z == -Inf] <- hi[z == -Inf]
  z
}

nextTimeMat <- function(mat,last=rep(ncol(mat),nrow(mat)),INC=F){

  # mat is n by time matrix, returns same shifted to left

  til <- cbind(c(1:nrow(mat)),last)

  x <- cbind(mat[,-1],mat[,ncol(mat)])
  x[ til ] <- mat[ til ]

  if(INC){
    inc <- x[ cbind(c(1:nrow(mat)),last-1) ] - x[ cbind(c(1:nrow(mat)),last-2) ]
    x[ til ] <- x[ til ] + inc
  }

  x
}

lastTimeMat <- function(mat,first=rep(1,nrow(mat)),INC=F,minInc=.001){  

  # mat is n by time matrix, returns same shifted to right
  # first repeats first value at the first time for obs i 

  nc  <- ncol(mat)
  tif <- cbind(c(1:nrow(mat)),first)
 
  x <- cbind(mat[,1],mat[,-ncol(mat)])
  x[ tif ] <- mat[ tif ]

  if(INC){   #increment first value

    f2 <- first+2
    f1 <- first+1
    
    f2[f2 > nc] <- nc
    f1[f1 >= f2] <- f2[f1 >= f2] - 1

    inc <- x[ cbind(c(1:nrow(mat)),f2) ] - x[ cbind(c(1:nrow(mat)),f1) ]
    inc[inc < minInc] <- minInc
    wf  <- which(first > 1)
    ff  <- x[ tif[wf,] ] - inc[wf]
    ff[ff < minInc] <- minInc
    x[ tif[wf,] ] <- ff

  }
    
  x
}


diffTimeMat <- function(mat,index,first = 0){  # increment matrix from obs by time matrix, pad last value

  md <- mat - lastTimeMat(mat,first) 
  vc <- md[index]
  list(dmat = md, vec = vc)
}


####################################################

tnorm.mvt <- function(avec,muvec,smat,lo=rep(-Inf,length(avec)),hi=rep(Inf,length(avec)),
                      whichSample=c(1:length(avec)),times=1){   

  # truncated multvariate normal
  # muvec is the vector of means
  # smat is the covariance matrix 
  # whichSample indicates which variables to sample

  if(length(lo) == 1)lo <- rep(lo,length(avec))
  if(length(hi) == 1)hi <- rep(hi,length(avec))

  for(j in 1:times){
   for(k in whichSample){

    tmp <- conditionalMVN(avec,muvec,smat,k)
    muk <- tmp$mu
    sgk <- tmp$vr

    if(length(muk) == 0)next

    avec[k] <- tnorm(1,lo[k],hi[k],muk,sqrt(sgk))
   }
  }
  avec
}

linearModel <- function(x,y,PLOT=F){  #objects for a linear regression
  
  n  <- length(y)
  p  <- ncol(x)
  b  <- solve(crossprod(x))%*%crossprod(x,y)  #coefficients 
  py <- x%*%b                                 #predictions
  
  if(PLOT){
    plot(y,py,xlab='Observed',ylab='Predicted')
    abline(0,1,lty=2)
  }
  
  s    <- (crossprod(y - py)/(n - p))[1]      #variance estimate
  varb <- s*solve(crossprod(x))               #parameter covar
  se   <- sqrt(diag(varb))                    #standard errors
  
  ci <- cbind( b - 1.96*se, b + 1.96*se )
  colnames(ci) <- c('.025','.975')
  
  pmat <- cbind(b,se,ci)
  colnames(pmat)[1] <- 'estimate'
  
  list(coeff = signif(pmat,4), shat = s, parVar = signif(varb,3))
}


priorMVReg <- function(x,y,sigma){

  #Minka (2001) Bayesian Linear Regression

  # V - error covariance (sigma)
  # m - columns in x
  # d - columns in y
  #

  m <- ncol(x)
  d <- ncol(y)

  xx  <- crossprod(x)
  xy  <- crossprod(x,y)
  ixx <- solve(xx)

  p1    <- solve(sigma)%*%t(xy)%*%ixx%*%xy
  p2    <- m*d
  alpha <- p2/(sum(diag(p1)) + p2) 

  priorAcov <- kronecker(sigma,ixx/alpha)

  list(alpha = alpha, priorAcov = priorAcov)
}


#######################################################

dmvnormLog <- function(x,mu,sigma){    #multivariate normal 

  #mv normal density

    if (is.vector(x))x <- matrix(x, ncol = length(x))

    distval <- mahalanobis(x, mu, sigma)
    logdet  <- sum(log(eigen(sigma, symmetric = TRUE, only.values = TRUE)$values))
    logretval <- -(ncol(x) * log(2 * pi) + logdet + distval)/2
    logretval
}


#############################################################

rwish <- function(df,SS){

  z  <- matrix(rnorm(df*nrow(SS)),df,nrow(SS))%*%chol(SS)
  crossprod(z)
}

#####################################################
riwish <- function(v,S){

  solve(rwish(v,solve(S)))
}


pmvnormCond <- function(q=0,xx,mu,smat,whichVar=c(1:nrow(smat)),log=F){   
  
  # conditional multvariate normal for Pr of each individual in vector conditioned on others
  # xx - matrix of values, n X p
  # mu - matrix of means, n X p
  # smat is the covariance matrix, p X p
  # whichVar indicates which variables to evaluate
  
  pj <- mu*0
  
  for(k in whichVar){
    tmp <- conditionalMVNVec(xx,mu,smat,k)
    muk <- tmp$mu
    sgk <- tmp$vr
    pj[,k] <- pnorm(q,muk,sd=sqrt(sgk),lower.tail=T,log.p=log)
  }
  pj
}




###################################
conditionalNorm <- function(x,mu,sigMat,cindex){  # bi- or trivariate
  
  if(ncol(x) == 2)return( conditionalBiVarNorm(x,mu,sigMat,cindex) )
  if(ncol(x) == 3)return( conditionalTriVarNorm(x,mu,sigMat,cindex) )
}
  
###########################################
conditionalTriVarNorm <- function(x,mu,sigMat,cindex){
  
  # cindex: conditional for this variable, condition on other variable
  # sigMat is nine columns
  
  nn <- nrow(x)
  nc <- length(cindex)
  r <- ncol(x)
  
  rs <- c(1:r^2)
  rm <- cm <- matrix(rs,r,r)
  tm <- rm
  tm[-cindex,] <- NA
  tm[,cindex]  <- NA
  rm[cindex,] <- rm[,cindex] <- NA
  cm[-cindex,] <- cm[,-cindex] <- NA
  
  notC <- sort(rm[is.finite(rm)])
  isC  <- sort(cm[is.finite(cm)])
  coC  <- sort(tm[is.finite(tm)])
  
  sinvMat <- invertcol3(sigMat)$I
  
  if(length(cindex) == 1){
    ct  <- c(1:r)[-cindex]
    p1  <- sigMat[,coC[1]] * sinvMat[,isC]
    p2  <- sigMat[,coC[2]] * sinvMat[,isC]
    mu1 <- mu[,cindex] + p1*(x[,ct[1]] - mu[,ct[1]]) + p2*(x[,ct[2]] - mu[,ct[2]])
    vr <- sigMat[,isC]
  }
  if(length(cindex) == 2){
    p1  <- sigMat[,coC[1]] * sinvMat[,isC[1]] + sigMat[,coC[2]] * sinvMat[,isC[2]]
    p2  <- sigMat[,coC[1]] * sinvMat[,isC[3]] + sigMat[,coC[2]] * sinvMat[,isC[4]]
    mu1 <- mu[,cindex[1]] + p1*(x[,-cindex] - mu[,-cindex])
    mu2 <- mu[,cindex[2]] + p2*(x[,-cindex] - mu[,-cindex])
    mu1  <- cbind(mu1,mu2)
    vr1 <- sigMat[,isC[1]] - p1*sigMat[,coC[1]]
    vr2 <- sigMat[,isC[2]] - p1*sigMat[,coC[2]]
    vr3 <- sigMat[,isC[3]] - p2*sigMat[,coC[1]]
    vr4 <- sigMat[,isC[4]] - p2*sigMat[,coC[2]]
    vr  <- cbind(vr1,vr2,vr3,vr4)
  }
  
  list(mu = mu1, vr = vr)
}
###########################################
conditionalBiVarNorm <- function(xx,mu,sigMat,cindex){

 # cindex: conditional for this variable, condition on other variable
 # sigMat is four columns

  isC  <- 1
  notC <- 4
  if(cindex == 2){
    isC <- 4
    notC <- 1
  }

  sinv <- 1/sigMat[,notC]

  p1  <- sigMat[,2]*sinv

  mu1 <- mu[,cindex] + p1*(xx[,-cindex] - mu[,-cindex])
  vr1 <- sigMat[,isC] - p1*sigMat[,2]

  list(mu = mu1, vr = vr1)
}
  ##################################3

conditionalMVNcdf <- function(up, xx, mu,sigma,k){  
  
  # xx and mu are n by p matrices for which we want the conditional CDF for k
  # k is index for conditional
  
  testv <- try(chol(sigma[-k,-k]),T)
  if(inherits(testv,'try-error')){
    return( list(mu = numeric(0), vr = numeric(0)) )
  }
  sin <- chol2inv(testv)
   
  p1  <- sigma[k,-k]%*%sin
  mu1 <- mu[,k] + t( p1%*%t(xx[,-k] - mu[,-k]) )
  vr1 <- sigma[k,k] - p1%*%sigma[-k,k]
  
  pnorm(up,mu1,sqrt(vr1))
}

condCov1 <- function(sigma,cdex){
  
  # cdex - condition for these vars
  # p    - condition on these vars
  
  S <- nrow(sigma)
  p <- c(1:S)[-cdex]   
  
  testv <- try(chol(sigma),T)
  # testv <- try(chol(sigma[-cdex,-cdex]),T)
  if(inherits(testv,'try-error')){
    if(SOLVESTOP) return( list(mu = numeric(0), vr = numeric(0)) )
    if(!SOLVESTOP)return( list(mu = mu[,cdex],vr = var(mu[,cdex],na.rm=T)) )
  }
  sinv <- chol2inv(testv)
  solve(sinv[-p,-p])
}

condCov2 <- function(sigma,cdex){
  
  S <- nrow(sigma)
  p <- c(1:S)[-cdex]
 
  testv <- try(chol(sigma[-cdex,-cdex]),T)
  if(inherits(testv,'try-error')){
    if(SOLVESTOP) return( list(mu = numeric(0), vr = numeric(0)) )
    if(!SOLVESTOP)return( list(mu = mu[,cdex],vr = var(mu[,cdex],na.rm=T)) )
  }
  sinv <- chol2inv(testv)
  p1  <- sigma[cdex,-cdex]%*%sinv
  sigma[cdex,cdex] - p1%*%sigma[-cdex,cdex]
}

#############################################################
directIndirectCoeffs <- function( snames=specs, xvector,
                                  bchain=chainList$betaChain,
                                  schain=chainList$sigmaChain, minvalue=.001,
                                  reverseSign=character(0),MEAN=T,
                                  keepNames=NULL, sum2zero=NULL,
                                  omitY = NULL, standardX = NULL, nsim=500){
  
  # if MEAN, then use means, otherwise median
  # indirect do not change with x, can choose not to calculate
  # sum2zero - multilevel factors are departure from mean 
  #          - a list of vectors, one for each multilevel factor, 
  #            e.g., "sum2zero <- list('hosts' = hostNames)", 
  #            where hostNames appear in colnames of bchain
  #indirFrom - effect from all others
  #indirTo   - effect on all others
 
  if(is.matrix(xvector))stop('xvector should be a single row vector')
  xnames <- names(xvector)
  
  S <- S1 <- length(snames)
  sindex <- c(1:S)
  knames <- snames
  
  nc <- nrow(bchain)
  
  if(length(omitY) > 0){
    
    wob <- grep(omitY,colnames(bchain))
    wos <- grep(omitY,colnames(schain))
    
    bchain[,wob] <- 0
    schain[,wob] <- 0
    sindex <- sindex[!snames %in% omitY]
    knames <- snames[sindex]
    S1     <- length(knames)
  }
  
  tmp <- as.vector( outer(snames,snames,paste,sep='_') )
  colnames(schain) <- tmp
  
  nspec <- length(snames)
  #indirect effects
  
  bchain <- reverseSignBetaChains(bchain,reverseSign=reverseSign)
  
  #  if(length(reverseSign) > 0){              #reverseSign means 1 - x
  #    for(j in 1:length(reverseSign)){
  #      wj <- grep(reverseSign[j],names(xx))
  #      if(length(wj) > 0){
  #        xx[wj][names(xx)[wj] != 'xeric'] <- 1 - xx[wj]
  #    }
  #  }
  
  ww  <- grep(':',xnames)
  main <- xnames
  if(length(ww) > 0)main <- xnames[-ww]
  main <- main[main != 'intercept']
  int  <- unique( unlist( strsplit(xnames[ww],':') ) ) 
  
  mainEffect <- matrix(NA,nspec,length(main))
  colnames(mainEffect) <- main
  rownames(mainEffect) <- snames
  intEffect  <- dirEffect <- indEffectFrom <- indEffectTo <- mainEffect
  mainSd <- dirSd <- intSd <- indSdFrom <- indSdTo <- mainEffect 
  
  
  gs <- sample(round(nc/5,0):nc,nsim,replace=T)
  
  
  pbar <- txtProgressBar(min=1,max=length(main),style=1)
  
  for(j in 1:length(main)){
    
  #  xj <- paste(main[j],snames,sep='_')
    
    ttt <- interactionsFromGibbs(mainx=main[j],bchain=bchain,schain=schain,
                                 specs=snames,xmnames=names(xvector),xx=xvector,
                                 omitY = omitY)
    maine <- ttt$main
    inter <- ttt$inter
    indirFrom <- indirTo  <- maine*0
    direct <- maine + inter
    
    if(MEAN){
      dmain  <- colMeans(maine)
      inte   <- colMeans(inter)
      dir    <- colMeans(direct)
    } else {
      dmain  <- apply(maine,2,median)
      inte   <- apply(inter,2,median)
      dir    <- apply(direct,2,median)
    }
    
    mainEffect[sindex,j] <- dmain
    intEffect[sindex,j]  <- inte
    dirEffect[sindex,j]  <- dir
    
    mainSd[sindex,j] <- apply(maine,2,sd)
    intSd[sindex,j]  <- apply(inter,2,sd)
    dirSd[sindex,j]  <- apply(direct,2,sd)
    
    kk <- 0
    
    setTxtProgressBar(pbar,j)
    
    for(s in sindex){
      
      # s effect on all others
      
      kk <- kk + 1
      
      wc   <- paste(knames[kk],knames[-kk],sep='_')
      sonk <- schain[,wc]
      wc   <- paste(knames[kk],knames[kk],sep='_')   # var(s)
      vs   <- schain[,wc] 
      
      wi <- paste(main[j],knames[kk],sep='_')
      wd <- which(colnames(direct) == wi)
      if(length(wd) == 0)wi <- paste(knames[kk],main[j],sep='_')
      indirFrom[,-kk] <- indirFrom[,-kk] + sonk/matrix(vs,nc,S1-1)*matrix(direct[,wi],nc,S1-1)
      
      # other effects on s
      
      wi <- paste(main[j],knames[-kk],sep='_')
      wd <- match(wi,colnames(direct))
      cc <- which(is.finite(wd)) 
      if( length(cc)== 0 )wi <- paste(knames[-kk],main[j],sep='_')
      
      ws  <- as.vector( outer(knames,knames,paste,sep='_') )
      
      for(g in gs){
        ss  <- matrix(schain[g,ws],S1,S1)
        indirTo[g,kk] <- matrix(ss[-kk,kk],1)%*%solve(ss[-kk,-kk])%*%direct[g,wi]
      }
    }
    
    if(MEAN){
      indirectFrom <- colMeans(indirFrom)
      indirectTo   <- colMeans(indirTo[gs,])
    }
    if(!MEAN){
      indirectFrom <- apply(indirFrom,2,median)
      indirectTo   <- apply(indirTo[gs,],2,median)
    }
    indEffectFrom[sindex,j] <- indirectFrom
    indSdFrom[sindex,j]     <- apply(indirFrom,2,sd)
    indEffectTo[sindex,j]   <- indirectTo
    indSdTo[sindex,j]       <- apply(indirTo[gs,],2,sd)
  }
    
  
  if(!is.null(keepNames)){
    wk <- which(rownames(mainEffect) %in% keepNames)
    mainEffect <- mainEffect[wk,]
    intEffect <- intEffect[wk,]
    dirEffect <- dirEffect[wk,]
    indEffectFrom <- indEffectFrom[wk,]
    indEffectTo   <- indEffectTo[wk,]
    mainSd    <- mainSd[wk,]
    dirSd     <- dirSd[wk,]
    indSdFrom     <- indSdFrom[wk,]
    indSdTo <- indSdTo[wk,]
  }
  
  if(!is.null(sum2zero)){
    for(j in length(sum2zero)){
      wm <- match(sum2zero[[j]])
      ss <- rowMeans(mainEffect[,wm])
      tmp <- sweep(mainEffect[,wm],1,ss,'-')
      mainEffect[,wm] <- tmp
    }
  }
  
  if(!is.null(standardX)){
    wx <- match(colnames(mainEffect),names(standardX))
    sx <- matrix(standardX[wx],nrow(mainEffect),length(wx),byrow=T)
    mainEffect <- mainEffect*sx
    intEffect <- intEffect*sx
    dirEffect <- dirEffect*sx
    indEffectFrom <- indEffectFrom*sx
    indEffectTo <- indEffectTo*sx
  }
  
  list(mainEffect = mainEffect, intEffect = intEffect, dirEffect = dirEffect,
       indEffectFrom = indEffectFrom, indEffectTo = indEffectTo, mainSd = mainSd, dirSd = dirSd,
       intSd = intSd, indSdFrom = indSdFrom, indSdTo = indSdTo)
}

standardizeColumns <- function(xx){
  
  #center and standardize columns of a matrix
  
  if(!is.matrix(xx)){
    return(xx - mean(xx,na.rm=T))/sd(xx,na.rm=T)
  } else {
    xm <- colMeans(xx,na.rm=T)
    xs <- apply(xx,2,sd,na.rm=T)
    xs <- (xx - matrix(xm,nrow(xx),ncol(xx),byrow=T))/matrix(xs,nrow(xx),ncol(xx),byrow=T)
  }
  xs
}

########################################################
conditionalPred <- function(bchains,schains,x,yy,predCol,
                            partition=NULL,typeNames='CA',
                            nsim=2000,burnin=1,standard=F,omitY = NULL){
  
  # yMu - conditional prediction
  # aMu - joint prediction
  # x   - if a matrix predicts the data, if a vector, predicts the scenario
  # yy  - list of response matrices; if more than one, returns log response ratio
  
  if(!is.list(yy))yy <- list('yy' == yy)
  ny <- length(yy)
  
  S <- ncol(yy[[1]])
  S1 <- S - 1

  
  if(length(omitY) > 0){
    
    wob <- grep(omitY,colnames(bchains))
    wos <- grep(omitY,colnames(schains))
    
  #  bchains[,wob] <- 0
    schains[,wos] <- 0
    
    schains[,paste(omitY,omitY,sep='_')] <- 1
    sindex <- sindex[!snames %in% omitY]
  }
  
  
  if(!is.matrix(x))x <- matrix(x,1)
  for(k in 1:ny)if(!is.matrix(yy[[k]]))yy[[k]] <- matrix(yy[[k]],1)
  
  Q <- ncol(x)
  S <- ncol(yy[[1]])
  n <- nrow(x)
  P <- length(predCol)
  sindex <- c(1:S1)
  
  if(length(typeNames) == 1)typeNames <- rep(typeNames,S)
  
  part <- F
  if(!is.null(partition))part <- T
  
  i <- sample(burnin:nrow(bchains),nsim,replace=T)
  
  logResRatio <- NULL
  
  pmat <- matrix(0,n,P)
  predList <- vector('list',ny)
  if(ny > 1)logResRatio <- vector('list',(ny-1))   # response relative to first yy
  
  for(k in 1:ny){
    
    predList[[k]] <- list(ypred = pmat, ypred2 = pmat, ypall = pmat,ypall2 = pmat, yprob = pmat)
    names(predList[[k]])   <- c('ypred','ypred2','ypall','ypall2','yprob')
    
    if('CC' %in% typeNames){
      yss <- yy[[k]]
      compCols <- which(typeNames == 'CC')
      ysum <- rowSums(yss[,compCols])
      yss[,compCols]   <- sweep(yss[,compCols],1,ysum,'/')
      yy[[k]] <- yss
    }
    
    if(k > 1)logResRatio[[k-1]] <- pmat
    
    predList[[k]]$y <- yy[[k]]
  }
  
  for(j in 1:nsim){
    
    bb <- matrix(bchains[i[j],],Q,S1)
    ss <- matrix(schains[i[j],],S1,S1)
    
    if(standard){
      bb <- bb/matrix(diag(ss),Q,S1,byrow=T)
      ss <- cov2cor(ss)
    }
    
    ss  <- ss[sindex,sindex]
    mmu <- x%*%bb[,sindex]
    
    for(k in 1:ny){
      
      tmp <- conditionalMVNRcpp(yy[[k]][,sindex],mmu,ss,cdex=predCol)
      
      mu <- tmp$mu
      vr <- tmp$vr
      
      if(ncol(mu) == 1)yc <- rnorm(n,mu,sqrt(vr))     #conditional
      if(ncol(mu) > 1) yc <- myrmvnorm(n,mu,vr)
      ym <- myrmvnorm(n,mmu,ss)[,predCol]             #marginal
      if(part){
        yc <- matrix( findInterval(yc,partition) - 1,n,P)
        ym <- matrix( findInterval(ym,partition) - 1,n,P)
      }
      
      yprob <- pnorm(yc,0,1)
      if(k == 1)yprob1 <- yprob
      
      predList[[k]]$ypred  <- predList[[k]]$ypred + yc        #conditional
      predList[[k]]$ypred2 <- predList[[k]]$ypred2 + yc^2
      predList[[k]]$ypall  <- predList[[k]]$ypall + ym        #marginal
      predList[[k]]$yall2  <- predList[[k]]$ypall2 + ym^2
      predList[[k]]$yprob  <- predList[[k]]$yprob + yprob
      
      if(k > 1)logResRatio[[k-1]] <- log(yprob/yprob1)
    }
  }
  
  for(k in 1:ny){
    
    yMu <- predList[[k]]$ypred/nsim
    ySe <- sqrt( predList[[k]]$ypred2/nsim - yMu^2 )
    
    predList[[k]]$condMu <- yMu
    predList[[k]]$condSe <- ySe
    
    aMu <- predList[[k]]$ypall/nsim                       #marginal
    aSe <- sqrt( predList[[k]]$ypall2/nsim - aMu^2 )
    
    predList[[k]]$jointMu <- aMu
    predList[[k]]$jointSe <- aSe
    
    predList[[k]]$condPr <-  predList[[k]]$yprob/nsim
  }
  
  list(predList = predList, logResRatio = logResRatio)
  
}


###########################################
conditionalMVNVec <- function(xx, mu, sigma, cdex, SOLVESTOP=T ){ 

  #  xx, mu are n by p matrices
  #  sigma is p by p or n by p (one row per x)
  #  cdex is a vector of elements < p
  #  if !SOLVESTOP inversion error returns mu[,cdex] and variance

  n <- nrow(xx)
  tiny <- min(diag(sigma))*.0001
  p   <- ncol(mu)
  if(ncol(xx) != p | ncol(sigma) != p)stop('different lengths in conditionalMVNVec')
  
  testv <- try(chol(sigma[-cdex,-cdex]),T)
  if(inherits(testv,'try-error')){
      if(SOLVESTOP) return( list(mu = numeric(0), vr = numeric(0)) )
      if(!SOLVESTOP)return( list(mu = mu[,cdex],vr = var(mu[,cdex],na.rm=T)) )
  }
  sinv  <- chol2inv(testv)
  
 # testv <- try(chol(sigma),T)
 # if(inherits(testv,'try-error')){
 #   if(SOLVESTOP) return( list(mu = numeric(0), vr = numeric(0)) )
 #   if(!SOLVESTOP)return( list(mu = mu[,cdex],vr = var(mu[,cdex],na.rm=T)) )
 # }
 # sinv <- chol2inv(testv)
 # vr1  <- solve(sinv[-p,-p])
  
  
  p1  <- sigma[cdex,-cdex]%*%sinv
  mu1 <- mu[,cdex] + t( p1%*%t(xx[,-cdex] - mu[,-cdex]) )
  vr1 <- sigma[cdex,cdex] - p1%*%sigma[-cdex,cdex]
  
  list(mu = mu1, vr = vr1)
}


#######################################

conditionalMVN <- function(xx, mu, sigma, cindex){  

  # x and mu are vectors, cindex is vector index for conditional

  tiny <- min(diag(sigma))*.0001
  nm   <- length(mu)
  if(length(xx) != nm)stop('x and mu different length in conditionalMVN')

  xx <- matrix(xx,nrow=1)
  mu <- matrix(mu,nrow=1)

  testv <- try(chol(sigma[-cindex,-cindex]),T)
  if(inherits(testv,'try-error')){
      return( list(mu = numeric(0), vr = numeric(0)) )
  }

  sin <- chol2inv(testv)
  p1  <- sigma[cindex,-cindex]%*%sin

  mu1 <- mu[cindex] +  p1%*%(xx[-cindex] - mu[-cindex]) 
  vr1 <- sigma[cindex,cindex] - p1%*%sigma[-cindex,cindex]

  list(mu = mu1, vr = vr1)
}
###################

processStates <- function(xchains,y,ADD=F,cols=1){
	
	xci <- apply(xchains,2,quantile,c(.5,.025,.975))
	
	yr <- range(xci)
	
	if(!ADD)plot(xci[1,],type='l',lwd=2,ylim=yr,col=cols)
	lines(xci[1,],type='l',lwd=2,ylim=yr,col=cols)
	for(j in 2:3)lines(xci[j,],lty=2,col=cols)
	points(y,col=cols[1])

        invisible(xci)
}

##########################################################
processPars <- function(xgb,xtrue=numeric(0),CPLOT=F,DPLOT=F,
                        sigOnly = F,burnin=1,xlimits = NULL){  

  #xg      - matrix of gibbs chains
  #xtrue   - true values (simulated data)
  #CPLOT   - if T, plot chains
  #DPLOT   - if T, plot density
  #burnin  - analyze chains > burnin
  #xlimits - xlimits for plot
  #sigOnly - plot only parameters that 95% CI does not include 0
  
  if(!is.matrix(xgb))xgb <- matrix(xgb,ncol=1)
  if(is.null(colnames(xgb)))colnames(xgb) <- paste('V',c(1:ncol(xgb)),sep='-')
  
  NOPARS <- F
  
  if(sigOnly){
    wi   <- grep('intercept',colnames(xgb))      #extract covariates for plotting
    btmp <- xgb
    if(length(wi) > 0){
    	btmp <- xgb[,-wi]
      if(length(xtrue) > 0)xtrue <- xtrue[-wi]
    }

    wq   <- apply(btmp,2,quantile,c(.025,.975),na.rm=T)  #extract parameters != 0
    wq   <- which(wq[1,] < 0 & wq[2,] > 0)
    
    if(length(wq) == ncol(btmp))NOPARS <- T
    if(NOPARS) warning('no significant pars to plot')
    if(length(wq) > 0 & !NOPARS){
      xgb  <- btmp[,-wq]
      if(length(xtrue) > 0)xtrue <- xtrue[-wq]
    }
   }

  if(!is.matrix(xgb))xgb <- as.matrix(xgb)
  if(burnin > 1){
  	     if(burnin > (nrow(xgb) + 100))stop("burnin too large")
  	     xgb <- xgb[-c(1:burnin),]
  }
  if(!is.matrix(xgb))xgb <- as.matrix(xgb)
  nc <- ncol(xgb)
  nf <- round(sqrt(nc),0)

  out <- t(rbind(apply(xgb,2,mean,na.rm=T),apply(xgb,2,sd,na.rm=T),
                 apply(xgb,2,quantile,c(.025,.975),na.rm=T)))
  if(!is.null(colnames(xgb)))rownames(out) <- colnames(xgb)
  colnames(out) <- c('estimate','se','0.025','0.975')
  if(length(xtrue) > 0){
    out <- cbind(out,xtrue)
    colnames(out) <- c('estimate','se','0.025','0.975','true value')
  }

  if(CPLOT | DPLOT)par(mfrow=c((nf+1),nf),mar=c(4,2,2,2))
  if(CPLOT & DPLOT)par(mfrow=c((nf+1),nc),mar=c(4,2,2,2))

  if(CPLOT & !NOPARS){
      for(j in 1:nc){
       plot(xgb[,j],type='l')
       abline(h=out[j,],lty=2)
       if(length(xtrue) > 0)abline(h=xtrue[j],col='red')
       abline(h = 0, col='grey',lwd=2)
       title(colnames(xgb)[j])
     }
  }
  xlims <- xlimits
  if(DPLOT & !NOPARS){
      for(j in 1:nc){
        xj <- density(xgb[,j])
        if(is.null(xlimits))xlims <- range(xj$x)
        plot(xj$x,xj$y,type='l',xlim=xlims)
        abline(v=out[j,],lty=2)
        if(length(xtrue) > 0)abline(v=xtrue[j],col='red')
        title(colnames(xgb)[j])
     }
  }
  list(summary = signif(out,4)
)

}

#####################
distmat <- function(xt,yt,xs,ys){
    xd <- outer(xt,xs,function(xt,xs) (xt - xs)^2)
    yd <- outer(yt,ys,function(yt,ys) (yt - ys)^2)
    t(sqrt(xd + yd)) 
}
#################################
parSummary <- function(x){

    if(!is.matrix(x)){
      xb <- c(mean(x),sd(x),quantile(x,c(.025,.975)))
    }
    if(is.matrix(x)){
      xb <- cbind(apply(x,2,mean,na.rm=T),apply(x,2,sd,na.rm=T),
            t(apply(x,2,quantile,c(.025,.975),na.rm=T)))
    }
    xb
}


multiLogitStates <- function(b,discStates,contStates,maxD){    # update of continuous states based on discrete states

      tiny <- 1e-20
      nn   <- length(discStates)
 
      discStates[is.na(discStates)] <- 0
      prob <- discStates*0 + 1

      if(maxD == 1)return(log(prob))

      wk <- which(discStates == 1 & is.finite(contStates),arr.ind=T)

      prob[wk] <- invlogt(b[1,],contStates[wk])

      if(maxD > 2){
        for(k in 2:(maxD-1)){
         wk       <- which(discStates == k & is.finite(contStates),arr.ind=T)
         prob[wk] <- invlogt(b[k,],contStates[wk]) - 
                     invlogt(b[(k-1),],contStates[wk])
        }
      }
      prob[discStates == maxD] <- 1 - invlogt(b[(maxD-1),],contStates[discStates == maxD])

    prob[prob < tiny] <- tiny
    log(prob)
}


#############################################

logit2Prob <- function(y){     #multivar logit to fractions

  if(!is.matrix(y)){
     z <- invlogit(y)
     return(cbind(z,1 - z))
  }

  zs   <- rowSums(exp(y))
  z1   <- 1/(1 + zs)
  zm   <- exp(y)/ (1 + zs)
  cbind(zm,z1)
}

prob2Logit <- function(z){     #fractions to multivar logit

  r <- ncol(z)
  if(r == 2){
    ss <- z[,1]/rowSums(z)
    return(log(ss/(1 - ss) ))
  }

  log(z[,-r]/(1 - rowSums(z[,-r])))
  
}

fisherR2Z <- function(r) {
  
  .5*(log(1 + r) - log(1 - r) )

}

compareCorMatrix <- function(r1,r2,n1,n2,alpha=.05){
  
  # zou (2007) Psycholog Methods 12:399
  # ci - r1 significantly less than or greater than r2
  # sumDis - total distance
  # sumDisST - in units of standard deviations
  
  f1 <- fisherR2Z(r1)
  f2 <- fisherR2Z(r2)
  
  zscore <- (f1 - f2)/sqrt( 1/(n1 - 3) + 1/(n2 - 3) )
  
  p <- 2*pnorm(-abs(zscore))
  
  zalpha <- qnorm(alpha/2)
  zd     <- zalpha*sqrt(1/(n1 - 3))
  
  l1 <- f1 + zd
  u1 <- f1 - zd
  l2 <- f2 + zd
  u2 <- f2 - zd
  
  l1 <- (exp(2*l1) - 1)/(exp(2*l1) + 1)
  u1 <- (exp(2*u1) - 1)/(exp(2*u1) + 1)
  l2 <- (exp(2*l2) - 1)/(exp(2*l2) + 1)
  u2 <- (exp(2*u2) - 1)/(exp(2*u2) + 1)
  
  dr <- r1 - r2
  
  lo <- dr - sqrt( (r1 - l1)^2 + (u2 - r2)^2 )
  hi <- dr + sqrt( (u1 - r1)^2 + (r2 - l2)^2 )
  
  dr <- r1 - r2
  
  wl <- which( lo < 0 & hi < 0)
  wh <- which( lo > 0 & hi > 0)
  
  # stdevs from zero
  stdev <- (dr - lo)/1.96
  from0 <- dr/stdev
  
  cmat <- matrix(0,nrow(r1),nrow(r1))
  cmat[wl] <- -1
  cmat[wh] <- 1
  
  loTri    <- upper.tri(cmat)
  sumPos   <- length( which(cmat[loTri] == 1) )
  sumNeg   <- length( which(cmat[loTri] == -1) )
  fraction <- (sumPos + sumNeg)/length(cmat[loTri])
  meanDisST <- sum( abs(from0[loTri]) )/length(cmat[loTri])
  meanDis   <- sum( abs(dr[loTri]) )/length(cmat[loTri])
  
  
  list(ci = cmat, pvalue = p, standDist = from0, sumPos = sumPos, sumNeg = sumNeg,
       fraction = fraction, meanDis = meanDis, meanDisST = meanDisST)
}
  
#################################
clusterPlot <- function(dcor=NULL,dist=NULL,main=' ',xlab='Species',method='complete',
                        cex=1,ncluster=2, add=F,
                        xlim=NULL, colCode = NULL,horiz=T,textSize=1){  
  
  #dcor is a correlation matrix
  #dist is a distance matrix

  
  if(!is.null(dist)){
    nn <- nrow(dist)
    diss <- as.dist( dist )
  }
  if(is.null(dist)){
    nn <- nrow(dcor)
    diss <- as.dist(cov2Dist(dcor))
  }

  tmp    <- hclust(diss,method)
  corder <- tmp$order
  ctmp   <- cutree(tmp,k=1:ncluster)
  
  wclus <- ctmp[,ncluster]
  clusterCol <- NULL
  
  clusterIndex <- ctmp[,ncluster]
  
  clusterList <- character(0)
  
#  mycols <- mapColors(ncluster)
  
  colF   <- colorRampPalette(c('black','blue','orange','brown','red'))
  mycols   <- colF(ncluster)
  
  if(is.null(colCode)){
    colCode <- mycols[ctmp[,ncluster]]
    names(colCode) <- rownames(ctmp)
  }
  
  colLab <- function(n) {
 
    if(is.leaf(n)) {
      a <- attributes(n)
      
      attr(n, "nodePar") <-
        c(a$nodePar, list(col = colCode[n[1]],lab.col = colCode[n[1]]))
    }
    n
  }
  
  tdendro <- as.dendrogram(tmp)
  dL      <- dendrapply(tdendro,colLab)
  
 # textSize <- exp(-.01*nn)
  
  new <- F
  if(add)new <- T
  par(new = new,cex=textSize)
  tmp <- plot(dL,nodePar=list(cex=.1,lab.cex=textSize),horiz=horiz,xlim=xlim)
  title(main)
  
  invisible(list( clusterList = clusterList, colCode = colCode, clusterIndex = clusterIndex,
                  corder = corder) )

}

dendrapply1 <- function (X, FUN, ...) 
{
  FUN <- match.fun(FUN)
  if (!inherits(X, "dendrogram")) 
    stop("'X' is not a dendrogram")
  Napply <- function(d) {
    r <- FUN(d, ...)
    if (!is.leaf(d)) {
      if (!is.list(r)) 
        r <- as.list(r)
      if (length(r) < (n <- length(d))) 
        r[seq_len(n)] <- vector("list", n)
      r[] <- lapply(d, Napply)
    }
    r
  }
  Napply(X)
}



dmvnormAR1 <- function(xx,psi=NULL,const,sinv=NULL,sig=NULL,ev=NULL){  
  
  #const - length nn vector of variance scalars
  #xx    - (x - mu)
  #sinv  - inverse AR1 covariance matrix
  #ev    - eigenvalues of sinv
  #must have either (sinv, ev) or (psi, variance)
  
  nn <- nrow(xx)
  nc <- ncol(xx)
  if(length(const) == 1)const <- rep(const,nn)
  
  if(is.null(sinv)){
    sinv <- invertAR1(nc,psi,sig)
    ev   <- eigen(sinv, only.values = TRUE)$values
  }
  
  s0 <- c(2:(nc-1))
  s1 <- c(1:(nc-2))
  s2 <- c(3:nc)
  
  x1 <- xx[,1]*(xx[,1] - psi*xx[,2])
  xn <- xx[,nc]*(xx[,nc] - psi*xx[,nc-1])
  
  r1 <- 1 + psi^2
  
  xi <- xx[,s0]*(xx[,s0]*r1 - (xx[,s1] + xx[,s2])*psi)
  
  expon <- rowSums( cbind(x1,xi,xn) )/sig
  
  distval <- expon/const
  logdet  <- -nc*log(const) - sum(log(ev))
  
  -(nc * log(2 * pi) + logdet + distval)/2
}

###########################################################


pmvnormAR1_approx <- function(neach,xx,psi,const,sig,lower=-Inf,upper=Inf,stages=ncol(xx)){
  
  # neach - number of samples per row of xx
  # xx    - (x - mu)
  # const - length nn vector of variance scalars
  
  xx <- xx[,stages]    # may not include large sizes, mostly empty
  nc <- length(stages)
  
  nn <- nrow(xx)
  
  oo <- makeAR1(nc,psi,ovar)$varmat
  
  
  mm <- nn*neach
  if(length(lower) == 1)lower <- matrix(lower,mm,nc)
  if(length(upper) == 1)upper <- matrix(upper,mm,nc)
  
  tmp <- myrmvnorm(mm,0,oo)
  ww  <- which(tmp > lower & tmp < upper,arr.ind=T)
  
  tmp <- tmp*0
  tmp[ww] <- 1
  
  tmp <- rowSums(tmp)
  tmp[tmp < nc] <- 0
  tmp[tmp != 0] <- 1
  
  ii <- rep(1:nn,each=neach)
  
  tmp <- byRcpp(tmp,ii,ii*0+1,matrix(0,nn,1),matrix(0,nn,1), fun='sum')/neach
  tmp
}
#######################3 

dmvnormZeroMean <- function(xx,smat=NULL,sinv=NULL){          #MVN density for mean 0
  
  if(!is.matrix(xx))xx <- matrix(xx,1)
  
  if(!is.null(sinv)){
    distval <- diag( xx%*%sinv%*%t(xx) )
    ev      <- eigen(sinv, only.values = TRUE)$values
    logdet  <- -sum(log(ev))
  }
  
  if(is.null(sinv)){
    testv <- try(chol(smat),T)
    if(inherits(testv,'try-error')){
       tiny <- min(abs(xx))/100 + 1e-5
       smat <- smat + diag(diag(smat + tiny))
       testv <- try(chol(smat),T)
    }
    covm <- chol2inv(testv)
    distval <- rowSums((xx %*% covm) * xx)
    ev <- eigen(smat, only.values = TRUE)$values 
    if(min(ev) < 0)ev <- nearPD(smat)$eigenvalues
    logdet   <- sum(log( ev ))
  }

    -(ncol(xx) * log(2 * pi) + logdet + distval)/2
}
########################################

getChainMeanVar <- function(chains){  #mean vector, covariance matrix from chains matrix

  chains <- row2Mat(chains)
  nn     <- nrow(chains)
  if(nn == 1){
    mub    <- chains
    vrb    <- diag(1,ncol(chains))
  }
  if(nn > 1){
     mub   <- colMeans(chains)
     vrb   <- cov(chains)
  }
  list(mu = mub, vr = vrb)
}
#####################################
getPost <- function(b,bchain){     #posterior estimate for marginal likelihood

  #gaussian
  # method of Chib 1995, JASA

 # b <- row2mat(b)

  tmp <- getChainMeanVar(bchain)
  mub <- tmp$mu
  vrb <- tmp$vr
  
  mu  <- matrix(mub,nrow(bchain),ncol(bchain),byrow=T) - bchain
 
  log( mean(exp( dmvnormZeroMean(mu,vrb) )) )

}

#######################################3
biVarMoments <- function(x1,x2,wt=1,PLOT = F, POINTS=F, color=1, pr = .95, ADD=F, lty=1,lwd=1,
                        xylab=c('x','y')){  
  
  require(cluster)

  #x1, x2 - vectors for variables, wt is weight


  if(length(pr) > 1){
     if(length(lty) == 1)lty = rep(lty,length(pr))
     if(length(lwd) == 1)lwd = rep(lwd,length(pr))
     if(length(color) == 1)color = rep(color,length(pr))
  }
  
  if(length(wt) == 1)wt <- rep(wt,length(x1))

  ww <- which(is.finite(x1) & is.finite(x2))

  x1 <- x1[ww]
  x2 <- x2[ww]
  wt <- wt[ww]

  w1 <- x1*wt
  w2 <- x2*wt
  m1 <- sum(w1)/sum(wt)
  m2 <- sum(w2)/sum(wt)

  v1  <- sum(wt*x1^2)/sum(wt) - m1^2
  v2  <- sum(wt*x2^2)/sum(wt) - m2^2
  c   <- sum(wt*(x1 - m1)*(x2 - m2))/sum(wt)

  for(k in 1:length(pr)){
      tmp <- list(loc = c(m1,m2), cov = matrix(c(v1,c,c,v2),2,2), d2 = qchisq(pr[k],1) )
      tmp <- predict.ellipsoid(tmp)
      if(PLOT & !ADD & k == 1)plot(tmp[,1],tmp[,2],type='l',col=color[k],lty=lty[k],
                              lwd=lwd[k],xlab=xylab[1], ylab=xylab[2])
      if(POINTS)points(x1,x2,cex=.3,col='grey')
      if(ADD | k > 1)lines(tmp[,1],tmp[,2],type='l',col=color[k],lwd=lwd[k],lty=lty[k])
  }
  
  invisible( list(loc = c(m1,m2), 
                    cov = matrix(c(v1,c,c,v2),2,2) ,
                    d2 = qchisq(pr[1],1), ellipse = tmp) )
}
#############################################################
pmake <- function(pars){  

  p   <- pars[1]                   #frequency of a
  f   <- pars[2]                   #inbreeding coefficient 
  paa <- p^2 + f*p*(1 - p)	  #frequency of p.aa
  pab <- 2*p*(1 - p)*(1 - f)	  #frequency of p.ab
  pbb <- (1 - p)^2 + f*p*(1 - p)  #frequency of p.bb
  c(paa,pab,pbb)
}

minf <- function(p){   #minimum inbreeding coefficient

  lof <- rep(-1,length(p))
  lop <- -(1 - p)/p
  hip <- -p/(1 - p)
  lof[p > (1 - p)] <- lop[p > (1 - p)]
  lof[p < (1 - p)] <- hip[p < (1 - p)]
  lof
}
pf.update <- function(p,f){  #log posterior

  multlik(c(p,f)) + dbeta(p,g1,g2,log=T) + 
           dnorm(f,priorF,priorFSD,log=T)
}

multlik <- function(pars){  #multinom likelihood
 
  p    <- pmake(pars)
  pmat <- matrix(rep(p,npop),3,npop)
  sum(y*log(pmat))

}

update_pf <- function(){

  propp  <- tnorm(1,.02,.98,pg,.002) #propose pa
  fl     <- minf(propp)              #lower f limit for pa
  propf  <- tnorm(1,fl,1,fg,.01)     #propose f > fl
  pnew   <- pf.update(propp,propf)
  pnow   <- pf.update(pg,fg)
  atmp   <- acceptMH(pnow,pnew,c(pg,fg),c(propp,propf),BLOCK=T)
  pg     <- atmp$x[1]
  fg     <- atmp$x[2]
  ac     <- atmp$accept

  list(pg = pg, fg = fg, accept = ac)
}

updateVariance <- function(yy,mu,s1=1,s2=1,lo=NULL,hi=NULL){
  
  # if yy is a matrix, one variance for each column
  

  k <- 1
  res <- (yy - mu)^2
  nn  <- length(yy)
  
  if(is.matrix(yy)){
    k  <- ncol(yy)
    nn <- nrow(yy)
    sr <- colSums(res)
  }else{
    sr <- sum(res)
  }
  
  u1 <- s1 + nn/2
  u2 <- s2 + .5*sr

  if(is.null(lo) & is.null(hi))return( 1/rgamma(k,u1,u2) )
  
  if(is.null(lo))lo <- 0

  rtrunc(k,lo,hi,u1,u2,'igamma')

}


rtrunc <- function(nn,lo,hi,p1,p2,FN){    #truncated using inv dist sampling
   #truncated 2 parameter distribution
  
  require(pscl)
  
  if(length(lo) == 1)lo <- rep(lo,nn)
  if(length(hi) == 1)hi <- rep(hi,nn)

  pf <- paste('p',FN,sep='')
  qf <- paste('q',FN,sep='')
  pf1 <- match.fun(pf)
  qf1 <- match.fun(qf)
  
  if(FN == 'igamma'){
    z1 <- 1 - pgamma(1/lo, p1, p2)
    z2 <- 1 - pgamma(1/hi, p1, p2)
  }
  if(FN != 'igamma'){
    z1 <- pf1(lo,p1,p2)
    z2 <- pf1(hi,p1,p2)
  }
  
  z  <- runif(nn,z1,z2)
  
  r  <- hi
  wr <- which(z > 0)
  
  if(length(wr) == 0)return(hi)
  
  if(FN == 'igamma')r[wr] <- 1/qgamma(1 - z[wr],p1[wr],p2[wr])
  if(FN != 'igamma')r[wr] <- qf1(z[wr],p1[wr],p2[wr])
  r[r < lo] <- lo[r < lo]
  r
}

fiaCode2spec <- function(fiaCode,genusSpec=F){
  
  # translate FIA codes to our codes
  
  codeTable <- read.table('/nfs/clark/clark.unix/allocationmodel/datafiles/treeCodesDuke.txt',header=T)
  
  mm <- match(fiaCode,codeTable[,'fiaCode'])
  
  nomatch <- which(is.na(mm))
  noCode  <- fiaCode[nomatch]
  
  specCode <- fiaCode
  specCode <- as.character( codeTable[mm,'code'] )
  specCode[is.na(specCode)] <- fiaCode[is.na(specCode)]
  
  if(genusSpec)specCode <- paste(codeTable[mm,'genus'],codeTable[mm,'species'])
  
  #  if(length(nomatch) > 0)specCode[nomatch] <- fiaCode[nomatch]
  
  list(specCode = specCode, noCode = noCode, allCodes = codeTable[,'code'])
}


pasteChar2End <- function(x,char2paste) paste(x,char2paste,sep='') #paste character to end

capFirstLetter <- function(xx) {     #capiltalize first letter of every word
   s <- unlist(strsplit(xx, " "))
   s <- paste(toupper(substring(s, 1, 1)), substring(s, 2),
               sep = "", collapse = " ")
   unlist(strsplit(s, " "))
}

upperFirstLetter <- function(xx) {     #capiltalize first letter of every word
  s <- unlist(strsplit(xx, " "))
  s <- paste(toupper(substring(s, 1, 1)), substring(s, 2),
             sep = "", collapse = " ")
  unlist(strsplit(s, " "))
}

lowerFirstLetter <- function(xx){
  s <- unlist(strsplit(xx, " "))
  s <- paste(tolower(substring(s, 1, 1)), substring(s, 2),
             sep = "", collapse = " ")
  unlist(strsplit(s, " "))
}

joinCharVec <- function(charVec,sep=''){
  
  if(length(charVec) == 1) return(charVec)
  xx <- charVec[1]
  for(j in 2:length(charVec))xx <- paste(xx,charVec[j],sep=sep)
  xx
}

capSubString <- function(x,string) {  #replace substring with caps
   s <- unlist(strsplit(x, string))
   ww <- which(nchar(s) == 0)
   s[ww] <- toupper(x[ww])
   s
}

.replaceString <- function(xx,now='_',new=' '){  #replace now string in vector with new
  
  ww <- grep(now,xx,fixed=T)
  if(length(ww) == 0)return(xx)
  
  for(k in ww){
    s  <- unlist( strsplit(xx[k],now,fixed=T) )
    ss <- s[1]
    if(length(s) == 1)ss <- paste( ss,new,sep='')
    if(length(s) > 1)for(kk in 2:length(s)) ss <- paste( ss,s[kk],sep=new)
    xx[k] <- ss
  }
  xx
}


wideFormat2waterBalance <- function( prec1, pet1, years, ylim = range( c(prec1, pet1) ),
                                     xlab = '', ylab = c( 'P, PET, (mm)'), 
                                     yaxt = rep( 's', nrow(prec1) ) ){
  sites <- rownames(prec1)
  nsite <- length(sites)
  
  if( length( yaxt ) == 1 )yaxt <- rep( yaxt, nsite )
  
  ym <- columnSplit( colnames(prec1), '_' )
  yi <- ym[,1]
  mi <- ym[,2]
  
  ww <- which( yi %in% years )
  prec1 <- prec1[,ww]
  pet1  <- pet1[, ww]
  yi   <- yi[ww]
  mi   <- mi[ww]
  
  siteIndex <- rep( sites, ncol(prec1) )
  monIndex  <- rep(mi, each = nsite)
  siteMonth <- list( site = siteIndex, month = monIndex )
  pr <- round( tapply( as.vector ( prec1 ), siteMonth, mean ), 1 )[sites,]
  pe  <- round( tapply( as.vector ( pet1 ), siteMonth, mean ), 1 )[sites,]
  
  for( i in 1:nsite ){
    plot( NA, xlim = c(0, 13), ylim = ylim, bty = 'n', xaxt = 'n', las = 1,
          ylab = ylab, yaxt = yaxt[i] )
    monthAxis()
    
    shadeThreshold( 1:12, pr[i,], tmin = pe[i,], tmax = NULL, 
                    border = 'darkblue', col = '#9ecae1', LINES = T )
    shadeThreshold( 1:12, pe[i,], tmin = pr[i,], tmax = NULL, 
                    border = '#bd0026', col = '#fed976', LINES = T )
    title( sites[i], cex.main = 1.2, cex.main = 1.5)
  }
}


shadeThreshold <- function( x, y, tmin = NULL, tmax = NULL, ylim = range( c(y, tmin, tmax) ),
                            border = 'brown', col = '#fed976', xaxt = 's', yaxt = 's', 
                            add = T, LINES = F ){
  
  # LINES draws vertical lines at data points
  
  if( length(y) == 1 ){
    if( !is.null( tmin ) )y <- rep( y, length(tmin) )
    if( !is.null( tmax ) )y <- rep( y, length(tmax) )
  }
  if( !is.null( tmin ) )
    if( length(tmin) == 1) tmin <- rep( tmin, length(y) )
  if( !is.null( tmax ) )
    if( length(tmax) == 1) tmax <- rep( tmax, length(y) )
  
  if( !add ){
    plot( x, y, type = 'l', xaxt = 'n', xlab = '', 
          xaxt = xaxt, yaxt = yaxt,
          ylab = '', col = 'white', lwd = 6, ylim = ylim,
          bty = 'n' )
  }
  
  if( !is.null( tmin ) ){   # shade above tmin
    pts <- intersectLines( x, y1 = y, y2 = tmin )
    
    if( length(pts$poly1) > 0 ){
      for( j in 1:length( pts$poly1 ) )
        polygon(pts$poly1[[j]][,1], pts$poly1[[j]][,2], border = border, col = col)
    }
  }
  if( !is.null( tmax ) ){  # shade below tmax
    
    pts <- intersectLines( x, y1 = y, y2 = tmax )
    
    if( length(pts$poly2) > 0 ){
      for( j in 1:length( pts$poly2 ) )
        polygon(pts$poly2[[j]][,1], pts$poly2[[j]][,2], border = border, col = col)
    }
  }
  if( LINES ){
    
    if( is.null(tmax) ){
      ww <- which( y > tmin )
      if( length(ww) > 0 )segments( x[ww], tmin[ww], x[ww], y[ww], col = border )
    }
    if( is.null(tmin) ){
      ww <- which( y < tmax )
      if( length(ww) > 0 )segments( x[ww], y[ww], x[ww], tmax[ww], col = border )
    }
  }
}


intersectLines <- function (x, y1, y2){
  
  # intersection points and polygons for y1 > y2 and y2 > y1
  
  DATE <- F
  
  if( inherits(x[1], 'Date') ){
    # assumes YYYY-MM-DD
    DATE <- T
    date <- x
    xx <- columnSplit( as.character( x ), '-' )
    x  <- as.numeric(xx[,1]) + as.POSIXlt(x)$yday/365
    #   origin <- as.POSIXlt( as.Date(paste( min( xx[,1] ), '-01-01', sep = '' )) 
  }
  
  tiny <- 1e-12
  y1 <- y1 + rnorm( length(y1), 0, tiny )
  
  n <- length(x)
  above <- y1 > y2
  intersectPts <- which(diff(above) != 0) 
  
  y1.diff <- y1[intersectPts+1] - y1[intersectPts]
  y2.diff <- y2[intersectPts+1] - y2[intersectPts]
  x.diff  <- x[intersectPts+1] - x[intersectPts]
  
  slope1 <- y1.diff/x.diff
  slope2 <- y2.diff/x.diff
  intercept1 <- y1[intersectPts] - slope1*x[intersectPts]
  intercept2 <- y2[intersectPts] - slope2*x[intersectPts]
  xPts <- ifelse(slope1 == slope2, NA, 
                 (intercept2-intercept1)/(slope1-slope2))
  yPts <- ifelse(slope1 == slope2, NA,
                 slope1*xPts+intercept1)
  
  jointPts <- which(y1 == y2)
  xPts <- c(xPts, x[jointPts])
  yPts <- c(yPts, y1[jointPts])
  
  ipoints <- cbind( xPts,  yPts)
  
  pt1  <- rbind( c(x[1], y1[1]), ipoints, c(x[n], y1[n]) )
  pt2  <- rbind( c(x[1], y2[1]), ipoints, c(x[n], y2[n]) )
  
  pt2 <- pt2[ !duplicated(pt1[,1]), ]
  pt1 <- pt1[ !duplicated(pt1[,1]), ]
  
  p1 <- p2 <- numeric(0)
  
  for( k in 2:nrow(pt1) ){
    
    w <- which( x >= pt1[k-1,1] & x <= pt1[k,1] )
    if( length(w) == 0 )next
    if( above[w[1]] ){
      pk <- rbind( pt1[k-1,], cbind(x[w], y1[w]), pt1[k,] )
      pb <- rbind( pt2[k-1,], cbind(x[w], y2[w]), pt2[k,] ) #bottom
      pk <- rbind( pk, pb[ nrow(pb):1, ] )
      pk <- pk[ !duplicated(pk[,1:2]), ]
      pk <- rbind( pk, pk[1,] )
      p1 <- append( p1, list( pk ) )
    }else{
      pk <- rbind( pt2[k-1,], cbind(x[w], y2[w]), pt2[k,] ) #top
      pb <- rbind( pt1[k-1,], cbind(x[w], y1[w]), pt1[k,] ) #bottom
      pk <- rbind( pk, pb[ nrow(pb):1, ] )
      pk <- pk[ !duplicated(pk[,1:2]), ]
      pk <- rbind( pk, pk[1,] )
      p2 <- append( p2, list( pk ) )
    }
  }
  if( length(p1) > 0 ){
    for( k in 1:length(p1)){  # all above y2 
      p1[[k]][ abs(p1[[k]][,2]) < 1e-10, 2] <- 0
      pk <- p1[[k]]
      pn   <- RANN::nn2( x, pk[,1], k = 1 )[[1]]
      ymin <- y2[ pn ]
      whi  <- which( pk[, 2] < ymin )
      if( length( whi ) == 0 )next
      pk[ whi, 2] <- ymin[ whi ]
      p1[[k]] <- pk
    }
  }
  if( length(p2) > 0 ){
    for( k in 1:length(p2)){  # all below y2 
      p2[[k]][ abs(p2[[k]][,2]) < 1e-10, 2] <- 0
      pk   <- p2[[k]]
      pn   <- RANN::nn2( x, pk[,1], k = 1 )[[1]]
      ymax <- y2[ pn ]
      whi  <- which( pk[, 2] > ymax )
      if( length( whi ) == 0 )next
      pk[ whi, 2] <- ymax[ whi ]
      p2[[k]] <- pk
    }
  }
  if( length( p1 ) > 1 ){
    p1 <- p1[ which( !sapply( p1, var )[3,] == 0 ) ]
  }else{
    if( length(p1) == 1)p1 <- p1[ which( !sapply( p1, var )[4] == 0 ) ]
  }
  if( length( p1 ) > 1 ){
    p2 <- p2[ which( !sapply( p2, var )[3,] == 0 ) ]
  }else{
    if( length(p2) == 1)p2 <- p2[ which( !sapply( p2, var )[4] == 0 ) ]
  }
  
  list(ipoints = ipoints, poly1 = p1, poly2 = p2) 
}
monthAxis <- function( at = c(1, 4, 7, 10 ), cex.axis = 1.2, tick = T ){
  
  mnames <- c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec' )
  
  axis( 1, at = at, labels = mnames[at], cex.axis = cex.axis, tick = tick ) 
  axis( 1, at = c(1,12), labels = F )
}



whichWords <- function(xx,numWords=1){  #return the first numWords in a character
  
  words <- matrix(character(0),length(xx),numWords)
  
  for(k in 1:length(xx)){
    tmp <- unlist( strsplit(xx[k],' ') )
    if(length(tmp) > numWords)tmp <- tmp[1:numWords]
    words[k,1:length(tmp)] <- tmp
  }
  words
}

joinWords <- function(xx,sep=''){  #join words by rows in xx
  
  if(ncol(xx) == 1)return(xx)
  
  tmp <- xx[,1]
  for(j in 2:ncol(xx))tmp <- paste(tmp,xx[,j],sep=sep)
  tmp
}
  


acceptMH <- function(p0,p1,x0,x1,BLOCK=F){   #accept for M, M-H
	# if BLOCK, then accept as a block,
	# otherwise, accept individually

  nz          <- length(x0)  #no. to accept
  if(BLOCK)nz <- 1
  
  a    <- exp(p1 - p0)       #acceptance PR
  z    <- runif(nz,0,1)
  keep <- which(z < a,arr.ind=T)
  
  if(BLOCK & length(keep) > 0)x0 <- x1
  if(!BLOCK)                  x0[keep] <- x1[keep]           
  ac <- length(keep)        

  list(x = x0, accept = ac)
}

predVsObs <- function(true,p,xlim=range(true),ylim=range(p,na.rm=T),xlab=' ',ylab=' ',
                      colors=rep(1,length(true)),lwd=2,add=F){ 
	
  #true  - length n vector of obs or true values
  #p - ng by n matrix of estimates
  
  if(!is.matrix(p))p <- matrix(p,ncol=1)
  
  nn <- length(true)
  y  <- apply(p,2,quantile,c(.5,.025,.975))
  
  if(!add)plot(true,y[1,],xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,col=colors,pch=3,lwd=lwd)
  points(true,y[1,],col=colors,pch=3,lwd=lwd)

  for(j in 1:nn)lines(c(true[j],true[j]),y[2:3,j],col=colors[j],lwd=lwd)
  abline(0,1,lty=2)
  
  invisible(y)
}



profilePlot <- function(x,y,proValues=quantile(x,seq(.1,.9,by=.2)),ylim=range(y),xlab=' ',ylab=' ',yscale=1){

  plot(x,y,xlab=xlab,ylab=ylab,ylim=ylim,cex=.3,bty='n',col='grey')

  nv <- length(proValues)


  wf <- which(is.finite(x) & is.finite(y))
  x  <- x[wf]
  y  <- y[wf]

  mids <- c(proValues[-nv] + diff(proValues)/2,max(x))
  mids <- c(proValues[1] - diff(proValues[1:2])/2,mids)

  tiny <- 1e-10
  huge <- 1 - tiny

  y[y < tiny] <- tiny
  y[y > huge] <- huge

  fi   <- findInterval(x,mids)
  fi[fi > nv] <- nv

  for(j in 1:nv){

    wj  <- which(fi == j)

    tmp <- hist(y[wj],breaks=quantile(y[wj],seq(0,1,length=20)) ,plot=F)
    yy  <- tmp$density
    yy  <- yy * yscale
    yy[!is.finite(yy)] <- 0
    yy  <- yy  + proValues[j]
    xx  <- tmp$mids

    yy  <- c(proValues[j],yy,proValues[j])
    xx  <- c(xx[1],xx,xx[1])
    
    lines(yy,xx,type='s',lwd=2)
  }
 # points(x,y,cex=.3,col='grey')
 
}


plotLabel <- function(label,location='topleft',cex=1.3,font=1,above=F,
                      below=F,bg=NULL, adj=0, line = 1){
  
  if(above){
    if(location == 'topright' & adj == 0)adj=1
    
    title(label, adj=adj, font.main = font, font.lab =font, cex=cex, 
          line = line, cex.main=cex)
    return()
  }
  if(below){
    if(location == 'bottomright' & adj == 0)adj=1
    mtext(label,side=1,adj=adj, outer=F,font.main = font, 
          font.lab =font,cex=cex)
    return()
  }
    
  if(is.null(bg)){
    tmp <- legend(location,legend=' ',bty='n')
  } else {
    tmp <- legend(location,legend=label,bg=bg,border=bg,text.col=bg,bty='o')
  }
  
  xt <- tmp$rect$left # + tmp$rect$w
  yt <- tmp$text$y
  
  pos <- 4
  tmp <- grep('right',location)
  if(length(tmp) > 0)pos <- 2
  
  XX <- par()$xlog
  YY <- par()$ylog
  
  if(XX)xt <- 10^xt
  if(YY)yt <- 10^yt
  
  text(xt,yt,label,cex=cex,font,pos=pos)
  
}


standardBetaChain <- function(bchain,schain,S,Q,nsim=2000,SENS=F,index=c(1:nrow(bchain))){
  
  #divide betas by sd
  #if SENS return sensitivity b%*%sigma%*%t(b)
  
  bnew <- bchain[1:nsim,]*0
  if(SENS) bnew <- matrix(NA,nsim,Q*Q)
  
  jj <- sample(index,replace=T)
  
  for(j in 1:nsim){
    
    ss <- matrix(schain[jj[j],],S,S)
    sd <- matrix( sqrt(diag(ss)),Q,S,byrow=T)
    bb <- matrix( bchain[jj[j],],Q,S)
    
    bs <- bb/sd
    
    if(!SENS)bnew[j,] <- bs
    if(SENS) bnew[j,] <- as.vector( bs%*%solve(ss)%*%t(bs) )
  }
  bnew
}

plot2wayCI <- function(mu1,mu2,ci1=NULL,ci2=NULL,se1=NULL,se2=NULL,xlab=' ',ylab=' ',
                          xlim=NULL,ylim=NULL,lwd=1,cex=1,LOG=F,barcol='black'){
  
  #mu1, mu2 - vectors
  #ci1, ci2 - n by 2 matrices
  
  if(is.null(ci1)){                    # if no CI, then +/- 1 se
    ci1 <- cbind(mu1 - se1,mu1 + se1)
    ci2 <- cbind(mu2 - se2,mu2 + se2)
  }
  
  if(is.null(xlim))xlim <- range(ci1)
  if(is.null(ylim))ylim <- range(ci2)
  
  if(LOG & xlim[1] < 0)xlim[1] <- .001
  if(LOG & ylim[1] < 0)ylim[1] <- .001
  
  ci1[ci1 < xlim[1]] <- xlim[1]
  ci2[ci2 < ylim[1]] <- ylim[1]
  
  nc <- length(mu1)
  
  if(!LOG)plot(mu1,mu2,xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim,cex=cex)
  if(LOG) plot(mu1,mu2,xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim,cex=cex,log='xy')
  
  ci1[ci1 < xlim[1]] <- xlim[1]
  ci2[ci2 < ylim[1]] <- ylim[1]
  
  for(j in 1:nc){
    lines( rep(mu1[j],2),ci2[j,],lwd=lwd,col=barcol)
    lines( ci1[j,],rep(mu2[j],2),lwd=lwd,col=barcol)
  }
  points(mu1,mu2, cex=cex)
  invisible(cbind(mu1,mu2))
}

plot2wayChain <- function(chain1,chain2,q=c(.5,.025,.975),xlab=' ',ylab=' ',
                       xlim=NULL,ylim=NULL,lwd=1,cex=1,barcol='black'){
  
  q1 <- apply(chain1,2,quantile,q,na.rm=T)
  q2 <- apply(chain2,2,quantile,q,na.tm=T)
  
  tmp <- plot2wayCI(mu1=q1[1,],mu2=q2[1,],ci1=t(q1[2:3,]),ci2=t(q2[2:3,]),
             xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim,lwd=lwd,cex=cex,barcol=barcol)
  invisible(tmp)
}


  

#######################################################
plotObsPred <- function(obs,yMean,ySE=NULL,nbin=NULL,nPerBin=NULL,breaks=NULL,
                        LOG=F,xlimit=NULL,ylimit=NULL,xlabel='Observed',ylabel='Predicted',
                        ptcol=NULL,boxPerc = .6826895, whiskerPerc = .95,
                        fill=NULL,add=F,box.col='black',wide=NULL,POINTS=T,
                        MEDIAN=T){
  
  aa <- (1 - boxPerc)/2
  boxQuant <- c(aa, 1 - aa )
  aa <- (1 - whiskerPerc)/2
  whiskerQuant <- c(aa, 1 - aa )
  
  if(is.null(ptcol)){
    ptcol <- 'black'
    ptcol <- 'grey'
    if(!is.null(nbin))ptcol <- 'grey'
  }
  if(length(ptcol) == 1)ptcol <- rep(ptcol,length(obs))
  
  if(is.null(xlimit))xlimit <- range(obs[is.finite(obs)],na.rm=T)
  if(is.null(ylimit))ylimit <- range(yMean[is.finite(yMean)],na.rm=T)
  
  xxx <- obs
  yyy <- yMean
  
  if(LOG){
    if(is.null(xlimit))xlimit <- range( obs[obs > 0],na.rm=T )
    if(is.null(ylimit))ylimit <- range( yMean[yMean > 0],na.rm=T )
    if(xlimit[1] <= 0)xlimit[1] <- .001
  }
  
  if(!POINTS){
    xxx <- xlimit[1]
    yyy <- ylimit[1]
  }
  
  if(!add){
    if(is.null(ylimit)){
      if(!LOG)plot(xxx,yyy,col=ptcol,cex=.3,xlab=xlabel,ylab=ylabel)
      if(LOG) plot(xxx,yyy,col=ptcol,cex=.3,xlab=xlabel,ylab=ylabel,log='xy')
    }
    if(!is.null(ylimit)){
      if(!LOG)plot(xxx,yyy,col=ptcol,cex=.3,xlab=xlabel,ylab=ylabel,xlim=xlimit,ylim=ylimit)
      if(LOG) plot(xxx,yyy,col=ptcol,cex=.3,xlab=xlabel,ylab=ylabel,xlim=xlimit,log='xy',ylim=ylimit)
    }
  }
  if(!is.null(ySE)){
    ylo <- yMean - 1.96*ySE
    yhi <- yMean + 1.96*ySE
    for(i in 1:length(obs))lines(c(obs[i],obs[i]),c(ylo[i],yhi[i]),col='grey',lwd=2)
  }
    
  if(is.null(breaks)){
    
    if(is.null(nbin))nbin <- 20
    
    br   <- range(obs[is.finite(obs)],na.rm=T)
    bins <- seq(br[1],br[2],length=nbin)
    if(LOG){
      oo <- min( obs[obs > 0],na.rm=T )
      br[1] <- .5*oo
      bins <- 10^seq(log10(br[1]),log10(br[2]),length=nbin)
    }
    if(!is.null(nPerBin)){
      nbb <- nPerBin/length(obs)
      nbb <- seq(0,1,by=nbb)
      if(max(nbb) < 1)nbb <- c(nbb,1)
      bins <- quantile(obs,nbb,na.rm=T)
      bins <- sort(unique(bins))
      nbin <- length(bins)
    }
  } else {
    bins <- breaks
    nbin <- length(bins)
  }
    
    if(is.null(wide))wide <- diff(bins)/2
    if(length(wide) == 1)wide <- rep(wide,nbin)
  
    for(k in 1:(nbin-1)){
      
      ok <- which(obs >= bins[k] & obs <= bins[k+1])
      qk <- which(is.finite(yMean) & obs >= bins[k] & obs <= bins[k+1])
      q  <- quantile(yMean[qk],c(.5,whiskerQuant[1],boxQuant[1],
                                 boxQuant[2],whiskerQuant[2]),na.rm=T)
      if(LOG){
        q[q <= 0] <- .001
      }
      ym <- mean(yMean[qk])
      xx <- mean(bins[k:(k+1)])
      if(!is.null(nPerBin) & MEDIAN)xx <- median(obs[ok],na.rm=T)
      points(xx,q[1],pch=3,col=box.col)
      yy <- q[c(2,5)]
      yy[1] <- max(yy[1],ylimit[1],na.rm=T) + .0000001
      yy[2] <- max(yy)
      lines(c(xx,xx),yy,lwd=2,col=box.col)
   
      yy1 <- q[3]
      yy1 <- max(yy1,ylimit[1],na.rm=T) + .00000001
      yy2 <- max(yy1,q[4])
      if(is.null(nPerBin)){
        
        minmax <- par('usr')[1:2]
        minx   <- max(c(minmax[1],xx-.5*wide[k])) + .00001
        maxx   <- min(c(minmax[2],xx+.5*wide[k])) - .00001
        rect(minx,yy1,maxx,yy2,col=fill,border=box.col)
        lines(c(minx,maxx),c(ym,ym),lwd=2,col=box.col)
      }
      if(!is.null(nPerBin)){
        qo <- quantile(obs[ok],c(.3,.7,.25,.75),na.rm=T)
        if(!MEDIAN)qo <- c(xx-.2*wide[k],xx+.2*wide[k],xx-.3*wide[k],xx+.3*wide[k])
        rect(qo[1],yy1,qo[2],yy2,col=fill,border=box.col)
        lines(c(qo[3],qo[4]),c(ym,ym),lwd=2,col=box.col)
      }
    }

  invisible( bins )
}


pars2List <- function(pars){    #pars is a named list, filled with values here
  
  for(k in 1:length(pars))pars[[k]] <- get(names(pars)[k])
  pars
}

#######################################################
deviance <- function(y,x,b,s=0,LIKE='norm'){
	    
    if(LIKE == 'norm')  dv <- dnorm(y,x%*%b,sqrt(s),log=T) 
    if(LIKE == 'pois')  dv <- dpois(y,exp(x%*%b),log=T)
    if(LIKE == 'binom') dv <- dbinom(y,1,invlogit(x%*%b),log=T)
    if(LIKE == 'mvnorm')dv <- dmvnormZeroMean(y - x%*%b,s)
    if(LIKE == 'multinom')dv <- multinomLike(y,x,b)$like
    -2*dv
}
#####################################33
dinvGamma <- function(x,a,b,log=FALSE){
	
	p <- a*log(b) - lgamma(a) - (a+1)*log(x) - b/x
	if(log)return(p)
	exp(p)
}
###########################################
dbivarNormFromCols <- function(yy,mu,S){    #bivariate norm, covariances for each y in columns, log density

  #yy  - n by 2
  #mu - n by 2 or scalar 0
  #S  - n by 4: S11, S12, S21, S22

  if(!is.matrix(yy)){
    yy <- matrix(yy,1)
    mu <- matrix(mu,1)
    S  <- matrix(S,1)
  }

  n <- nrow(yy)
  if(length(mu) == 0)mu <- yy*0

  yy <- yy - mu

  invS <- invertcol2(S)
  z    <- yy[,1]^2*(invS[,1] + invS[,2])+ yy[,2]^2*(invS[,3] + invS[,4])

  ldet <- log(S[,1]*S[,4] - S[,2]*S[,3])

  -(n*log(2*pi) + ldet + z)/2
}

###########################################
dtrivarNormFromCols <- function(yy,mu,S){    #trivariate norm, covariances for each y in columns, log density
  
  #y  - n by 3
  #mu - n by 3 or scalar 0
  #S  - n by 9: S11, S12, S31, S21, S22, ...
  
  n <- nrow(y)
  if(length(mu) == 0)mu <- y*0
  
  yy <- yy - mu
  
  tmp <- invertcol3(S)
  invS <- tmp$I
  D    <- tmp$D
  
  y1 <- cbind( rowSums(yy*invS[,1:3]), rowSums(yy*invS[,4:6]), rowSums(yy*invS[,7:9]) )
  z  <- rowSums(y1*yy)
  
  -(n*log(2*pi) + log(D) + z)/2
  
}


#############################################################

invertcol2 <- function(S){   
    # inverse of matrix supplied and returned as nX4 matrix: var1, cov12, cov12, var2

  dt  <- S[,1]*S[,4] - S[,2]*S[,2]
  I1  <- S[,4]
  I2  <- S[,1]
  I12 <- -S[,2]
  I21 <- -S[,2]

  cbind(I1,I12,I21,I2)/dt
}

#############################################################

invertcol3 <- function(S){   
  # inverse of matrix supplied and returned as nX9 matrix: 
  # var1, cov12, cov13, cov21, var2, cov31, cov31, cov32, var3
  
  dt  <- S[,1]*S[,5]*S[,9] + S[,2]*S[,6]*S[,7] + S[,3]*S[,4]*S[,8] - 
        (S[,3]*S[,5]*S[,7] + S[,1]*S[,6]*S[,8] + S[,2]*S[,4]*S[,9])
    
  I1  <- S[,5]*S[,9] - S[,6]*S[,6]
  I2  <- S[,1]*S[,9] - S[,3]*S[,3]
  I3  <- S[,1]*S[,5] - S[,2]*S[,2]
  I12 <- S[,3]*S[,6] - S[,2]*S[,9]
  I13 <- S[,2]*S[,6] - S[,3]*S[,5]
  I23 <- S[,2]*S[,3] - S[,1]*S[,6]
  
  list(I = cbind(I1,I12,I13,I12,I2,I23,I13,I23,I3)/dt, D = dt)
}



############################################################

initialStatesSS <- function(y){

 if(!is.matrix(y))y <- matrix(y,ncol=1)
  r <- ncol(y)

  n    <- length(y)
  time <- c(1:n)
  wm   <- which(is.na(y))
  notMiss <- c(1:n)
  
  x <- y

  if(length(wm) > 0){
   notMiss <- notMiss[-wm]
   x[wm]   <- predict(smooth.spline(time[-wm],y[-wm]),time)$y[wm]
  }

  list(x = as.vector(x), notMiss = notMiss, miss = wm)
}


############################################################
condProb <- function(mu,V,y){     #bivariate conditional on y for 2nd variables

  cm <- mu[1] + V[1,2]/V[2,2]*(y - mu[2])
  cV <- V[1,1] - V[1,2]/V[2,2]*V[2,1]

  list(muY = cm, cV = cV)
}
###########################################

multiLogitStates <- function(b,discStates,contStates,maxD){    # continuous based on discrete states

      tiny <- 1e-20
      nn   <- length(discStates)
 
      discStates[is.na(discStates)] <- 0
      prob <- discStates*0 + 1

      if(maxD == 1)return(log(prob))

      wk <- which(discStates == 1 & is.finite(contStates),arr.ind=T)

      prob[wk] <- invlogt(b[1,],contStates[wk])

      if(maxD > 2){
        for(k in 2:(maxD-1)){
         wk <- which(discStates == k & is.finite(contStates),arr.ind=T)
         prob[wk] <- invlogt(b[k,],contStates[wk]) - 
                                  invlogt(b[(k-1),],contStates[wk])
        }
      }
      prob[discStates == maxD] <- 1 - invlogt(b[(maxD-1),],contStates[discStates == maxD])

    prob[prob < tiny] <- tiny
    log(prob)
}

############################################
breaks2pars <- function(cmat,breaks){        #break points and coefficents to ordinal multinomial logit pars

  nd <- length(breaks)
  c0 <- rep(0,nrow(cmat))

  for(k in 1:nd){

    qk <- 0
    if(k > 1){
      for(j in 1:(k-1)){
       ej <- exp(c0[j] + cmat[j,2]*breaks[k])
       qk <- qk + ej/(1 + ej)
      } 
    }
    D     <- -log(1/(.5 + qk) - 1)
    cc0   <- D - cmat[k,2]*breaks[k]
    if(k > 1){
      if(cc0 < cmat[k-1,1]){
        cc0 <- cmat[k-1,1] + .5
        cmat[k,2] <- (D - cc0)/breaks[k]
      }
    }
    cmat[k,1] <- cc0
  }
  cmat
}
#################################################
proposeBreak <- function(cmat,br,brLims){      #propose new breakpoints given limits brLims

  c1 <- cmat[,1]
  c2 <- cmat[,2]

  nd <- nrow(cmat)
  brNew <- tnorm(nd,brLims[1,],brLims[2,],br,.1)

  if(nd == 1){
     c1New <- tnorm(1,1,100,cmat[1],.1)
     c2New <- (log(.5) + log(1 + exp(c1New + cmat[2]*brNew)) - c1New)/brNew
     cnew <- matrix(c(c1New,c2New),1,2)
  }

  if(nd > 1)cnew <- getCmat(cmat,brNew)

  list(cmat = cmat,br = brNew)
}
##########################################
pars2p <- function(cmat,h){               #multinomial logit pars and scale h to Pr for multinom logit

  nd <- nrow(cmat)
  nh <- length(h)

  c1 <- matrix(cmat[,1],nh,nd,byrow=T)
  c2 <- matrix(cmat[,2],nh,nd,byrow=T)
  hh <- matrix(h,nh,nd)
  eh <- exp(c1 + c2*hh)

  theta <- matrix(0,nh,nd+1)
  sumt  <- rep(0,nh)

  for(k in 1:nd){

     tk   <- eh[,k]/(1 + eh[,k]) - sumt
     theta[,k] <- tk 
     sumt <- sumt + tk
  }
  theta[,nd+1] <- 1 - rowSums(theta)

  theta
}


######################################################
plotLogit <- function(lims,breaks,h,cmat){  #plot multinomial ordinal logit

  tmp  <- pars2p(cmat,h)

  plot(h,tmp[,1],type='l',lwd=2,xlab='latent health scale',ylab='Probabilty')

  for(j in 1:ncol(lims)){
 #   polygon(c(lims[,j],rev(lims[,j])),c(0,0,1,1),border=NA,col='grey')
    lines(h,tmp[,j],col=j,lwd=2)
    l1   <- 0
    if(j > 1)l1 <- breaks[j-1]
    midx <- (l1 + breaks[j])/2
    text(midx,.8,j,col=j,cex=1.4)
  }
  lines(h,tmp[,j+1],col=j+1,lwd=2)
  midx <- ( breaks[j] + 100 )/2
  text(midx,.8,j+1,col=j+1,cex=1.4)
  abline(h=.5,lty=2)
}   	
###########################

plotSetup <- function(xtic,ytic,xvals = xtic, yvals = ytic, xlabel=' ',ylabel=' ',
                      endFactor=c(.05,.05),
                      fc = NULL,lcHor=NULL,lcVer=NULL){

  xr <- range(xtic)
  yr <- range(ytic)

  xlimit <- c(xr[1] - diff(xr)*endFactor[1],xr[2] + diff(xr)*endFactor[1])
  ylimit <- c(yr[1] - diff(yr)*endFactor[2],yr[2] + diff(yr)*endFactor[2])

  plot(-100,-100,xlim=xlimit,ylim=ylimit,ylab=ylabel,
       cex=1.2,xlab=xlabel,xaxt='n',yaxt='n')

  if(!is.null(fc))rect(xlimit[1],ylimit[1],xlimit[2],ylimit[2],col=fc,border=NA)

  axis(1,at=xtic,labels=xvals,cex.lab=1.2)
  axis(2,at=ytic,labels=yvals,cex.lab=1.2)

  if(!is.null(lcVer))abline(v=xtic,lwd=3,col=lcVer)
  if(!is.null(lcHor))abline(h=ytic,lwd=3,col=lcHor)
}


############################### set up scale for equal dimensions
mapSetup <- function(xlim,ylim,scale=NULL,widex=10.5,widey=6.5){  
  
  #scale is x per inch
  
  if(is.null(scale))scale <- 1

  px   <- diff(xlim)/scale
  py   <- diff(ylim)/scale
  
  if(px > widex){
    dx <- widex/px
    px <- widex
    py <- py*dx
  }
  if(py > widey){
    dx <- widey/py
    py <- widey
    px <- px*dx
  }
    
  par(pin=c(px,py))
  
  invisible( c(px,py) )

}


getQuant <- function(data,dim,q){

  quantile(data,dim,q,na.rm=T)
}

cov2Dist <- function(sigma){ #distance induced by covariance
	
	n <- nrow(sigma)
	matrix(diag(sigma),n,n) + matrix(diag(sigma),n,n,byrow=T) - 2*sigma
}


myBarPlot <- function(dataMat,widthFactor=1,xlab=' ',ylab=' ',ylim=NULL,barLabels=NULL,
                      textFactor=1,stack=F){
  
  #dataMat - one row per bar, columns are bar segments
  
  nc <- ncol(dataMat)
  nr <- nrow(dataMat)
  
  if(is.null(ylim))ylim=c(range(rowSums(dataMat,na.rm=T)))
  
  wide <- widthFactor
  
  plot(c(0,0),c(0,0),xlim=c(0,nr),cex=.1,ylim=ylim,xlab = xlab,ylab=ylab)
  
  ystart <- numeric(0)
  
  if(!stack){
    rr    <- range(dataMat)
    space <- diff(rr)*.2
    tmp   <- apply(dataMat,2,max,na.rm=T)
    ystart <- cumsum( tmp + space )
  }
  
  for(j in 1:nr){
    
    xj <- c(j - wide,j + wide)
    
    y1 <- 0
    
    for(k in 1:nc){
      
      y2 <- y1 + dataMat[j,k]
      rect(xj[1],y1,xj[2],y2,col=k,border=k)

      y1 <- y2
      if(!stack)y1 <- ystart[k]
    }
    if(!is.null(barLabels))text(xj[1],y2+.1*diff(range(ylim)),
                                barLabels[j],srt=90,pos=4,cex=textFactor)
  }
  
  invisible(ystart)
}

####################################################

inBox <- function(xx,bounds=matrix(c(-Inf,Inf),2,ncol(x))){
  
  # find rows of x within bounds for each column
  
  ww <- numeric(0)
  
  for(j in 1:ncol(bounds)){
    
    ww   <- c(ww,which(xx[,j] < bounds[1,j] | xx[,j] > bounds[2,j]))
  }
  
  ww <- unique(ww)
  
  keep <- c(1:nrow(xx))[-ww]
  
  list(x = xx[keep,], keep = keep)
}


stackedBoxPlot <- function( stackList,stackSd=character(0),
                            ylim=NULL,sortBy = NULL, barnames=NULL,
                            col=rep(NULL,length(stackList)),
                            border=rep(NA,length(stackList)),
                            decreasing=T, nsd=1.96, cex=1,
                            legend=NULL,scaleLegend=.1){
  
  # sortBy - if length 1 indicates which variable in stackList to sort by
  #        - if a vector it is the order to plot
  # nds    - no. standard deviations for whiskers
  
  nn  <- length(stackList)
  ord <- c(1:length(stackList[[1]]))
  
  xx <- 0:(length(ord)-1)
  
  if(is.null(ylim)){
    
    ymax <- ymin <- 0
    
    for(j in 1:nn){
      ymax <- ymax + max( c(0,stackList[[j]]) )
      ymin <- ymin + min( c(0,stackList[[j]]) )
    }
      
    ylim <- c(ymin,ymax)
    
    yscale <- diff(ylim)*.4
    ylim[1] <- ylim[1] - yscale
    ylim[2] <- ylim[2] + yscale
  }
  
  
  if(!is.null(sortBy)){
    
    if(length(sortBy) > 1){
      ord <- sortBy
    } else {
      ord <- order( stackList[[sortBy]], decreasing = decreasing)
    }
    if(!is.null(barnames))barnames <- barnames[ord]
  }
  
  add <- F
  
  offset <- offsetPos <- offsetNeg <- rep(0,length(stackList[[1]]))
  
  if(is.null(col))col <- c(1:nn)
  
  for(j in 1:nn){
    
    xj <- stackList[[j]][ord]
    names(xj) <- NULL
    
    wp <- which(xj > 0)
    wn <- which(xj < 0)
    
    offset[wp] <- offsetPos[wp]
    offset[wn] <- offsetNeg[wn]
    
    hj <- xj 
    
    barplot(height= hj,offset=offset,xlim=c(0,1.2*length(ord)),ylim=ylim,
            col=col[j],border=border[j],add=add)
    
    ww <- grep(names(stackList)[j],names(stackSd))
    if(length(ww) > 0){
      xz <- xx + .5
      xz <- xz*1.2
      
      tall <-  nsd*stackSd[[ww]]
      y1   <-  hj + offset + tall
      y2   <-  hj + offset - tall
      
      for(i in 1:length(ord)){
        lines(xz[c(i,i)],c(y1[i],y2[i]),lwd=6,col='white')
        lines(c(xz[i]-.1,xz[i]+.1),y1[c(i,i)],lwd=6,col='white')
        lines(c(xz[i]-.1,xz[i]+.1),y2[c(i,i)],lwd=6,col='white')
        
        lines(xz[c(i,i)],c(y1[i],y2[i]),lwd=2,col=col[j])
        lines(c(xz[i]-.1,xz[i]+.1),y1[c(i,i)],lwd=2,col=col[j])
        lines(c(xz[i]-.1,xz[i]+.1),y2[c(i,i)],lwd=2,col=col[j])
      }
    }
      
      if(j == 1)add <- T
    
    offsetPos[wp] <- offsetPos[wp] + hj[wp]
    offsetNeg[wn] <- offsetNeg[wn] + hj[wn]
    
    if(j == nn & !is.null(barnames)){
      xx[wp] <- xx[wp] + .3
      xx     <- xx*1.2 + .5
      yy     <- offsetNeg
      yy[wn] <- offsetPos[wn]
      yy     <- yy + diff(ylim)*.1
      yy[wp] <- yy[wp] - diff(ylim)*.15
      pos    <- yy*0 + 2
      pos[wn] <- 4
      text(xx,yy,barnames,srt=90.,pos=pos,cex=cex)
    }
  } 
  
  if(!is.null(legend)){
    
    dy <- diff(ylim)*scaleLegend
    dx <- 1.2
    x1 <- length(ord)*.02 + 1
    y1 <- ylim[1]
    pos <- 4
    if(legend == 'topright'){
      x1  <- length(ord)
      y1  <- ylim[2]
      dy  <- -dy
      dx <- -dx
      pos <- 2
    }
    if(legend == 'topleft'){
      y1  <- ylim[2]
      dy  <- -dy
    }
    if(legend == 'bottomright'){
      x1  <- length(ord)
      dx <- -dx
      pos <- 2
    }
    for(j in 1:length(stackList)){
      y2 <- y1 + dy
      rect(x1,y1,x1 + 1,y2,col=col[j],border=border[j])
      y1 <- y2
      colj <- col[j]
      if(colj == 'white')colj <- border[j]
      text(x1 + dx,y1 - dy/2,names(stackList)[[j]],col=colj,pos=pos,cex=cex)
    }
  }
  
  invisible( ord )
}  
  

myBoxPlot <- function(mu,limits,ORDER = F,myorder=NULL,xvalues=c(1:length(mu)),
                      boxcol='black',xlab=' ',ylab=' ',
                      xtic=c(1,length(mu)),ytic=NULL,xvals=xtic,yvals=ytic,
                      plabels=NULL,
                      shadelo=NULL,shadehi=NULL,add=F,widthScale=1){

  #mu has median
  #limits - a list with lo,hi as paired columns, inner to outer
  #shadehi - shade above this value
  #shadelo - shade below this value

  n  <- length(mu)
  xi <- xvalues
  or <- c(1:n)
  if(ORDER & is.null(myorder))or <- order(mu,decreasing=T)
  if(!is.null(myorder))or <- myorder
  
  if(length(boxcol) == 1)boxcol <- rep(boxcol,n)

  if(!is.list(limits))limits <- list(limits = limits)
  x  <- limits
  nq <- length(x)

  if(is.null(ytic))ytic <- signif(range(mu),1)

  if(!add){
    plotSetup(xtic,ytic,xvals=xvals,yvals=yvals,xlabel=xlab,ylabel=ylab)
  }

  if(!is.null(shadehi))rect(xtic[1],shadehi,max(xtic),max(ytic),col='bisque1',border=NA)
  if(!is.null(shadelo))rect(xtic[1],min(ytic),max(xtic),shadelo,col='bisque1',border=NA)

  points(xi,mu[or],pch=3,col=boxcol)

 # ji <- c(1:nq)/diff(range(xi))*3

  ji <- .2*diff(xi)
  ji <- c(ji[1],ji)

  jk <- seq(.4,.6,length=nq)

  for(i in 1:n){

    maxi <- ytic[1]

    for(k in 1:nq){

      yi <- x[[k]][or[i],]
      if(is.na(yi[1]))next
      
      ki <- jk[nq - k + 1]*ji[i]

      x1 <- xi[i] - ki*widthScale
      x2 <- xi[i] + ki*widthScale

      rect(x1,yi[1],x2,yi[2],lwd=1,col=boxcol[i],border=boxcol[i])
      if(k == nq)lines(c(xi[i],xi[i]),yi,lwd=1,col=boxcol[i])
      if(max(yi) > maxi)maxi <- max(yi)
    }

    if(maxi > max(ytic))maxi <- max(ytic) - .1*diff(range(ytic))
    if(!is.null(plabels))text(i-.4,maxi,plabels[or[i]],srt=90,pos=4,col=boxcol[i])
  }

  invisible(or)  #return order

}

trim <- function(xx,p){     #p - lower,upper %tile

  xq <- quantile(xx,p,na.rm=T)
  
  xx[xx < xq[1]] <- xq[1]
  xx[xx > xq[2]] <- xq[2]

  xx
}

ordTraitsFromWts <- function(yWt,ordTraits){
  
  # yWt - n by S species weights
  # ordTraits - S by p ordinal traits
  # returns n by p modal ordinal values
  
  if(!is.matrix(ordTraits))ordTraits <- matrix(ordTraits)
  
  n <- nrow(yWt)
  s <- ncol(yWt)
  
  ii <- rep(c(1:n),s)
  omat <- matrix(NA,n,ncol(ordTraits))
  
  for(j in 1:ncol(ordTraits)){
    
    mm <- matrix(0,n,max(ordTraits[,j]) )
    jj <- as.vector( matrix(ordTraits[,j],n,s,byrow=T) )
    tmp <- byRcpp(as.vector(yWt),ii,jj,mm,mm,fun='sum')
    
    omat[,j] <- apply(tmp,1,which.max)
  }
  omat
}

#################################################
weightAve <- function(x,wt){  #x vector for variable, wt is weight

  w1 <- x*wt
  m1 <- sum(w1)/sum(wt)

  v1  <- sum(wt*x^2)/sum(wt) - m1^2

  c(m1,v1)
}

cov2cor <- function(covmat){  #covariance matrix to correlation matrix

  d    <- nrow(covmat)
  di   <- diag(covmat)
  s    <- matrix(di,d,d)
  covmat/sqrt(s*t(s))
}

cor2cov <- function(sigvec,cormat){ #correlation matrix and variance vector to covariance

  d <- length(sigvec)
  s <- matrix(sigvec,d,d)
  cormat*sqrt(s*t(s))
}

  


heatGrid <- function(x,groupIndex=NULL,legendTick=c(0,1),htScale=1){

  nr   <- nrow(x)
  nc   <- ncol(x)

  ncol <- 100
  colseq <- heat.colors(ncol)

  scale <- seq(min(x,na.rm=T),max(x,na.rm=T),length=ncol)

  x[is.na(x)] <- 0

  ww   <- as.matrix(expand.grid(c(1:nr),c(1:nc)))

  icol <- findInterval(x,scale,all.inside=T)
  coli <- colseq[icol]

  bwide <- rep(1,nrow(ww))
  bhigh <- rep(1/nr/5,nrow(ww))*htScale

  symbols(ww[,2],ww[,1]+1,rectangles=cbind(bwide,bhigh),xlim=c(0,nc+2),ylim=c(0,nr+2),
              fg=coli,bg=coli,inches=F,axes=F,xlab=' ',ylab=' ')
  
  cc <- rownames(x)
  if(!is.null(cc))text(rep(nc+1,nc),1 + c(1:nr),cc,cex=.6)
  rr <- colnames(x)
  if(!is.null(rr))text(1:nc,rep(nr+2,nc),rr,srt=45,cex=.6)

  if(!is.null(groupIndex)){
    groups <- unique(groupIndex)
    ns     <- length(groups)

    for(k in 1:ns){

       wk <- which(groupIndex == groups[k])
       points(rep(1 + .2*k,length(wk)),wk,col='white',pch=7,cex=.5)
    }
  }

  colorLegend(c(.1,.4),nr*c(.9,1),ytic=legendTick,scale=scale,cols=colseq,labside='left')
}

colorWhite2Brown <- function(ncol){
    
    colF   <- colorRampPalette(c('white','yellow','orange','red','brown'))
    colF(ncol)
}

colorWhite2Red <- function(ncol){
  
  colF   <- colorRampPalette(c('white','yellow','orange','red'))
  colF(ncol)
}

colorWheat2Red <- function(ncol){
  
  colF   <- colorRampPalette(c('wheat','yellow','orange','red'))
  colF(ncol)
}

colorWheat2Brown <- function(ncol){
  
  colF   <- colorRampPalette(c('wheat','orange','brown'))
  colF(ncol)
}

colorWhite2Black <- function(ncol){
  
  colF   <- colorRampPalette(c('white','grey','grey','grey','black','black'))
  colF(ncol)
}

colorGrey2Black <- function(ncol){
  
  colF   <- colorRampPalette(c('grey','black'))
  colF(ncol)
}


plotPredictAspect <- function(obsSA,predSA,slopePlot=T,plotType='grid',bubbleSize=.002){
  
  # slope/aspect matrices each have two columns, slope, aspect, in degrees
  #    generated by u2slopeAspect
  
  nPerBin <- nrow(obsSA)/20
  
  if(slopePlot){
    plotObsPred(obsSA[,1],predSA[,1],nPerBin=nPerBin,xlab='Degrees',ylab='Predicted',
                POINTS=F,ylim=c(0,70),fill='grey')
    plotLabel('Slope',above=T,cex=2)
    abline(0,1,lty=2)
  }
  
  if(plotType == 'grid'){
    ngrid <- 70
    aseq <- seq(-180,180,length=ngrid)
    xygrid <- as.matrix( expand.grid(aseq,aseq) )
    tmp <- t( points2grid(obsSA[,2],predSA[,2],xygrid,wt=obsSA[,1])$mat )
    
    tmp <- apply(tmp,2,rev)
    
    corPlot(tmp,tri='both',squarePlot=T,specLabs=F,colorGrad = light2darkColors,
            LEGEND=F,widex=5.5,widey=2.5)
    axis(1,at=seq(0,ngrid,length=5),labels=c('S','W','N','E','S') )
    axis(2,at=seq(0,ngrid,length=5),labels=c('S','W','N','E','S') )
    abline(0,1,lty=2)
    abline(ngrid,1,lty=2)
    abline(-ngrid,1,lty=2)
    plotLabel('Aspect',above=T,cex=1)
    
  } else {
    ww <- which(obsSA[,1] > 5)
    maxVal <- 2*max(obsSA[ww,1])
    symbols(obsSA[ww,2],predSA[ww,2],circles=obsSA[ww,1]*bubbleSize,inches=F,
            xlab='Direction',ylab=' ',xaxt='n',yaxt='n',fg='black',bg='grey')
    axis(1,at=c(-180,-180/2,0,180/2,180),labels=c('S','W','N','E','S') )
    axis(2,at=c(-180,-180/2,0,180/2,180),labels=c('S','W','N','E','S') )
    abline(0,1,lty=2)
    abline(360,1,lty=2)
    abline(-360,1,lty=2)
    plotLabel('Aspect on slopes > 5 degrees',above=T,cex=1)
  }
}

#####################################################
corPlot <- function(cmat,slim=NULL,diag0=F,plotScale=1,
                    makeColor=NULL,textSize=NULL,
                    corLines=T,tri='lower',colorGrad = NULL,
                    cex=1,specLabs=T,squarePlot=T,LEGEND=T,
                    widex=10.5,widey=6.5,add=F,new=T){  #correlation or covariance matrix

  # makeColor - list of matrices of indices for boxes
  #   names of matrices are colors
  # if(diag0)diag(cmat) <- 0
  # tri - 'lower','upper', or 'both'
  # colorGrad - mapColors, heatColors, darkColors
  # squarePlot makes symbols square
  # new means NOT NEW 

  if(tri == 'upper')cmat[lower.tri(cmat)] <- 0
  if(tri == 'lower')cmat[upper.tri(cmat)] <- 0

  dy  <- nrow(cmat)
  dx  <- ncol(cmat)
  d <- dx
  xtext <- rep(c(1,100),dx/2)
  if(length(xtext) < d)xtext <- c(xtext,1)

  if(d < 20)xtext <- xtext*0 + 1

  xtext <- xtext*0 + 1

  ncol <- 200
  colF   <- colorRampPalette(c('blue','green','white',
                               'yellow','brown'))
  colseq <- colF(ncol)

  if(is.null(slim))slim = range(cmat)
  slim  <- signif(slim,1)
  scale <- seq(slim[1],slim[2],length.out=ncol)
  
  if(slim[1] < 0 & slim[2] > 0){
    dp <- slim[2] - 0
    dm <- 0 - slim[1]
    ncol <- 200
    
    colseq <- colF(ncol)
    if(dp < dm)colseq <- colseq[101 + c(-100:round(dp/dm*100))]
    if(dp > dm)colseq <- colseq[ round((1 - dm/dp)*100):200 ]
    ncol  <- length(colseq)
  }

  ww   <- as.matrix(expand.grid(c(1:dy),c(1:dx)))  # note reverse order

  if(tri == 'upper'){
    ww  <- ww[ww[,1] <= ww[,2],]
    ww  <- ww[order(ww[,1]),]
  }
  if(tri == 'lower'){
    ww  <- ww[ww[,1] >= ww[,2],]
    ww  <- ww[order(ww[,1]),]
  }

  icol <- findInterval(cmat[ww],scale,all.inside=T)
  coli <- colseq[icol]

  if(diag0)coli[ww[,1] == ww[,2]] <- 'white'
  
  ss <- max(c(dx,dy))/5/plotScale
  
  if(squarePlot)mapSetup(c(0,dx),c(0,dy),ss,widex=widex,widey=widey)
  
  if(new)par(new = new)
  
  symbols(ww[,2],dy - ww[,1] + 1,squares=rep(1,nrow(ww)),
          xlim=c(0,dx+4),ylim=c(0,dy+4),
          fg=coli,bg=coli,inches=F,xlab=' ',ylab=' ',xaxt='n',yaxt='n',add=add)
  
  if(!is.null(makeColor)){
    
    for(k in 1:length(makeColor)){
      mm <- makeColor[[k]]
      if(length(mm) == 0)next
      if(tri == 'upper')mm <- mm[mm[,1] <= mm[,2],]
      if(tri == 'lower')mm <- mm[mm[,1] >= mm[,2],]
      ss <- matrix(0,dy,dx)
      ss[mm] <- 1
      wk <- which(ss[ww] == 1)
      ccc <- names(makeColor)[[k]]
      symbols(ww[wk,2],dy - ww[wk,1]+1,squares=rep(1,length(wk)),
              fg=ccc,bg=ccc,inches=F,xaxt='n',yaxt='n',add=T)
    }
  }

  if(corLines & tri == 'lower'){
    for(kk in 1:d){
      kb <- kk - .5
      ke <- d - kk + .5
      
      if(kk <= d)lines(c(kb,kb),c(0,ke),col='grey',lwd=1.5)          #verticle
      if(kk > 1) lines( c(.5,kb),c(ke,ke),col='grey',lwd=1.5)        #horizontal
      if(kk > 1) lines(c(kb,kb+.5),c(ke,ke+.5),col='grey',lwd=1.5)    #diagonal
    }
  }
  rect(0,-1,d+1,.5,col='white',border=NA)
  
  if(is.null(textSize))textSize <- exp(-.02*ncol(cmat))
  labels   <- rev(rownames(cmat))
  if(!specLabs)labels <- F
  
  if(tri == 'lower' & specLabs)text( c(d:1)+.1*xtext, c(1:d)+.5, rev(colnames(cmat)),pos=4,srt=45,cex=textSize)
  if(tri == 'both'){
    if(specLabs)text( c(dx:1)-.1*xtext, xtext*0+dy+.8, rev(colnames(cmat)),
                      pos=4,srt=55,cex=textSize)
    par(las = 1)
    axis(side=2,at=c(1:dy),labels=labels,tick=F,lwd=0,pos=.5,cex.axis=textSize)
    par(las = 0)
  }
  
  labside <- 'right'
  
  if(LEGEND)colorLegend(c(dx+1,dx+2),c(.5*dy,dy),ytick=c(slim[1],0,slim[2]),
                        scale,cols=colseq,labside=labside,
                        endLabels=range(slim),text.col='black')
}
##################################################

zeroOneScale <- function(x) (x - min(x,na.rm=T))/(max(x,na.rm=T) - min(x,na.rm=T))

######################################################

mapContoursUSA <- function(xx,yy,z,zlevs=NULL,yname=NULL){ #must have lon,lat in main program

  map('state',interior=F,col='grey',xlim=c(-97,-65))
  map('state',boundary=F,col='grey',add=T)

  values2contour(xx,yy,z,zlevs,lwd=1,col='black',add=T)
  if(!is.null(yname))title(yname)
 
}
######################################################

arrowField <- function(xy0,xy1=NULL,directionLong=NULL,angle=20,col=1){

  #xy0, xy1, 2 columns each
  #directionLong in radians, length

  if(is.null(xy1) & is.null(directionLong))stop('either xy1 or directionLong must be given')

  if(!is.null(directionLong)){
     xy1     <- xy0*0
     xy1[,1] <- xy0[,1] + directionLong[,2]*cos(directionLong[,1])
     xy1[,2] <- xy0[,2] + directionLong[,2]*sin(directionLong[,1])
     long    <- sqrt( (xy1[,1] - xy0[,1])^2 + (xy1[,2] - xy0[,2])^2)
  }

  lwide <- 2*zeroOneScale(long)

  arrows(xy1[,1],xy1[,2],xy0[,1],xy0[,2],length=.05*long,angle,col,lwd=lwide)
}

#################################################################3
mapArrowsUSA <- function(xx,yy,direction,long,yname=NULL){


  map('state',interior=F,col='grey',xlim=c(-97,-65))
  map('state',boundary=F,col='grey',add=T)
      

    wf <- which(is.finite(direction) & is.finite(long))

    lo <- xx[wf]
    la <- yy[wf]

    arrowField(cbind(lo,la),directionLong=cbind(direction,long)[wf,],angle=20,col=1)
 
  if(!is.null(yname))title(yname)

  map('state',interior=F,add=T,lwd=4,col='white')
  map('state',interior=F,add=T)

 # map('state',boundary=F,col='white',lwd=3,add=T)
  map('state',boundary=F,col='grey',add=T)
 # title(mapname)

}



regMap <- function(x,y,z=NULL,fg=rep(1,length(x)),bg=rep(1,length(x)),
                   xl=range(x),yl=range(y),lineCol='black',axes=T,lwd=1,
                   zminmax=NULL,scaleSym=1,POINTS=F,
                   IMAGE=F,ADD=F,mapscale=1,county=F,
                   stateVec=c('north carolina','south carolina','tennessee',
                              'georgia','virginia'),region='dummy'){
  
  # if IMAGE z is a matrix with length(x) rows, length(y) columns)
  # if !IMAGE, xx, yy, and z are all same length
  require(maps)
 # require(maptools)
  require(mapdata)
  # require(sp)
  # require(rworldmap)
  
  # if(region == 'andes'){
  #   data(countryExData)
  #   sPDF <- joinCountryData2Map( countryExData, 
  #                                joinCode = "ISO3", 
  #                                nameJoinColumn = "ISO3V10")
  #   mapCountryData( sPDF, nameColumnToPlot="BIODIVERSITY" )
  
  #   con <- url('http://gadm.org/data/rda/CHE_adm1.RData')
  #   print( load(con))
  #   close(con)
  #   spplot(gadm)
  # }
  
  if(!ADD){
    mapSetup(xl,yl,scale=mapscale)
  }
  
  colSea   <- colorRampPalette(c('black','darkblue','blue','lightblue','white'))
  
  cols <- c(colSea(5),terrain.colors(100))
  
  zr <- max(z)
  bb <- c( -10000,seq(-300,-5,length=5), seq(-2,160,length=50), seq(170,zr+1,length=50))
  
  if(!is.null(zminmax)){
    bb <- seq(zminmax[1],zminmax[2],length=106)
  }
  
  if(IMAGE){
    image(x,y,z,col=cols,breaks=bb,add=ADD,xlab=' ',ylab=' ',xlim=xl,ylim=yl) 
    ADD <- T
  }
  
  #  map('county','north carolina',boundary=F,col='grey',add=T)
  
  if(!is.null(stateVec) & county & region != 'andes'){
    map('county',stateVec,boundary=T,col='grey',lwd=lwd,
        add=ADD,xlim=range(x),ylim=range(y))
    ADD <- T
  }
  if(region != 'andes')map('state',add=ADD,col=lineCol,lwd=lwd,xlim=xl,ylim=yl)
  
  #  map('worldHires','Peru',add=ADD,col=lineCol,lwd=lwd,xlim=xl,ylim=yl)
  #  map('county',stateVec,interior=F,col='black',add=T)
  if(axes)map.axes()
  
  if(IMAGE)return()
  
  if(POINTS){
    if(is.null(z)) points(x,y,col=bg)
    if(!IMAGE & !is.null(z) & length(x) == length(y)){
      symbols(x,y,circles=z*scaleSym,fg=fg,bg=bg,inches=F,add=T)
    }
  }
}


#######################################3


mapPointsUSA <- function(xx,yy,...,sym='values',cols='black',whiteOutDist=NULL,
                         fill='none',yname=NULL,scale=NULL,legendScale=NULL,KRIG=F){  

  #scale - symbol size
  #legendScale - range of values for color scale

  map('state',interior=F,col='grey',xlim=c(-97,-65))
  map('state',boundary=F,col='grey',add=T)

  zz  <- list(...)

  if(length(zz) > 1){
    if(length(sym) == 1) sym  <- rep(sym,length(x))
    if(length(cols) == 1)cols <- rep(cols,length(x))
    if(length(fill) == 1)fill <- rep(fill,length(x))
  }

  if('ramp' %in% cols | 'ramp' %in% fill){

    ncol   <- 100
    colseq <- mapColors(ncol)

  }

  for(j in 1:length(zz)){

    xj <- zj <- zz[[j]]
    wf <- which(is.finite(xj))
    xj <- xj[wf]

    if(is.null(legendScale)) colScale <- seq(min(xj,na.rm=T),max(xj,na.rm=T),length=ncol)
    if(!is.null(legendScale))colScale <- seq(legendScale[1],legendScale[2],length=ncol)

    lo <- xx[wf]
    la <- yy[wf]

    if(sym[j] == 'ones'){
       w0 <- which(xj == 1)
       xj <- rep(.1,length(w0))
       ll <- la[w0]
       lo <- lo[w0]
    }

    fj <- cols[j]
    bj <- fill[j]
    if(bj == 'none')bj <- NA

    if(cols[j] == 'ramp' | fill[j] == 'ramp'){
      fj <- bj <- colseq[findInterval(xj,colScale)]
    }

    if(KRIG){
        tmp <- surf.gls(2,expcov,x=lo,y=la,xj,d=.7)
        ps  <- prmat(tmp,xl=min(lo),xu=max(lo),yl=min(la),yu=max(la),n=50)
        xy  <- as.matrix(expand.grid(ps$x,ps$y))
        zk  <- as.vector(ps$z)
        zk[zk < min(colScale)] <- colScale[1] + .1
        zk[zk > max(colScale)] <- max(colScale) - .1
        fj <- bj <- colseq[findInterval(zk,colScale)]
        symbols(xy[,1],xy[,2],squares=rep(.5,nrow(xy)),inches=F,fg=fj,bg=bj,add=T)
        zj <- ps$z

        if(!is.null(whiteOutDist)){   #whiteout areas further than this distance from points
            xy1  <- xy[,1]
            xy2  <- xy[,2]
            tmp  <- distmat(xx,yy,xy1,xy2)
            mind <- apply(tmp,1,min,na.rm=T)
            ww   <- which(mind > whiteOutDist)
            symbols(xy1[ww],xy2[ww],squares=rep(.5,length(ww)),inches=F,fg='white',bg='white',add=T)
        }
    }

    if(!KRIG){
      if(is.null(scale)) xj <- zz[[j]]
      if(!is.null(scale))xj <- rep(scale,length(zz[[j]]))
      symbols(xx,yy,circles=xj,inches=F,fg=fj,bg=bj,add=T)
      zj <- xj
    }

    if(is.null(legendScale))legendScale <- signif(range(zj,na.rm=T),1)
    
    tscale <- legendScale
    ttic   <- tscale
    if(tscale[1] < 0 & tscale[2] > 0)ttic   <- c(tscale[1],0,tscale[2])
    tscale <- seq(tscale[1],tscale[2],length=ncol)

    polygon(c(-75,-65,-65,-75),c(25,25,36,36),col='white',border=NA)
    colorLegend(c(-70,-69),c(26,35),ytick=ttic,scale=tscale,colseq,labside='left')
  }

  if(!is.null(yname))title(yname)

  map('state',interior=F,add=T,lwd=4,col='white')
  map('state',interior=F,add=T)

 # map('state',boundary=F,col='white',lwd=3,add=T)
  map('state',boundary=F,col='grey',add=T)
 # title(mapname)

}



binMultDim <- function(xyData,xyGrid=NULL,nbin=rep(10,ncol(xyData)),
                       xyRange=NULL,LOG=F,allIn=T){
  
  #xyData  - matrix, columns are locations in space
  #xyGrid  - list of grid points in each dimension
  #          if xyGrid not provided it will be created from nbin
  #xyRange - 2 by d matrix, 1st row low values, 2nd row high values
  
  tiny <- min(xyData[xyData > 0],na.rm=T)
  
  d    <- ncol(xyData)
  nn   <- nrow(xyData)
  
  if(is.null(xyRange))xyRange <- apply(xyData,2,range,na.rm=T)
  
  if(LOG & min(xyRange) <= 0){
    warning("negative values omitted on log scale")
    xyData[xyData < tiny] <- tiny
    xyRange[xyRange < tiny] <- tiny
  }
    
  if(is.null(xyGrid)){  # must have nbin vector
    
    xyGrid <- numeric(0)
    for(j in 1:d){
      if(!LOG)jbins <- seq(xyRange[1,j],xyRange[2,j],length=nbin[j])
      if(LOG) jbins <- 10^seq(log10(xyRange[1,j]),log10(xyRange[2,j]),length=nbin[j])
      xyGrid <- append(xyGrid,list( jbins ) )
    }
    
  } else {
    if(ncol(xyData) != length(xyGrid))stop('xyGrid must be a list with one vector per column in xyData')
  }
  
  bins <- numeric(0)
  
  for(j in 1:d){
    jbin <- findInterval(xyData[,j],xyGrid[[j]])
    bins <- cbind(bins,jbin)
    
  }
  colnames(bins) <- paste('d',c(1:d),sep='-')
  
  if(allIn){
    for(j in 1:d){
      bins[bins[,j] < 1,j] <- 1
      bins[bins[,j] > length(xyGrid[[j]]),j] <- 1
    }
  }
  
  nb <- apply(bins,2,max)
  
  index <- as.data.frame(bins)
  
  tmp <- aggregate(rep(1,nn),by=index,FUN=sum)
  
  if(!allIn){
    tt <- matrix( apply(tmp,2,max), nrow(tmp),ncol(tmp),byrow=T )
    ww <- which( tmp > tt | tmp < 1, arr.ind=T)
    if(length(ww) > 0){
      wr <- unique(ww[,1])
      tmp <- tmp[-wr,]
    }
  }
  
  binWidth <- binWideMat <- numeric(0)
  
  for(j in 1:d){
    dx <- diff(xyGrid[[j]])
    dx <- c(dx,dx[length(dx)])
    
    binWidth   <- append(binWidth,list(dx))
    binWideMat <- cbind(binWideMat,dx[tmp[,j]])
    tmp[,j]    <- xyGrid[[j]][tmp[,j]]
  }
  
  list(data2bins = bins, gridSums = tmp, xyGrid = xyGrid, 
       binWidth = binWidth, binWideMat = binWideMat)
}



mapMask <- function(xx,yy,dx=NULL,dy=NULL,whiteOutDist=1,col='white'){    
  
  #mask parts of map with no data, grid density (dx,dy)
  #add to exising map
  #xy is an expand.grid of unmasked pixels
  
  require(RANN)
  
  rx <- range(xx,na.rm=T)
  ry <- range(yy,na.rm=T)
  
  if(is.null(dx)){
    dx <- diff(rx)/50
    dy <- diff(rx)/50
  }
  
  
  xnew <- seq(rx[1],rx[2],by=dx)
  ynew <- seq(ry[1],ry[2],length=length(xnew))
  xxy   <- as.matrix( expand.grid(xnew,ynew) )
  nr   <- nrow(xxy)
  
  nx <- 0
  wd <- whiteOutDist*1.5
  
  while(nx < min(c(5000,nr/4))){
    tmp <- nn2(cbind(xx,yy),xxy,k = 1)$nn.dists
    xy  <- xxy[tmp > wd,]
    nx  <- nrow(xy)
    if(is.null(nx))return()
    wd  <- wd/1.3
  }
  
  symbols(xy[,1],xy[,2],squares=rep(dx,nx),inches=F,fg=col,bg=col,add=T)
  
}


checkDesign <- function(x,intName='intercept',xflag=':'){  # name of intercept column

  # xflag - indicates that variable is an interaction
  
  p <- ncol(x)
  
  if(ncol(x) < 3){
    message('no design check, x has 2 columns')
    return( list(VIF = 0, correlation = 1, rank = 2, p = 2, isFactor=character(0)) )
  }
    
  if(is.null(colnames(x))){
    colnames(x) <- paste('x',c(1:p),sep='_')
  }
  xrange      <- apply(x,2,range)
  wi          <- which(xrange[1,] == 1 & xrange[2,] == 1)
  if(length(wi) > 0)colnames(x)[wi] <- 'intercept'
  
  wx <- grep(xflag,colnames(x))
  wi <- which(colnames(x) == 'intercept')
  wi <- unique(c(wi,wx))

  xname <- colnames(x)
  
  wmiss <- which(is.na(x),arr.ind=T)
  
  if(length(wmiss) > 0){
    rowTab <- table( table(wmiss[,1]) )
    colTab <- table(wmiss[,2])
  }
    
  VIF <- rep(NA,p)
  names(VIF) <- xname
  
  isFactor <- character(0)
  
  for(k in 1:p){

    if(xname[k] %in% xname[wi])next
    
    notk <- xname[xname != xname[k] & !xname %in% xname[wi]]

    ykk <- x[,xname[k]]
    xkk <- x[,notk]

    tkk <- summary(lm(ykk ~ xkk))$adj.r.squared
    VIF[k] <- 1/(1 - tkk)
    
    xu <- sort( unique(x[,k]) )
    tmp <- identical(c(0,1),xu)
    if(tmp)isFactor <- c(isFactor,xname[k])
  }

  VIF <- VIF[-wi] 

  corx <- cor(x[,-wi])
  rankx <- qr(x)$rank
  corx[upper.tri(corx,diag=T)] <- NA
  
  findex <- rep(0,p)
  
  findex[xname %in% isFactor] <- 1
  
  designTable <- list('table' = rbind( round(VIF,1),findex[-wi],round(corx,2)) )
  rownames(designTable$table) <- c('VIF','factor',xname[-wi])
  
  designTable$table <- designTable$table[-3,]
  
  if(p == rankx)designTable$rank <- paste('full rank:',rankx,'= ncol(x)')
  if(p < rankx) designTable$rank <- paste('not full rank:',rankx,'< ncol(x)')

  list(VIF = round(VIF,1), correlation = round(corx,2), rank = rankx, p = p,
       isFactor = isFactor, designTable = designTable)
}


############################################3
points2angle <- function(xy0,xy){            
  
  #(x,y) columns in xy from reference (x,y) in xy0, returns angle

  h  <- distmat(xy0[1],xy0[2],xy[,1],xy[,2])
  dx <- xy[,1] - xy0[1]
  dy <- xy[,2] - xy0[2]

  list(theta = atan2(dy,dx), dist = h)
}
#################################################
localSlope <- function( xk, yk, zk, direction=NULL, na.rm=F ){ 
  
  # if not null, direction is theta
  # theta in radians, 0 to the right, pi to the left
  # positive theta is up, negative theta is down
  # theta points uphill, aspect points downhill
  # aspect is 0, 360 up (north), 90 right (east), 180 down (south)
  
  if(na.rm){
    ww <- which(is.finite(xk) & is.finite(yk) & is.finite(zk))
    xk <- xk[ww]
    yk <- yk[ww]
    zk <- zk[ww]
  }

  x0 <- mean(xk)
  y0 <- mean(yk)
  z0 <- mean(zk)

  xk <- xk - x0
  yk <- yk - y0
  zk <- zk - z0

  u1 <- u2 <- 1

  X <- cbind(xk,yk)

  if(length(xk) == 2){
    fx <- diff(zk)/diff(xk)
    fy <- diff(zk)/diff(yk)
  }

  if(length(xk) > 3){
    b <- invMat(crossprod(X),NEARPD=T)%*%crossprod(X,zk)
    fx <- b[1,]
    fy <- b[2,]
  }

  if(is.null(direction)){     #maximum derivative
    grade <- sqrt(fx^2 + fy^2)
    theta <- atan2(fy,fx)
  }

  if(!is.null(direction)){    #derivative in direction theta
     theta <- direction
     u1    <- cos(direction)
     u2    <- sin(direction)
     grade <- fx*u1 + fy*u2
  }
  
  deg <-  radians2degrees(pi/2 - theta) + 180
  if(deg > 360)deg <- deg - 360

  list(theta = theta, grade = grade, aspect = deg)
}

radians2degrees <- function(rads){
  
  rads/2/pi*360
}


mapColors <- function(ncol){

    colF   <- colorRampPalette(c('darkblue','blue','lightblue',
                                 'green','lightgreen',
                                 'yellow','orange','red','brown'))
    colF(ncol)
}

heatColors <- function(ncol){
  
  colF   <- colorRampPalette(c('darkblue','blue','brown','red'))
  colF(ncol)
}

darkColors <- function(ncol){
  
  colF   <- colorRampPalette(c('blue','green','red','brown','black'))
  colF(ncol)
}


light2darkColors <- function(ncol){
  
  colF   <- colorRampPalette(c('white','orange','brown','black'))
  colF(ncol)
}

degrees2radians <- function(degrees){
  
  degrees/360 * 2 * pi
}


clim4FIA <- function(lat,lon,plVec,yrVec,nplot=max(plVec),
                     monTemp=NULL,monPrec=NULL,monThermIndex=c(1:12),
                     start=1969,end=2012,GINDEX=F,ONEYRVEC=F,LST=T, 
                     path="climate/"){
  
  # if ONEYRVEC there is one yrVec that applies to all plots
  # if !ONEYRVEC there is a vector containing a year for each lat,lon
  
  tmp <- inClimate(lat=lat,lon=lon,start=start,end=end,LST=LST,
                   PADLAST=T, path=path)
  temp <- tmp$temp
  prec <- tmp$prec
  
  endy  <- max( round( as.numeric(colnames(temp),0 ) ) )
  npp <- max(plVec)
  nyy <- endy - start + 1
  yiIndex <- start:endy
  
  gindex <- numeric(0)
  
  monthNames <- getMonthNames()
  
  regClim <- matrix(NA,length(plVec),2)
  colnames(regClim) <- c('temp','prec')
  
  if(GINDEX){
    yvec   <- rep(yiIndex,each=12)
    mvec   <- rep(1:12,nyy)
    ccols  <- as.character(round( yvec + mvec*.01,2 ))
    
    wshort <- which(nchar(ccols) == 6)
    if(length(wshort) > 0)ccols[wshort] <- paste(ccols[wshort],'0',sep='')
    ccols  <- ccols[ccols %in% colnames(temp)]
    tmp    <- monthlyPHr(yi=yvec,mi=mvec,
                         tempMatrix=temp[,ccols],precMatrix=prec[,ccols],lat=lat)
    gindex <- tmp$degHrPosYr
    gindex <- gindex[,yvec - start + 1]
    dindex <- tmp$degHrNegYr
    dindex <- dindex[,yvec - start + 1]
    regClim <- matrix(NA,length(plVec),4)
    colnames(regClim) <- c('temp','prec','therm','deficit')
  }
  
  
  lastYear <- max( floor( as.numeric(colnames(temp)) ) )
  
  monTempIndex <- monPrecIndex <- c(1:12)
  
  if(!is.null(monTemp)){
    moVec <- numeric(0)
    s1 <- seq(1,nchar(monTemp),by=3)
    s2 <- s1 + 2
    monTempIndex <- match( substring(monTemp,s1 , s2 ), monthNames)
  }
  if(!is.null(monPrec)){
    moVec <- numeric(0)
    s1 <- seq(1,nchar(monPrec),by=3)
    s2 <- s1 + 2
    monPrecIndex <- match( substring(monPrec,s1 , s2 ), monthNames)
  }
  
  cCols <- matrix( unlist( strsplit(colnames(temp),'[.]') ), ncol=2,byrow=T)
  yrCol <- as.numeric(cCols[,1])
  moCol <- as.numeric(cCols[,2])
  
  for(j in 1:npp){
    
    wj <- which(plVec == j)   #yrs for plot j
    nw <- length(wj)
    if(nw == 0)next
    
    if(ONEYRVEC) yjj <- yrVec
    if(!ONEYRVEC)yjj <- yrVec[wj]
    
    if(nw == 1){       #only one row for plot
      
      yjj[yjj > endy] <- endy
      
      wt <- which( yrCol %in% yjj & moCol %in% monTempIndex )
      wp <- which( yrCol %in% yjj & moCol %in% monPrecIndex )
      wg <- which( yrCol %in% yjj & moCol %in% c(1:12) )
      
      tj <- temp[j,wt]
      pj <- prec[j,wp]
      
      meanTemp <- byIndex(tj,yrCol[wt],mean,na.rm=T)
      meanPrec <- byIndex(pj,yrCol[wp],sum,na.rm=T)
      
      regClim[wj,'temp']    <- mean( meanTemp, na.rm=T)
      #    regClim[wj[k],'temp'] <- meanTemp[length(meanTemp)]
      
      regClim[wj,'prec'] <- mean( meanPrec, na.rm=T)
      #    regClim[wj[k],'prec'] <- meanPrec[length(meanPrec)]
      
      if(GINDEX){
        #      wg <- which( yiIndex %in% yjj  )
        #      gj <- gindex[j,wg]
        
        gj <- gindex[j,wt]
        gj <- byIndex(gj,yrCol[wp],mean,na.rm=T)
        regClim[wj,'therm']    <- mean(gj)
        
        dj <- dindex[j,wt]
        dj <- byIndex(dj,yrCol[wp],mean,na.rm=T)
        regClim[wj,'deficit']    <- mean(dj)
        
        #      regClim[wj[k],'therm'] <- gj[length(gj)]
      }
      
      next
    }
    
    if(nw > 1){
      
      for(k in 2:nw){
        
        kk <- wj[k-1]:(wj[k] - 1)
        yj <- yrVec[kk]
        yj[yj >= lastYear] <- lastYear
        if(length(yj) > 1){
          if((yj[2] - yj[1]) > 1)yj <- yj[1]:yj[2]
        }
        
        wt <- which( yrCol %in% yj & moCol %in% monTempIndex )
        wp <- which( yrCol %in% yj & moCol %in% monPrecIndex )
        wg <- which( yrCol %in% yj & moCol %in% c(1:12) )
        
        tj <- temp[j,wt]
        pj <- prec[j,wp]
        
        meanTemp <- byIndex(tj,yrCol[wt],mean,na.rm=T)
        meanPrec <- byIndex(pj,yrCol[wp],sum,na.rm=T)
        
        regClim[kk,'temp']    <- mean( meanTemp, na.rm=T)
        regClim[wj[k],'temp'] <- meanTemp[length(meanTemp)]
        
        regClim[kk,'prec'] <- mean( meanPrec, na.rm=T)
        regClim[wj[k],'prec'] <- meanPrec[length(meanPrec)]
        
        if(GINDEX){
          #     yj[yj > end] <- end
          #     wg <- which( yiIndex %in% yj  )
          #     gj <- gindex[wj[k],wg]
          
          gj <- gindex[j,wg]
          gj <- byIndex(gj,yrCol[wg],mean,na.rm=T)
          regClim[kk,'therm']    <- gj
          regClim[wj[k],'therm'] <- gj[length(gj)]
          
          dj <- dindex[j,wg]
          dj <- byIndex(dj,yrCol[wg],mean,na.rm=T)
          regClim[kk,'deficit']    <- dj
          regClim[wj[k],'deficit'] <- dj[length(dj)]
        }
      }
    } #end if nw > 1
    
  }
  
  if(!is.null(monTemp)){
    colnames(regClim)[colnames(regClim) == 'temp'] <- paste('temp',monTemp,sep='-')
  }
  if(!is.null(monPrec)){
    colnames(regClim)[colnames(regClim) == 'prec'] <- paste('prec',monPrec,sep='-')
  }
  
  list(regClim = regClim, gindex = gindex)
}

inClimate <- function(lat=plotLat, lon=plotLon,start=1969,end=2012,LST=T,PADLAST=F,
                      path="climate/"){
  
  #LST     - use modis LST
  #PADLAST - pad missing values at end with most recent available
  
  require(RANN)
  
  se <- paste(start,end,sep='_')
  pname <- paste('ppt_800m',se,'monthly.RData',sep='_')
  pname <- paste(path,pname,sep='')
  print(pname)
  load(pname)
  # ppt.dat <- out
  precDat    <- ppt.dat[,-c(1:2)]
  precPlots  <- as.numeric( rownames(ppt.dat) )
  precLonLat <- ppt.dat[,c('lon','lat')]
  
  if(!LST){
    rname <- paste('tmax_800m',se,'monthly.RData',sep='_')
    sname <- paste('tmin_800m',se,'monthly.RData',sep='_')
    rname <- paste(path,rname,sep='')
    sname <- paste(path,sname,sep='')
    
    print(rname)
    load(rname)
    #  tmax.dat <- out
    load(sname)
    #  tmin.dat <- out
    
    temp <- tmax.dat
    temp[,-c(1:2)] <- ( tmax.dat[,-c(1:2)] + tmin.dat[,-c(1:2)] )/2
    tempLonLat <- temp[,c('lon','lat')]
    tmp        <- nn2(tempLonLat[,c('lat','lon')], cbind(lat,lon),k = 1  )$nn.idx
    temp <- tempDat <- as.matrix(temp[tmp,-c(1:2)])
  }
  
  if(LST){
    
    load( paste( path,'LST_estimates/calibLST.Rdata',sep='') )
    
    tempPlots  <- rownames( calibLST )
    tempLonLat <- calibLST[,c('lon','lat')]
    tempDat    <- calibLST[,-c(1:2)]
    
    tmp  <- nn2(tempLonLat[,c('lat','lon')], cbind(lat,lon),k = 1  )$nn.idx
    temp <- as.matrix(tempDat[tmp,])
  }
  
  wrr <- 1
  
  while(length(wrr) > 0){
    tmp  <- nn2(ppt.dat[,c('lat','lon')], cbind(lat,lon),k = 1  )$nn.idx
    prec <- as.matrix(precDat[tmp,])
    
    ps <- rowSums(prec*0 + 1)
    wrr <- which(is.na(ps))
    
    if(length(wrr) > 0)precDat <- precDat[-tmp[wrr],]
  }
  
  mcol <- match(colnames(temp),colnames(prec))  #temp has more data (columns) than prec
  ww   <- which(is.na(mcol))
  if(length(ww) > 0){
    tmp <- matrix(NA,nrow(prec),length(ww))
    colnames(tmp) <- colnames(temp)[ww]
    prec <- cbind(prec,tmp)
  }
  
  wshort <- which(nchar(colnames(temp)) == 6)
  if(length(wshort) > 0)colnames(temp)[wshort] <- pasteChar2End(colnames(temp)[wshort],'0')
  wshort <- which(nchar(colnames(prec)) == 6)
  if(length(wshort) > 0)colnames(prec)[wshort] <- pasteChar2End(colnames(prec)[wshort],'0')
  
  if(PADLAST){
    prec <- padClimateData(prec,startYr = 2011)
    temp <- padClimateData(temp,startYr = 2011)
  }
  
  list(temp = temp, prec = prec) 
  
}

padClimateData <- function(data,startYr=2010,endYr){
  
  # pad missing values at end with most recent available
  # colnames - "1969.01" "1969.02" "1969.03" "1969.04" "1969.05"...
  # one row per location
  # startYr - earliest yr to pad
  
  numericNames <- as.numeric(colnames(data))
  yrCol <- round( numericNames,0 )
  moCol <- round( 100*( numericNames - yrCol ), 0 )
  ncol  <- length(yrCol)
  
  if(missing(endYr))endYr <- max( yrCol )
  
  maxYrData <- yrCol[ncol]
  cols4Last <- numeric(0)         
  
  if(moCol[ncol] != 12){               #remaining months of last year
    mm        <- 12 - moCol[ncol]
    moreCols <- matrix(NA,nrow(data),mm) 
    mnames   <- substr(as.character(c(1:12)/100),2,4)[(moCol[ncol]+1):12]
    colnames(moreCols) <- paste(maxYrData,mnames,sep='')
    data <- cbind(data,moreCols)
  }
  if(maxYrData < endYr){              #remaining years
    moreCols <- matrix(NA,nrow(data),12*(endYr - maxYrData))
    mnames   <- substr(as.character(c(1:12)/100),2,4)
    cnames   <- as.vector( t( outer( (maxYrData+1):endYr,mnames,paste,sep='') ) )
    colnames(moreCols) <- cnames
    data <- cbind(data,moreCols)
  }
  wshort <- which(nchar(colnames(data)) == 6)
  if(length(wshort) > 0)colnames(data)[wshort] <- pasteChar2End(colnames(data)[wshort],'0')
  
  wyr   <- max( grep(startYr + 1,colnames(data)) ) 
  wmiss <- which(!is.finite(data),arr.ind=T)
  wmiss <- wmiss[wmiss[,2] > wyr,]
  
  while(length(wmiss) > 0){
    wcols <- sort(unique(wmiss[,2]))
    yrMiss <- colnames(data)[wcols]
    yrB4   <- as.character( as.numeric(yrMiss) - 1 )
    wshort <- which(nchar(yrB4) == 6)
    if(length(wshort) > 0)yrB4[wshort] <- pasteChar2End(yrB4[wshort],'0')
    
    yrB4   <- match(yrB4,colnames(data))
    wfrom  <- yrB4[ match( colnames(data)[wmiss[,2]],yrMiss ) ]
    data[wmiss] <- data[ cbind(wmiss[,1],wfrom) ]
    wmiss <- which(!is.finite(data),arr.ind=T)
    if(length(wmiss) == 2)wmiss <- matrix(wmiss,1,2)
  }
  data
}



#############################################

climGradient <- function( xx, yy, z1, z2, dthreshold=1, mapname=NULL,
                          climGrad=F, climDir=F, vegDir=F, ABS=T, DYDX=T, DTHETA=F){  
  
  # how closely does the gradient in z2 align with the gradient in z1?
  # directional gradient in z2 relative to maximum gradient in z1
  # climGrad - map of climate gradient derivative
  # climDir  - map of climate direction
  # vegDir   - vegetation, rather than climate direction
  # abs value- when there are multiple species
  
  ncol <- 100
  cols <- mapColors(ncol)
  
  if(!is.matrix(z2))z2 <- as.matrix(z2)
  r   <- ncol(z2)
  
  out <- matrix(NA,length(xx),5)
  colnames(out) <- c('theta','grad','thetaV','dtheta','dydx')
  
  distance <- distmat(xx,yy,xx,yy)
  
  for(i in 1:length(xx)){
    
    wi <- unique(which(distance[i,] < dthreshold))
    if(length(wi) < 4)wi <- order(distance[i,])[1:4]
    tmp   <- localSlope(xx[wi],yy[wi],z1[wi])
    theta <- tmp$theta
    dx    <- tmp$grade
    
    wz <- which(z2[wi,] == 0)  #zeros do not occur at this location
    cdir <- theta
    
    gvec <- vslp <- rep(NA,r)
    
    for(j in 1:r){
      zj  <- z2[wi,j]
      #   wj  <- which(zj != 0)
      #   if(length(wj) < 4)next
      if(vegDir)cdir <- NULL
      #   tmp <- localSlope(xx[wi[wj]],yy[wi[wj]],zj[wj],direction=cdir)
      
      vslp[j] <- localSlope(xx[wi],yy[wi],zj)$theta
      tmp  <- localSlope(xx[wi],yy[wi],zj,direction=cdir)
      gvec[j] <- tmp$grade
    }
    
    dydx <- gvec/dx
    if(ABS)dydx <- abs(dydx)
    vslp    <- mean(vslp,na.rm=T)
    dtheta  <- atan2(sin(theta - vslp),cos(theta - vslp))
    out[i,] <- c(theta,dx,vslp,dtheta,mean(dydx,na.rm=T))
  }
  
  qGrad <- c(.01,.99)
  qdydx <- c(.01,.99)
  
  # tmap <- cos( out[,'theta'] )
  tmap <- out[,'theta']
  
  gmap <- zeroOneScale( trim( out[,'grad'],qGrad) )
  
  dtmp <- trim( out[,'dydx'],qdydx)
  # dtmp <- log10(dtmp)
  dtmp[dtmp == -Inf] <- min(dtmp[dtmp != -Inf],na.rm=T)
  
  dmap <- zeroOneScale( dtmp )
  dmap[out[,'dydx'] == 0] <- min(dmap,na.rm=T)
  
  mmap <- zeroOneScale( abs(out[,'dtheta']) )
  
  
  # dmap <- zeroOneScale( trim( out[,'dydx'],qdydx) )
  
  mscale <- signif( c(0,pi),2)
  tscale <- signif(c(-pi,pi),2)
  gscale <- signif(quantile(out[,'grad'],qGrad,na.rm=T),1)
  dscale <- signif(quantile(out[,'dydx'],qdydx,na.rm=T),1)
  
  ttic   <- c(tscale[1],0,tscale[2])
  gtic   <- gscale
  dtic   <- dscale
  mtic   <- mscale
  
  if(gtic[1] < 0 & gtic[2] > 0)gtic <- c(gtic[1],0,gtic[2])
  if(dtic[1] < 0 & dtic[2] > 0)dtic <- c(dtic[1],0,dtic[2])
  
  tscale <- seq(tscale[1],tscale[2],length=ncol)
  gscale <- seq(gscale[1],gscale[2],length=ncol)
  dscale <- seq(dscale[1],dscale[2],length=ncol)
  mscale <- seq(mscale[1],mscale[2],length=ncol)
  
  if(climDir){
    mapPointsUSA(xx,yy,tmap,sym='values',cols='ramp',whiteOutDist=1.4,
                 fill='none',yname=paste(mapname,'theta',sep=' '),KRIG=T)
    polygon(c(-75,-65,-65,-75),c(25,25,36,36),col='white',border=NA)
    
    colorLegend(c(-70,-69),c(26,35),ytick=ttic,scale=tscale,cols,labside='left')
  }
  if(climGrad){
    mapPointsUSA(xx,yy,gmap,sym='values',cols='ramp',whiteOutDist=1.4,
                 fill='none',yname=paste(mapname,'Gradient',sep=' '),KRIG=T)
    polygon(c(-75,-65,-65,-75),c(25,25,36,36),col='white',border=NA)
    colorLegend(c(-70,-69),c(26,35),ytick=gtic,scale=gscale,cols,labside='left')
  }
  
  if(DYDX){
    mapnew <- paste(mapname,'dxdy',sep=' ')
    
    mapPointsUSA(xx,yy,dmap,sym='values',cols='ramp',whiteOutDist=1.4,
                 fill='none',yname=mapnew,KRIG=T)
    polygon(c(-75,-65,-65,-75),c(25,25,36,36),col='white',border=NA)
    colorLegend(c(-70,-69),c(26,35),ytick=dtic,scale=dscale,cols,labside='left')
  }
  if(DTHETA){
    mapnew <- paste(mapname,'del theta',sep=' ')
    
    mapPointsUSA(xx,yy,-(mmap-1),sym='values',cols='ramp',whiteOutDist=1.4,
                 fill='none',yname=mapnew,KRIG=T)
    polygon(c(-75,-65,-65,-75),c(25,25,36,36),col='white',border=NA)
    colorLegend(c(-70,-69),c(26,35),ytick=rev(-mtic),scale=mscale,cols,labside='left')
  }
  
  invisible(out)
}

################################################3333

associationPlot <- function(z,groupIndex=c(1:nrow(z))){       

  #z - n rows by c columns give abundance for each species
  #groupIndex - length n

  allnames <- sort(unique(groupIndex))
  ns       <- length(allnames)
  nf       <- ncol(z)

  ncol <- 50
  scale <- seq(-1,1,length.out=ncol)
  
  freqY <- numeric(0)

  nr <- round(ns/2,0) + 1
  par(mfrow=c(nr,2),bty='n',mar=c(2,2,3,1))

  for(i in 1:ns){

      yi <- z[groupIndex == allnames[i],]

      yc <- cor(yi)

      corPlot(yc,slim=c(-1,1),diag0=T)
      
      title(allnames[i])
  }
  dev.print(device=postscript,file='associationPlot.ps',width=6,horizontal=F)
}

  
############################################################

colorLegend <- function(xx,yy,ytick=NULL,scale=seq(yy[1],yy[2],length=length(cols)),
                        cols,labside='right', text.col=NULL, bg=NULL, 
                        endLabels=NULL, endCols=NULL, cex = 1){  
  # xx = (x1,x2), y = (y1,y2)
  # bg is color of border
  
  nn <- length(scale)
  ys <- seq(yy[1],yy[2],length=nn)
  
  for(j in 1:(length(scale)-1)){
    rect(xx[1],ys[j],xx[2],ys[j+1],col=cols[j],border=NA)
    if(!is.null(ytick))text(xx[1],ys[j],ytick[j], pos=2, cex = cex)
  }
  if(!is.null(ytick))text(xx[1],ys[j+1],ytick[j+1], pos=2, cex = cex)
  
  if(!is.null(bg))rect(xx[1],yy[1],xx[2],yy[2],border=bg,lwd=3, cex = cex)
  
  if(!is.null(endLabels)){ 
    cx <- 'black'
    if(!is.null(endCols))cx <- cols[c(1,nn)]
    if(!is.null(text.col))cx <- text.col
    if(!is.null(text.col))cx <- text.col
    if(labside == 'right')text(diff(xx)+c(xx[2],xx[2]),yy,endLabels,col=cx, cex = cex)
    if(labside == 'left')text(c(xx[1],xx[1]),yy,endLabels,pos=2,col=cx, cex = cex)
  }
}

####################################################
bmultiProp <- function(r,k,b = matrix(0,k,r-1),pBVar=diag(.1,k*(r-1)),
                       loB=NULL,hiB=NULL){  
	
    b <- as.vector(b)
    cvec <- myrmvnorm(1,t(b),pBVar)
    cc   <- matrix(cvec,k,r-1)
    
    if(!is.null(loB)){                     
    	cvec <- tnorm.mvt(b,b,pBVar,loB,hiB,times=1)
      cc   <- matrix(cvec,k,r-1)
    }

  list(cc = cc, cvec = cvec)
}

#################################################################
simX <- function(n,loX,hiX){                #generate design matrix

  k <- length(loX)
  x <- matrix(1,n,k)
  for(j in 1:k)x[,j] <- runif(n,loX[j],hiX[j])
  x
}
####################################################
simY <- function(x,b,LIKE,r = 1,size=rep(1,nrow(x)),sigma = 0,Effort = 1){     #simulate response

  u <- x%*%b

  if(LIKE == 'norm')return( rnorm(length(u),u,sqrt(sigma)) )
  if(LIKE == 'pois'){
  	 u <- exp(u)*Effort
  	 return( rpois(nrow(x),u) )
  } 
  if(LIKE == 'binom'){
  	 u <- invlogit(u)
  	 return( rbinom(nrow(x),1,u) )
  }
  if(LIKE == 'multinom'){
     zs <- rowSums(exp(u))
     z1 <- 1/(1 + zs)
     zm <- exp(u)/ (1 + zs)
     u  <- cbind(zm,z1)
     return( myrmultinom(size,u) )
  }
  if(LIKE == 'mvnorm'){
    u <- myrmvnorm(n,u,sigma)
    return( u )
  }
  if(LIKE == 'mvnorm-multinom'){
    u <- myrmvnorm(n,u,sigma)
    zs <- rowSums(exp(u))
    z1 <- 1/(1 + zs)
    zm <- exp(u)/ (1 + zs)
    u  <- cbind(zm,z1)
    return( myrmultinom(size,u) )
  }
  numeric(0)
}
#########################################

simMVData <- function(LIKE,k,r){
	
   xnames <- c('intercept',paste('x',c(1:(k-1)),sep='-') )
   p      <- length(xnames)
   x      <- matrix(runif(n*k,-4,4),n,k)
   x[,1]  <- 1
   colnames(x) <- xnames
   
   b           <- matrix(runif(k*ns,-2,1),k,r)
   specnames   <- paste('S',c(1:r),sep='-')
   colnames(b) <- specnames
   rownames(b) <- xnames
   sigma <- diag(runif(r,0,.1/r),r)
   rownames(sigma) <- specnames
   colnames(sigma) <- specnames

   mu <- x %*% b
   y     <- myrmvnorm(n,mu,sigma)
   sigma <- crossprod(y - x%*%b)/(n-k)
   y     <- myrmvnorm(n,mu,sigma)
   colnames(y) <- specnames
   
   if(LIKE == 'mvnorm-Pois')    z <- matrix(rpois(length(y),exp(y)),n,r,byrow=F)
   if(LIKE == 'mvnorm-multinom'){
   	
   	  p1 <- matrix(rbeta(n*r,.2,.2),n,r)
   	  p1[cbind(c(1:n),sample(c(1:r),n,replace=T))] <- 10
   	  p1 <- p1/matrix(rowSums(p1),n,r)
   	  zz <- myrmultinom(10,p1)
   	  
   	  b  <- bInitMNomLogit(x,zz,size)
   	
   	  sigma <- sigma[1:(r-1),1:(r-1)]
     y  <- myrmvnorm(n, x %*% b,sigma)
   	  zs <- rowSums(exp(y))
     z1 <- 1/(1 + zs)
     zm <- exp(y)/ (1 + zs)
     y  <- cbind(zm,z1)
   	  z  <- myrmultinom(size,y)
   	  plot(y,jitter(z))
   	}
   
   list(x = x, b = b, y = y, z = z, s = sigma, specnames = specnames)
}

gibbsLoop <- function(LIKE,ng,x,y,b,priorB=b*0,priorVB=diag(1000,length(b)),
                      sigma = 0,z = numeric(0),burnin=1,
                      loB=NULL,hiB=NULL){
                      	
  tiny <- 1e-6
  
  k  <- ncol(x)
  r <- br <- 1                       #responses
  if(is.matrix(y)){
    r <- ncol(y)
    if(LIKE == 'multinom')br <- r - 1
  }
  n  <- nrow(x)
  kx <- length(b)
  if(!is.matrix(b))b <- as.matrix(b)
  
  priorIVB <- solve(priorVB)
  
  bgibbs <- matrix(NA,ng,kx)
  if(is.null(rownames(b))) rownames(b) <- paste('q',c(1:k),sep='-')
  if(is.null(colnames(b))) colnames(b) <- paste('s',c(1:br),sep='-')
  colnames(bgibbs) <- as.vector(outer(rownames(b),colnames(b),paste,sep='_') )

  sgibbs <- matrix(NA,ng,max(1,length(sigma)))
  

  
  if(sigma[1] == 0)sgibbs <- numeric(0) #variances
  
  sg <- sigma
  
  if(is.matrix(sigma)){                  #Wishart prior
    sg        <- prior.W
    colnames(sgibbs) <- outer(rownames(sigma),rownames(sigma),paste,sep='_') 
    colnames(bgibbs) <- as.vector(outer(colnames(x),colnames(b),paste,sep='_'))
  }
  
  pBVar <- solve(crossprod(x))

  pred <- pred2 <- rep(0,nrow(x)*r)

  bg <- b
  yg <- y
  y1 <- yg
  
  if(LIKE == 'pois') bg <- pBVar%*%crossprod(x,log(y + .1))
  if(LIKE == 'binom'){
  	yl <- y
  	yl[yl == 0] <- .1
  	yl[yl == 1] <- .9
  	bg <- pBVar%*%crossprod(x,logit(yl))
  }
  
  if(LIKE == "mvnorm-multinom"){
  	y1   <- prob2Logit(yg)
  }
  if(LIKE == 'multinom'){
  	size <- rowSums(y)
  	pBVar   <- diag(.0001,k*(r-1))
  	priorB  <- b*0
  	priorVB <- diag(100,k*(r-1))
  	priorIVB <- solve(priorVB)
  }
  
  dev <- 0
  np  <- 0
  
  pbar <- txtProgressBar(min=1,max=ng,style=1)

  for(g in 1:ng){

    bg <- bUpdateGibbs(LIKE,x,y1,bg,priorB,priorIVB,loB,hiB,sigma=sg,pBVar)
    
    if(sigma[1] > 0 & length(sigma) == 1){
    	 sg <- updateSigma(y1,x%*%bg)
    }
    
    if(length(grep('mvnorm',LIKE)) > 0){
    	 sinv <- wishsamp(x,y1,bg)
    	 sg   <- solve(sinv)
    }
    
    if(LIKE == "mvnorm-multinom"){
    	y1 <- ysampMvnormMultinom(x,y1,z,bg,sg)$y
    	yg <- logit2Prob(y1)
    }
    
    dev <- dev + sum(deviance(y1,x,bg,sg,LIKE))

    if(g > burnin){
      py    <- as.vector(simY(x,bg,LIKE,r,size,sg))
      pred  <- pred + py
      pred2 <- pred2 + py^2
      np    <- np + 1
    }

    bgibbs[g,] <- bg
    if(length(sgibbs) > 0)sgibbs[g,] <- sg
    
    if(g %in% c(100,200,500,1000)){
    	pBVar <- .5*cov(bgibbs[20:g,]) + diag(tiny,ncol(bgibbs))
    }
    
    setTxtProgressBar(pbar,g)
  }

  ymean <- pred/np
  yse   <- sqrt(pred2/np - ymean^2)
  
  if(r > 1){
  	ymean <- matrix(ymean,n,r)
  	yse   <- matrix(yse,n,r)
  }
  bmean <- colMeans(bgibbs)
  bmean <- matrix(bmean,nrow(bg),ncol(bg))
  smean <- numeric(0)
  if(length(sg) > 1) smean <- matrix(colMeans(sgibbs),nrow(sg),ncol(sg))
  if(length(sg) == 1)smean <- mean(sgibbs)
  
  meanDev <- dev/ng
  
  pd  <- meanDev - sum(deviance(y1,x,bmean,smean,LIKE))
  dic <- 2*pd + meanDev
  
  colnames(ymean) <- colnames(yse) <- colnames(y)

  list(bgibbs = bgibbs,sgibbs = sgibbs, ymean = ymean, yse = yse, dic = dic)
}
####################################################
simY <- function(x,b,LIKE,r = 1,size=rep(1,nrow(x)),sigma = 0,Effort = 1){     #simulate response

  u <- x%*%b

  if(LIKE == 'norm')return( rnorm(length(u),u,sqrt(sigma)) )
  if(LIKE == 'pois'){
  	 u <- exp(u)*Effort
  	 return( rpois(nrow(x),u) )
  } 
  if(LIKE == 'binom'){
  	 u <- invlogit(u)
  	 return( rbinom(nrow(x),1,u) )
  }
  if(LIKE == 'multinom'){
     zs <- apply(exp(u),1,sum)
     z1 <- 1/(1 + zs)
     zm <- exp(u)/ (1 + zs)
     u  <- cbind(zm,z1)
     return( myrmultinom(size,u) )
  }
  if(LIKE == 'mvnorm'){
    u <- myrmvnorm(n,u,sigma)
    return( u )
  }
  if(LIKE == 'mvnorm-multinom'){
    u <- myrmvnorm(n,u,sigma)
    zs <- apply(exp(u),1,sum)
    z1 <- 1/(1 + zs)
    zm <- exp(u)/ (1 + zs)
    u  <- cbind(zm,z1)
    return( myrmultinom(size,u) )
  }
  numeric(0)
}
##########################################################

sensIntercept <- function(bgibbs,xnames,ynames,hostNames){     

  #sensitivity coeffs for multinomial regression
  #bgibbs - multinomial logit par chains
  #xnames - length-k vector of names for input x 
  #ynames - length-(r-1) vector of names for response y

  r <- length(ynames)
  k <- length(xnames)

  fullnames <- as.vector(outer(xnames,ynames[-r],paste,sep='-')) #columns of bgibbs

  sbgibbs <- numeric(0)

  speccol <- matrix(unlist(strsplit(fullnames,'-')),ncol=2,byrow=T)[,2]
  trtcol  <- matrix(unlist(strsplit(fullnames,'-')),ncol=2,byrow=T)[,1]
  
  hostCols <- which(trtcol %in% hostNames)

  for(h in 1:(r-1)){
    wh1 <- which(speccol == ynames[h])
    wh2 <- wh1[ wh1 %in% hostCols ]
    wh0 <- wh1[ !wh1 %in% wh2 ]
    fmean  <- rowMeans(bgibbs[,wh2])          #mean over all hosts
    int    <- sweep(bgibbs[,wh2],1,fmean)     #host effects as departure from mean
    int    <- cbind(fmean,bgibbs[,wh0],int)
    colnames(int)[1] <- paste('int',ynames[h],sep='-')
    sbgibbs <- cbind(sbgibbs,int)
  }
  sbgibbs
}
#######################################################
multinomLike <- function(y,x,b){  #log likelihood multinomial logit
	
  tiny <- 1e-20
  huge <- 1 - tiny
  
     z2 <- logit2Prob(x%*%b)
  
     z2[z2 < tiny] <- tiny
     z2[z2 > huge] <- huge
     list(like = y*log(z2), theta = z2)
}
####################################################

bmultiProp <- function(r,k,b = matrix(0,k,r-1),pBVar=diag(.1,k*(r-1)),loB=NULL,hiB=NULL){  
	
    bvec <- as.vector(b)
    cvec <- myrmvnorm(1,t(bvec),pBVar)
    cc   <- matrix(cvec,nrow(b),ncol(b))
    
    if(!is.null(loB)){                         #if lob and hib available, use tnorm.mvt
    	cvec <- tnorm.mvt(bvec,bvec,pBVar,loB,hiB,times=4)
      cc   <- matrix(cvec,nrow(b),ncol(b))
    }

  list(cc = cc, cvec = cvec)
}

####################################################

bUpdateMNom <- function(x,y,b,priorB,priorIVB,
              pBVar=diag(.1,nrow(b)*(ncol(b)-1)),loB=NULL,hiB=NULL){

  bvec <- as.vector(b)
  priorVB <- solve(priorIVB)

  tmp  <- bmultiProp(ncol(b)+1,nrow(b),b,pBVar,loB,hiB)
  cc   <- tmp$cc
  cvec <- tmp$cvec
  
  pnow <- multinomLike(y,x,b)$like
  pnew <- multinomLike(y,x,cc)$like

  pnow <- sum(pnow) + mydmvnorm(bvec,as.vector(priorB),priorVB,log=T)
  pnew <- sum(pnew) + mydmvnorm(cvec,as.vector(priorB),priorVB,log=T)

  a <- exp(pnew - pnow)
  z <- runif(1,0,1)
  if(z < a)b <- cc
  b
}


updateWishart <- function(yy,predy,priorS,priorSdf,INVERSEONLY=F){
  
  n  <- nrow(yy)
  ss <- crossprod(yy - predy) + priorS*priorSdf
  df <- n + priorSdf
  
  testv <- try( solve(ss) )
  
  if(inherits(testv,'try-error')){
    diag(ss) <- diag(ss) + diag(ss)*.1
    testv       <- solve(ss)
  }
  sinv <- rwish(df,testv)
  
  if(INVERSEONLY)return(sinv)
  sig  <- solve(sinv)
  list(sigma = sig, sinv = sinv)
}

getScWishMat <- function(SS,df,delta,priorO){
  
  di    <- diag(1/diag(delta))
  tt    <- priorO + di%*%SS%*%di
  IO    <- rwish(df,solve(tt)) 
  O     <- solve(IO)
  sigma <- delta%*%O%*%delta
  
  list(omega = O, omegaInv = IO, sigma = sigma)
}
  
updateScInvWish <- function(yz,predy,delta=diag(1,ncol(yz)),
                            priorO=diag(1,ncol(yz)),priorOdf=(1+ncol(yz)),
                            priorDmu=rep(1,ncol(yz)),priorDvar=rep(10,ncol(yz)),
                            varLo=rep(1e-8,ncol(yz)),varHi=rep(10000,ncol(yz)) ){ 
  
  #varBound is permissible range of variances
  
  ww <- which(diag(delta) < 1e-5)
  if(length(ww) > 0)delta <- delta*10
  
#  ww <- which(diag(delta) > varBound[2])
#  if(length(ww) > 0)diag(delta)[ww] <- varBound[2]
  
  
  k  <- ncol(yz)
  n  <- nrow(yz)
  SS <- crossprod(yz - predy)
  df <- n + priorOdf
  
 # del <- diag(exp(priorDmu))  #prior
  
  tmp <- getScWishMat(SS,df,delta,priorO)
  IO  <- tmp$omegaInv
  O   <- tmp$omega
  sig <- tmp$sigma
  
  vars <- diag(sig)   # vars <- diag(SS)
  
  loDel <- sqrt( varLo/diag(O) )
  hiDel <- sqrt( varHi/diag(O) )
  
  ww <- which(loDel < 1e-6)
  if(length(ww) > 0)loDel[ww]  <- 1e-6
  ww <- which(hiDel < loDel)
  if(length(ww) > 0)hiDel[ww]  <- loDel[ww]*1.01
  
  psd <- diag(delta)/20
  
  propD <- diag( tnorm(k,loDel,hiDel,diag(delta),psd),k )
  
  ts <- IO*SS                
  t0 <- ts
  diag(t0) <- 0
  rs1 <- rowSums(t0/matrix(diag(delta),k,k,byrow=T))/diag(delta)
  rs2 <- rowSums(t0/matrix(diag(propD),k,k,byrow=T))/diag(propD)
  
  tsv <- diag(ts)/2/vars
  
  pnow <- -(n+1)*log(diag(delta)) - tsv - rs1 - ( (log(diag(delta)) - priorDmu)^2)/2/priorDvar
  pnew <- -(n+1)*log(diag(propD)) - tsv - rs2 - ( (log(diag(propD)) - priorDmu)^2)/2/priorDvar
  
  pnow <- sum(pnow)
  pnew <- sum(pnew)
  
  if( runif(1,0,1) < exp(pnew - pnow) ) diag(delta) <- diag(propD)
  
  sig <- delta%*%O%*%delta
  
  list(delta = delta, sigma = sig)
}


updateInvWish <- function(yy,predy,v=2,sigma,aa,A=rep(1000,ncol(yy))){   
  #Huang and Wand, Bayesian Analy (2013)
  
  k  <- ncol(yy)
  nn <- nrow(yy)
  SS <- crossprod(yy - predy)
  df <- nn + v + k
  
  sinv <- solve(sigma)
  
  aa <- 1/rgamma(k,(v + k)/2,v*diag(sinv) + 1/A^2)
  
  sinv  <- rwish(df,solve(SS + 2*v*diag(aa))) 
  smat  <- solve(sinv)
  
  list(sigma = smat, sinv = sinv, aa = aa)
}


################################################
bUpdateNorm <- function(xx,yy,b,
                priorB=matrix(ncol(xx)*0,ncol(xx)),priorIVB=diag(1/100,ncol(xx)),
                loB=NULL,hiB=NULL,sigma){

  V <- solve( crossprod(xx)/sigma + priorIVB )
  v <- crossprod(xx,yy)/sigma + priorIVB%*%priorB
  if( is.null(loB) & is.null(hiB) )return( t( myrmvnorm(1,t(V%*%v),V) ) )
  
  tnorm.mvt(V%*%v,V%*%v,V,loB,hiB)
}
 

###########################################33
bUpdateGibbs <- function(LIKE,x,y,b,priorB,priorIVB,
                         loB=NULL,hiB=NULL,sigma = 0,pBVar=0){

  if(LIKE == 'norm')    return( bUpdateNorm(x,y,b,priorB,priorIVB,loB,hiB,sigma) )
  if(LIKE == 'multinom')return( bUpdateMNom(x,y,b,priorB,priorIVB,pBVar,loB,hiB) )
  if(LIKE %in% c('mvnorm','mvnorm-multinom'))return( bUpdateMVNorm(x,y,b,sigma) )

  b <- matrix(b,length(b),1)
  if( is.null(loB)) c <- t(myrmvnorm(1,t(b),pBVar))   #proposal
  if(!is.null(loB)) c <- tnorm.mvt(b,b,pBVar,loB,hiB,times=1)

  znow <- x%*%b
  znew <- x%*%c

  if(LIKE == 'pois'){
     pnow <- dpois(y,exp(znow),log=T)
     pnew <- dpois(y,exp(znew),log=T)
  }
  if(LIKE == 'binom'){
     pnow <- dbinom(y,1,invlogit(znow),log=T)
     pnew <- dbinom(y,1,invlogit(znew),log=T)
  }

  priorVB <- solve(priorIVB)
  
  pnow <- sum(pnow) + mydmvnorm(t(b),priorB,priorVB,log=T)
  pnew <- sum(pnew) + mydmvnorm(t(c),priorB,priorVB,log=T)

  a <- exp(pnew - pnow)
  z <- runif(1,0,1)
  if(z < a)b <- c
  b
}
####################################################
updateSigma <- function(y,mu,s1=1,s2=1){

  u1 <- s1 + length(y)/2
  u2 <- s2 + .5*sum( (y - mu)^2 )
  1/rgamma(1,u1,u2)
}

######################################
wishsamp <- function(x,yy,b){   #sample from Inv Wishart

   r    <- ncol(b)
   scp  <- crossprod((yy - x %*% b))
   vmat <- solve(scp + prior.W*prior.WDF)
   v2   <- ns + prior.WDF
   stmp <- myrmvnorm(v2,matrix(0,v2,r),vmat)
   crossprod(stmp)
}
###################################
bUpdateMVNorm <- function(x,yy,bg,sigma,alpha=0,lo=NULL,hi=NULL,X=NULL){  # update b's for mvnorm

  r      <- ncol(yy)
  k      <- ncol(x)

  if(is.null(X))X <- crossprod(x)
  
  bigv   <- solve(X)      #multivariate scale invariate prior (Minka)
  smallv <- crossprod(x,yy)
                 
  mu     <- bigv%*%smallv/(alpha + 1)
  vaa    <- kronecker(sigma,bigv/(alpha + 1))   
  if(is.null(lo)) bg <- matrix( rmvnorm(1,as.vector(mu),sigma=vaa) ,k,r,byrow=F)
  if(!is.null(lo))bg <- matrix(tnorm.mvt(as.vector(mu),as.vector(mu),vaa,lo,hi,times=1),k,r)
  
  colnames(bg) <- colnames(yy)
  rownames(bg) <- colnames(x)
  bg
}

#################################

allomConvert <- function(xx,specCodes,
                         allomFile='/nfs/clark/clark.unix/allocationmodel/datafiles/allometricCoeffs.txt',
                         codeColumn='species',
                         int='bmassInt',slope='bmassSlope',
                         defaultSpec='other', invert=F){
  
  # x          - vector or matrix of values to convert (e.g., diam)
  # specCodes  - vector equal to nrow(x) 
  # codeColumn - column in allomFile containing specCodes
  # defaultSpec- value to use if specCodes missing from allomFile (e.g., 'other')
  #              there is a row in allomFile with defaultValue in codeColumn
  
  if(!is.matrix(xx))xx <- matrix(xx,ncol=1)
  if(!is.matrix(specCodes))specCodes <- matrix(specCodes,ncol=1)
  
  coeff  <- read.table(allomFile,header=T)
  smatch <- match(specCodes,coeff[,codeColumn])
  
  ws <- which(coeff[,'species'] == 'other')
  
  wmiss <- which(is.na(smatch))
  if(length(wmiss) > 0){
    ss <- sort( unique( specCodes[wmiss] ) )
    sc <- 'missing from allometry:'
    for(k in 1:length(ss))sc <- paste(sc,ss[k],sep=', ')
    message( sc )
    smatch[wmiss] <- which(coeff[,'species'] == 'other')
  }
  
  c1 <- coeff[smatch,int]
  c2 <- coeff[smatch,slope]
  c1[is.na(c1)] <- coeff[ws,int]
  c2[is.na(c2)] <- coeff[ws,slope]
  
  tmp <- 10^( c1 + c2*log10(xx) )
  if(invert)tmp <- 10^( (log10(xx) - c1 )/ c2 )
  tmp
}



diam2mass <- function(d,alloB,alloL){  #allo has c(int,slope) on log10 scale

  b <- 10^alloB[,1] *d^alloB[,2]
  l <- 10^alloL[,1] *d^alloL[,2]
  l[is.finite(l) & l > b] <- b[is.finite(l) & l > b]*.9
  cbind(b,l)

}

stemMass2diam <- function(stem,leaf,alloB){

  10^( (log10(stem + leaf) - alloB[,1])/alloB[,2] )
} 

mass2diam <- function(mass,allo){

  10^((log10(mass) - allo[,1])/allo[,2])
}


#################################
ysampMvnormMultinom <- function(x,y,z=NULL,b,sigma){  

  #sample y on MVlogit scale
  #z - counts, same dimensions as y

  r  <- ncol(y)
  k  <- nrow(b)
  n  <- nrow(y)
  
  propy <- matrix(rnorm(length(y),y,.01),n,r,byrow=F)
  
  pnow <- rep(0,n)
  pnew <- rep(0,n)
  
  for(i in 1:n){
    pnow[i] <- mydmvnorm(y[i,],(x %*% b)[i,], sigma,log=T)
    pnew[i] <- mydmvnorm(propy[i,],(x %*% b)[i,], sigma,log=T)
   }

  zs    <- rowSums(exp(y))
  z1    <- 1/(1 + zs)
  zm    <- exp(y)/ (1 + zs)
  znow  <- cbind(zm,z1)

  zs    <- rowSums(exp(propy))
  z1    <- 1/(1 + zs)
  zm    <- exp(propy)/ (1 + zs)
  znew  <- cbind(zm,z1)

  if(!is.null(z)){
    pnow <- pnow + z*log(znow)
    pnew <- pnew + z*log(znew)
  }

  pnow <- rowSums(pnow)
  pnew <- rowSums(pnew)

  a  <- exp(pnew - pnow)
  zz <- runif(length(a),0,1)
  y[zz < a,] <- propy[zz < a,]
  accept <- length(zz[zz < a])

  list(y = y, a = accept)
}

####################################################
getSens <- function(xv,bchain){

  #names(xv) are repeated for each class in bgibbs
  #ncol(bgibbs) = (r - 1)*length(xv)

  kk <- length(xv)

  nsim <- 2000
  sens <- matrix(NA,nsim,ncol(bchain))
  colnames(sens) <- colnames(bchain)

  for(j in 1:nsim){

    gj  <- sample(ng,1)
    bgj <- matrix(bchain[gj,],kk,r-1)
    sens[j,] <- as.vector( multinomSens(bgj,xvals=xv) )
  }
  sens
}

###############################################
plotSens <- function(sens,ytic,xnames,ynames,
                     SIGONLY=F){  #plot multinomial sensitivities

# sens - output from getSens

  colF <- colorRampPalette(c('darkblue','blue','green','yellow','orange','red'))

  if('int' %in% xnames)xnames <- xnames[xnames != 'int']

  r <- length(ynames)
  k <- length(xnames)

  fullnames <- as.vector(outer(xnames,ynames[-r],paste,sep='-')) #columns of bgibbs

  tmat <- processPars(sens)$summary 
  tt  <- matrix( unlist(strsplit(rownames(tmat),'-')),nrow(tmat),2,byrow=T )
  w0  <- grep('int',rownames(tmat))
  if(length(w0) > 0){
    tmat <- tmat[-w0,]
    tt   <- tt[-w0,]
  }

  par(mfrow=c(1,1),bty='n')

  xtic <- c(0:k)+.5

  plotSetup(xtic,ytic,xvals=rep(' ',length(xtic)),ylabel='Sensitivity')
  
  ylvalue <- ytic[1]+.05*diff(range(ytic))

  text(c(1:k),ylvalue,xnames,srt=90,pos=3,cex=1)

  specindex <- matrix(unlist(strsplit(fullnames,'-')),ncol=2,byrow=T)[,2]

  jseq <- seq(-1,1,length=r)*.2
  cols <- colF(r-1)
  
  keep <- numeric(0)
  plusMinus <- matrix(0,length(xnames),2)
  colnames(plusMinus) <- c('negative','positive')
  rownames(plusMinus) <- xnames

  for(j in 1:(r-1)){
    ww <- which(tt[,2] == ynames[j])
    tj <- tmat[ww,]
    tx <- tt[ww,1]
    ty <- tt[ww,2]
    wj <- which(tj[,'0.025'] < 0 & tj[,'0.975'] > 0)  # not different from zero
    wlo <- which(tj[,'0.975'] < 0)
    whi <- which(tj[,'0.025'] > 0)
    plusMinus[tx[wlo],1] <- plusMinus[tx[wlo],1] + 1
    plusMinus[tx[whi],2] <- plusMinus[tx[whi],2] + 1
    if(SIGONLY){
      if(length(wj) == nrow(tj))next
      if(length(wj) > 0){
        tj <- tj[-wj,]
        tx <- tx[-wj] 
        ty <- ty[-wj]
      }
    }
    if(!is.matrix(tj))tj <- matrix(tj,1)
    for(jk in 1:nrow(tj)){
      xv <- match(tx[jk],xnames)
      xx <- jseq[j] + xv
      lines( c(xx,xx), tj[jk,3:4],lwd=3,col=cols[j])
      points( xx,tj[jk,1],col=cols[j],pch=3)
    }
    keep <- c(keep,j)
  }	
  plo <- paste( round(100*plusMinus[,1]/rowSums(plusMinus),0),'%',sep='')
  phi <- paste( round(100*plusMinus[,2]/rowSums(plusMinus),0),'%',sep='')
  text(.3+c(1:k),ylvalue,plo,srt=90,cex=.8,pos=3)
  text(c(1:k),rep(max(ytic),k),phi,srt=90,cex=.8)
  
  
#  legend('topleft',ynames[keep],text.col=cols[keep],cex=.7,bty='n')
  dev.print(device=postscript,file='multinomModel.ps',width=6,horizontal=F)
  
  invisible( list(sensMat = tmat) )

}

#######################################################
multinomSens <- function(bb,x=NULL,xvals=NULL){  #sensitivity coefficients multinomial logit
	
  tiny <- 1e-20
  huge <- 1 - tiny
  kk   <- nrow(bb)
  rr   <- ncol(bb)

  if(is.null(xvals))xvals <- matrix(colMeans(x),1)

  z     <- xvals%*%bb
  zs    <- rowSums(exp(z))
  zm    <- exp(z)/ (1 + zs)
  theta <- matrix(zm,kk,rr,byrow=T)
  bb*theta*(1 - theta)
}
####################################################

reverseSignFIA <- function(xx,vname){
  if(vname %in% c('plotBA', 'xeric','deficit') ){
    xx <- -xx
  } else {
    xx <- 1 - xx
  }
  xx
}

par2sens <- function(pchain,vname,xrangej,reverseSign,iSymbol='X'){
  
  #parChain is a MCMC matrix of beta values 
  #iSymbol is interaction, either 'X' or ':'
  #none are factors
  
  
  wc <- grep(vname,colnames(pchain))       #all
  wi <- grep(iSymbol,colnames(pchain))     #interactions
  wm <- wc[!wc %in% wi]                    #main effect
  
  xi <- wi[wi %in% wc]
  
 # xi <- matrix( unlist( strsplit(colnames(pchain)[wi],iSymbol) ) , ncol=2,byrow=T)[,1]
  xx <- pchain[,wm]
  if(vname %in% reverseSign)xx <- reverseSignFIA(xx,vname)
 
 if(length(xi) > 0){
  for(ii in 1:length(xi)){
    ix <- pchain[,xi[ii]] * .5
    if(xi[ii] %in% reverseSign)ix <- reverseSignFIA(ix,vname)
    xx <- xx + ix
  }
 }
  
  #incorporate range
  xx <- xx/(xrangej[2,vname] - xrangej[1,vname])
  xx
}
##########################################

reverseSignBetaChains <- function(bchain,reverseSign=character(0),iSymbol='X'){
  
  # when sign is reversed, names in reverseSign
  # iSymbol in columns with interactions
  
  if(length(reverseSign) == 0)return(bchain)
  
  for(k in 1:length(reverseSign)){
    
    wk <- wn <- grep(reverseSign[k],colnames(bchain))
    
    nk <- matrix( unlist(strsplit(colnames(bchain)[wk],'_')),ncol=2,byrow=T)[,1]
    wi <- grep(iSymbol, nk)
    
    if(length(wi) > 0){
      
      wn <- wn[-wi]
      ii <- matrix( unlist(strsplit(nk[wi],iSymbol)),ncol=2,byrow=T)
      
      w1 <- which(ii[,1] %in% reverseSign & !ii[,2] %in% reverseSign)
      w2 <- which(ii[,2] %in% reverseSign & !ii[,1] %in% reverseSign)
      ww <- unique( c(w1,w2) )
      if(length(ww) > 0)bchain[,wk[wi[ww]]] <- -bchain[,wk[wi[ww]]]
    }
    bchain[,wn] <- -bchain[,wn]
  }
  
  bchain
}


plotDensColumns <- function(var2plot=NULL,chain2plot,allChains=chainList,chainLength=1,
                            vnames=xnames,MAINEFFECT=T, intSign = ':',xtic=NULL,labs=NULL,
                            htFraction=.5,textSize=1,ORD=NULL,COLCODE=NULL,
                            reverseSign=character(0),omitSpecs=character(0),FILL=F, nd=512){
  
  #if !MAINEFFECT then interactions added to main effect
  
  quadCols <- numeric(0)

  wvar2 <- character(0)
  wvar <- var2plot
  wd   <- grep( paste(var2plot,'-',sep='') ,vnames)  # a climate variable with months
  if(length(wd) > 0)wvar <- vnames[wd[1]]          # if so, new name
  
  reverseM <- F
  if(wvar %in% reverseSign)reverseM <- T
  hasQuad <- hasInt <- F
  isInt <- F

  xx <- grep(intSign,vnames)           #some are interactions
  x2 <- grep('2',vnames)           #some are quadratic
  xi <- grep(intSign,var2plot)         #variable is an interaction

  wint <- grep(var2plot,vnames[xx])  #variable to plot has interaction terms
  w2   <- grep(var2plot,vnames[x2])  #variable to plot has quadratic terms

  if(length(wint) > 0)hasInt  <- T
  if(length(w2) > 0)  hasQuad <- T
  if(length(xi) > 0)  isInt   <- T
  if(isInt)hasInt <- hasQuad <- F

  
  xtAll <- ytAll <- numeric(0)
  ymax  <- numeric(0)
  xrange <- numeric(0)
  
  modeBySpec <- numeric(0)
  allSpecs   <- character(0)
  
  for(jj in 1:chainLength){
    
    ww <- grep('Chain',names(allChains))
    if(length(ww) > 0){
      chainj <- allChains
    } else {
      chainj <- allChains[[jj]]
    }
    
    namej  <- names(chainj)
    namej  <- unlist( strsplit(namej,'Chain') )
    
    wc     <- which(namej == chain2plot)
    mat    <- chainj[[wc]]
    cnames <- colnames(mat)
    cm     <- cnames[grep(wvar,cnames)]   #names for main effects
    
    if(length(cm) == 0){
      xtAll <- append(xtAll,list(numeric(0)))
      ytAll <- append(ytAll,list(numeric(0)))
      modeBySpec <- append(modeBySpec,list(numeric(0)))
      next
    }
    
    if(hasQuad){
      wvar2 <- paste(var2plot,'2',sep='')
      ww    <- sort( grep(wvar,cnames) )
      cm    <- cnames[ww]
    }
    mat    <- mat[,cm]
    cm     <- colnames(mat)
    
    if(hasInt){
      wvar2 <- paste(var2plot,'2',sep='')
      ww    <- sort( grep(wvar,cnames) )
      cm    <- cnames[ww]
    }
      
    
    if(length(omitSpecs) > 0){
      woc <- numeric(0)
      for(kk in 1:length(omitSpecs)){
        woc <- c(woc,grep(omitSpecs[kk],colnames(mat)))
      }
      mat <- mat[,-woc]
      cm  <- colnames(mat)
    }
    
    intCols  <- sort( unique( c( grep( paste(wvar,intSign,sep=''),cm), 
                                 grep( paste(intSign,wvar,sep=''),cm ) ) ) )
    if(length(wvar2) > 0)quadCols <- grep( wvar2,cm)
    mainCols <- grep(wvar,cm)
    mainCols <- mainCols[!mainCols %in% c(intCols,quadCols)]
    matM    <- mat[,mainCols]  # only main effects
    
    # if(reverseM)matM <- -matM
    
    if(isInt){
      xxn <- unlist(strsplit(var2plot,intSign)) 
    }
    
    if(!MAINEFFECT & hasInt){      # full effect: add interactions
      wi <- xx[wint]
      for(j in wi){
        xi <- unlist( strsplit(vnames[j],intSign) )
        viname <- xi[xi != wvar] 
        
        ii <- 1
        if(viname %in% reverseSign)ii <- -1 
        icol  <- grep(vnames[j],cm)
        if(length(icol) > 0)matM <- matM + mat[,icol]*.5*ii # assume interaction term is at mean value (.5)
      }
    }
    if(hasQuad){
      wi <- x2[w2]
      for(j in wi){
        qcol  <- grep(vnames[j],cm)
        matM <- matM + mat[,qcol]*2*.5
      }
    }
    
    mat <- matM
    cm <- colnames(mat)
    
  #  if(is.null(labs))
    labs  <- matrix( unlist( strsplit(cm,'_') ),ncol=2,byrow=2)[,2]
    nj    <- length(labs)
    
    rr <- range(matM,na.rm=T)
    if(!is.finite(rr[1]) & !is.finite(rr[2])){
      xtAll <- append(xtAll,list(numeric(0)))
      ytAll <- append(ytAll,list(numeric(0)))
      modeBySpec <- append(modeBySpec,list( numeric(0) ))
      next
    }
    
    tmp    <- chains2density(matM,labs=labs,reverseM=reverseM, nd=nd)   
    xt     <- tmp$x
    yt     <- tmp$y
    ymax   <- max( c(ymax,quantile(yt,.8)) ) 
    xrange <- range( c(xrange,tmp$xrange))
    
    wmax <-  apply(yt,1,which.max )
    xmax <- xt[ cbind(1:nrow(xt),wmax) ]
    names(xmax) <- names(wmax)
    modeBySpec <- append(modeBySpec,list( xmax ))
    
  #  rownames(xt) <- rownames(yt) <- paste('r',jj,rownames(xt),sep='-')
    xtAll <- append(xtAll,list(xt))
    ytAll <- append(ytAll,list(yt))
    
    allSpecs <- sort(unique(c(allSpecs,names(xmax))))
  }
  
  
  minMaxSpec <- matrix(NA,length(allSpecs),2)
  rownames(minMaxSpec) <- allSpecs
  
  for(s in 1:length(allSpecs)){
    
    minmax <- c(Inf,-Inf)
    
    for(j in 1:chainLength){
      
      wsj <- which( names(modeBySpec[[j]]) == allSpecs[s] )
      if(length(wsj) == 0)next
      
      if( minmax[1] > modeBySpec[[j]][wsj] )minmax[1] <- modeBySpec[[j]][wsj]
      if( minmax[2] < modeBySpec[[j]][wsj] )minmax[2] <- modeBySpec[[j]][wsj]
    }
    minMaxSpec[s,] <- minmax
  }
    
  
  if(is.null(xtic)){
    sc     <- diff(xrange)
    x1     <- xrange[1] - sc/3
    x2     <- xrange[2] + sc/3
    xtic   <- signif(seq(xrange[1],xrange[2],length=4),3)
  }
  
  if(is.null(ORD)){
    ordSpec <- order( minMaxSpec[,2] )   #order by minimum value
    if(min(xtic) <= 0 & abs(xtic[1]) > max(abs(xtic[-1])))ordSpec <- order( minMaxSpec[,1] )
    ORD <- allSpecs[ordSpec]
  }
  
  if( length(xtAll[[1]]) == 0 ){
    xtmp <- ytmp <- numeric(0)
    for(j in 2:length(xtAll)){
      xtmp <- append(xtmp,list(xtAll[[j]]))
      ytmp <- append(ytmp,list(ytAll[[j]]))
    }
    xtAll <- xtmp
    ytAll <- ytmp
    COLCODE <- COLCODE[-1]
  }
                     
  mainTitle <- paste( demLabel(chain2plot),var2plot )
  
  tmp <- plotPosteriorOrder(xx=xtAll[[1]],yy=ytAll[[1]],x1=xtAll,y1=ytAll,
                            xtic=xtic,
                            ymax=ymax,textSize=textSize,ORD=ORD,COLCODE=COLCODE,
                            MAIN=mainTitle,htFraction=htFraction,
                            FILL=FILL, nd=nd)

  #  tmp <- plotPosteriorOrder(xx=xt,yy=yt,xtic=xtic,
#                          ymax=ymax,textSize=textSize,ORD=ORD,COLCODE=COLCODE,
#                          MAIN=paste(chain2plot,var2plot),htFraction=htFraction,FILL=FILL)
  invisible(tmp)

}

plotChainDensity <- function(cMat,vnames=NULL,ncolPlot=2,xlab=' ',ylab=' ',
                             LOG=F,inColor=F){
  
  if(is.null(vnames)){
    mp <- 1
  } else {
    mp  <- length(vnames)/ncolPlot
  }
  nrr <- floor(mp) + 1
  
  graphics.off()
  par(mfrow=c(nrr,ncolPlot),bty='n')
  
  knames <- vnames
  if(is.null(vnames))knames <- 'all'
  
  for(k in 1:length(knames)){
    
    kk <- knames[k]
    if(is.null(vnames))kk <- NULL
    
    if(!is.null(kk)){
      gg <- grep(kk,colnames(cMat))
      if(length(gg) == 0)next
    }
    
    tmp <- chains2density(cMat,varName=kk)
    xt  <- tmp$x
    yt  <- tmp$y
    chainMat <- tmp$chainMat
    
    cols <- rep('black',nrow(xt))
    if(inColor){
      colF   <- colorRampPalette(c('black','brown','orange'))
      cols <- colF(nrow(xt))
    }
    
    nn <- nrow(chainMat)
    
    xrange <- quantile(xt,c(.05,.95))
    xlim <- range(xt)
    if(!LOG) plot(10,10,xlim=xlim,ylim=c(0,1.8*max(yt)),xlab=xlab,ylab=ylab,cex=.01)
    if(LOG)  plot(10,10,xlim=xlim,ylim=c(0,1.8*max(yt)),xlab=xlab,ylab=ylab,cex=.01,log='x')
    
    plotLabel(kk,above=T)
    
    j1 <- 1
    if(knames[1] == 'intercept')j1 <- 2
    
    for(j in j1:nrow(xt)){
      
      cj <- cumsum(yt[j,])
      cj <- cj/max(cj)
      ci <- xt[j, findInterval(c(.02,.98),cj) ]
      
      label <- rownames(xt)[j]
      
      wm <- which.max(yt[j,])
      lines(xt[j,],yt[j,],lwd=2)
      lines(range(xt[j,]),c(0,0),col=cols[j],lwd=2)
      if(ncol(chainMat) < 25)text(xt[j,wm],1.1*yt[j,wm],label,srt=55,pos=4,cex=1)
    }
  }
}
  
chains2density <- function(chainMat,labs=NULL,reverseM=F,varName=NULL, nd=512){
  
  #assumes column names are varName or 'something_varname'
  
  #chainMat - MCMC output [samples,chains]
  
  chNames <- colnames(chainMat)
  
  if(!is.null(varName)){
    wc <- grep(varName,colnames(chainMat))
    if(length(wc) == 0)stop('varName not found in colnames(chainMat)')
    
    ww <- grep('_',colnames(chainMat))
    if(length(ww) > 0){
      tmp <- matrix( unlist(strsplit(colnames(chainMat),'_')),ncol=2,byrow=T) 
      wc <- which(tmp[,2] == varName)
      if(length(wc) == 0)wc <- which(tmp[,1] == varName)
    }
    chainMat <- chainMat[,wc]
    if(!is.matrix(chainMat))chainMat <- matrix(chainMat,ncol=1)
    colnames(chainMat) <- chNames[wc]
  }
  
  nj <- ncol(chainMat)
  
  clab <- colnames(chainMat)
  if(is.null(labs) & !is.null(clab))labs <- clab
  if(is.null(labs) & is.null(clab)) labs <- paste('v',c(1:nj),sep='-')
  
  xt <- yt <- matrix(NA,nj,nd)
  rownames(xt) <- rownames(yt) <- labs
  
  xrange <- signif(range(chainMat),2)
  
  for(j in 1:nj){
    
 #   lj  <- labs[j]
    xj  <- chainMat[,j]
    tmp <- density(xj,n = nd, cut=0, na.rm=T)
    xt[j,]  <- tmp$x
    yt[j,]  <- tmp$y
    
  }
  yymax <- max(yt,na.rm=T)
  
  if(reverseM){
    xt <- -t( apply(xt,1,rev) )
    yt <- t( apply(yt,1,rev) )
  }
  
  list(x = xt, y = yt, xrange = xrange, ymax = yymax, chainMat = chainMat)
}


makeMap <- function( xx = mastData[,'lon'],yy = mastData[,'lat'], z,
                     mapx = c(-135, -55), mapy = c(25, 50), zlevs = NULL,
                     ngrid = 80, maskx=.12, masky=.12, whiteOutDist=.2,
                     units = ' ', POWER = 1, posValues = F, 
                     xlegend = NULL, ylegend = NULL, 
                     COLS = 'other', REVCOL = F, add = F,
                     zquantile = .999, STATES=T, LEGEND = T, cex = .9,
                     notMask = NULL, ZEROLINE = F, ONZERO = T ){
  
  require(mapdata)
  require(maptools)
  require(ggplot2)
  require(rgeos)
  
  wl <- which(xx >= mapx[1] & xx <= mapx[2] &
                yy >= mapy[1] & yy <= mapy[2] & 
                is.finite(z))
  xx <- xx[wl]
  yy <- yy[wl]
  z  <- z[wl]
  
  # exclude oceans
  elev <- subETOPO5(xx, yy)
  
  lgrid <- as.matrix( expand.grid(elev$x,elev$y) )
  wo    <- which(elev$z < 0, arr.ind=T)
  omask <- cbind( elev$x[wo[,1]],elev$y[wo[,2]]  )
  
  # plot(elev$x[wo[,1]],elev$y[wo[,2]], cex=.1)
  
  kk <- RANN::nn2( omask, cbind(xx, yy),  k = 1)$nn.dists # far to an ocean point
  wl <- which(kk > .05)
  xx <- xx[wl]
  yy <- yy[wl]
  z  <- z[wl]
  
  
  z  <- z^POWER
  zz <- z
  
  # erratic high values
  
  zq <- c(1 - zquantile, zquantile)
  
  if(posValues){
    w0 <- which(zz > 0)
    qz <- quantile(zz[w0], zquantile)
    zz[zz > qz] <- qz
  }else{
    qz <- quantile(zz, zq)
    zz[zz < qz[1]] <- qz[1]
    zz[zz > qz[2]] <- qz[2]
  }
  
  if( is.null(zlevs) ){
    
    zlevs <- quantile(zz, seq(0, 1, length=50))
    
    if(posValues)zlevs <- zlevs[zlevs > 0]
    
    zlevs <- sort(unique(zlevs))
    
    if(length(zlevs) < 10){
      zlevs <- quantile(zz, seq(0, 1, length=1000))
      zlevs <- sort(unique(zlevs))
    }
    
    while(length(zlevs) > 40){
      zlevs <- zlevs[seq(1, length(zlevs), by=2)]
    }
    if( max(zlevs) <= max(zz) )zz[zz > max(zlevs)] <- max(zlevs)
    if( min(zlevs) >= min(zz) )zz[zz < min(zlevs)] <- min(zlevs)
    if( posValues & min(zlevs) < 0 )zlevs <- c(0, zlevs[zlevs > 0])
  }
  
  nz <- length(zlevs)
  
  if(COLS == 'greens'){
    cseq <- c('#edf8fb','#b2e2e2','#66c2a4','#2ca25f','#006d2c')
  }else{
    cseq <- rev( c('#543005','#8c510a','#bf812d','#dfc27d','#f6e8c3','#f5f5f5',  # odd number needed for white
                   '#c7eae5','#80cdc1','#35978f','#01665e','#003c30') )
  }
  if(REVCOL)cseq <- rev(cseq)
  
  
  colM <- colorRampPalette( cseq ) 
  cols <- colM(nz)
  
  kk <- which( zlevs < 0 )
  if(length(kk) == 0 )ONZERO <- F
  
  if(ONZERO){  # white is on zero

    w0 <- max( which(zlevs < 0) )
    ff <- w0/nz
    hf <- nz/2
    
    if(ff < .5){
      start <- round( hf*(1 - ff) )
      cols  <- cols[start:nz]
    }
    colM <- colorRampPalette( cols ) 
  }
  
  
  cols <- colM(nz)
  
  tcols <- .colorSequence(slim = range(zlevs), colorGrad=colM, ncol=nz)
  # cols  <- tcols$cols
  scale <- tcols$scale
  
  z2 <- max(zlevs) + diff(zlevs)[nz-1]
  zlevs <- c(zlevs, z2)
  cols <- c(cols, cols[nz])
  nzz  <- length(zlevs) - 1
  
  maps::map(boundary=T, col = .getColor('tan', .5),
            lwd=2,xlim=mapx, ylim=mapy, add = add)
  surf <- values2contour(xx=xx,yy=yy, z=zz,
                         xlim = mapx, ylim = mapy,
                         nx=ngrid,ny=ngrid,col=cols[1:nzz], zlevs=zlevs,
                         lwd=.1,add=T,fill=T)
  xycol <- findInterval(zz, zlevs)
  
  if(ZEROLINE)
    contour(surf, levels= 0,lwd= 1,lty= 1 ,col= 'brown',add=T,
            drawlabels=F)
  
  
  if( !is.null(notMask) ){
    xm <- seq(notMask[1,1],notMask[1,2], length=20)
    ym <- seq(notMask[2,1],notMask[2,2], length=20)
    xg <- expand.grid(xm,ym)
    xx <- c(xx, jitter(xg[,1], 40) )
    yy <- c(yy, jitter(xg[,2], 40) )
  }
  
  csxy <- maps::map(xlim = mapx, ylim = mapy, plot = F)
  cx <- csxy$x
  cy <- csxy$y
  cw <- which(is.finite(cx) & is.finite(cy) )
  
  mapMask( xx = c(xx, cx[cw]), yy = c(yy, cy[cw]), 
           dx = maskx, dy = masky, whiteOutDist, col='white')
  
  maps::map(boundary=T,add=T,interior=F,col='white',lwd=4)
  maps::map(boundary=T,add=T,interior=F, col = .getColor('tan', .5),lwd=2)
  maps::map(boundary=T,add=T,interior=T, col = .getColor('tan', .5),lwd=1)
  if(STATES)maps::map('state',boundary=T,add=T,interior=T,col='tan',lwd=1.5)
  
  if(LEGEND){
    lseq  <- seq(1,length(zlevs)-1, length=5)
    ytick <- zlevs[lseq]      
    
    ytry <- signif(ytick^(1/POWER),1)
    
    if( length(unique(ytry)) < length(ytick) )ytry <- signif(ytick^(1/POWER),2)
    ytick <- ytry
    
    
    if(is.null(xlegend)){
      dx <- .02*diff(par('xaxp')[1:2])
      xlegend <- par('xaxp')[1] + dx*c(1.5, 3)
    }
    
    if(is.null(ylegend)){
      dy <- .5*diff(par('yaxp')[1:2])
      ylegend <- par('yaxp')[1] + dy*c(.35, .85)
    }
    
    print( c(xlegend, ylegend) )
    
    dy <- diff(ylegend)
    
    colorLegend(xlegend, ylegend, ytick=ytick,
                scale=zlevs[lseq],cols=cols[lseq],labside='right', cex = cex)
    text( mean(xlegend), ylegend[2] + .33*dy, units, pos=2, cex = cex)
  }
  
  invisible( surf )
}



#######################################
plotPosteriorOrder <- function(xx,yy, x1 = NULL, y1 = NULL, 
                               xtic=NULL,ymax=NULL,xlab=' ',
                               ALTLABEL=F,textSize=1,ORD=NULL,
                               COLCODE=NULL,
                               MAIN=NULL,htFraction=1,
                               X5=F,trimx=0,CI=.95,FILL=F, nd=512){    

  #  yy - each row is a density
  #  xx - each row is the value for the density in y
  #  x1, y1 - lists of additional densities to be plotted with x,y
  #  names are rownames for y
  #  COLCODE  - names used to match with y
  #  htFraction < 1 makes overlap with next density
  #  X5 - plot 5x exaggeration
  #  trimx/2 is fraction of tail to trimx
  
  wsord   <- apply(yy,1,which.max)
  word    <- order( xx[cbind(1:nrow(xx),wsord)] )
  specOrd <- rownames(yy)[word]
  
  if(length(y1) > 0){
    
    wxx <- numeric(0)
    
    for(j in 1:length(y1)){
      
      xj      <- x1[[j]]
      yj      <- y1[[j]]
      
      wsx   <- matrix(xj[apply(yj,1,which.max)],nrow=1)
      colnames(wsx) <- rownames(xj)
      
      wxx <- appendMatrix(wxx,wsx,SORT=T)
    }
    
    minBySpec <- apply(wxx,2,min,na.rm=T)
    word      <- order(minBySpec)
    specOrd   <- names(minBySpec)[word]
      
  }
  
  specOrd <- rev(unique(specOrd))
  yvalw   <- seq(0,ymax,length=3)
  nc      <- length(specOrd)
  
  xrw     <- signif(range(xx),1)
  yrw     <- round(htFraction*quantile(yy,.95),-1)
  if(yrw == 0)yrw <- signif(quantile(yy,.98),1)

  
  
  if(!is.null(ORD)){
    specOrd <- ORD
    word    <- match(specOrd,rownames(yy))
  }
  
  tail1 <- (1 - CI)/2
  tail2 <- 1 - tail1
  
  postInt <- matrix(0,length(specOrd),3)
  rownames(postInt) <- specOrd
  colnames(postInt) <- c(tail1,tail2,'neg/pos')

 # if(!is.null(x1)){
 #   nc      <- nc + nrow(x1)
 #   wnord   <- apply(y1,1,which.max)
 #   xrw     <- range(c(xx,x1))
 #   yrw     <- round(htFraction*max(rbind(yy,y1)),-2) 
 #   word    <- order( c(xx[cbind(1:nrow(xx),wsord)],x1[cbind(1:nrow(x1),wnord)]) )
 #   specOrd <- c(rownames(yy),rownames(y1))[word]
 # }

 # totalHt <- yrw*nc
  
#  yrw     <- htFraction*ymax
  totalHt <- yrw*nc

  postext <- 4

  if(is.null(ymax))ymax <- yrw


  if(!is.null(COLCODE)){
  #  colSpecs <- COLCODE[specOrd]
    colSpecs <- COLCODE[word]
    if(length(x1) > 0){
      repSpecs <- COLCODE[1:length(x1)]
      colSpecs <- rep(COLCODE[1],length(x1))
    }
  }
  if(is.null(COLCODE)){
    colS     <-  colorRampPalette( c("darkblue","orange") )
    colSpecs <- colS(max(1,length(y1)))
    repSpecs <- colS(max(1,length(y1)))
  }
#  print(colSpecs)
  ytic <- seq(0,totalHt,by=yrw)
  if(length(ytic) < (nc-1))ytic <- c(ytic,max(ytic) + diff(ytic)[1])

  if(is.null(xtic)){
    xinc <- signif(diff(xrw)/5,1)
    xtic <- signif( seq(xrw[1],xrw[2],by=xinc),1)
  }
  jj   <- rev(ytic)[-1]

################## color box if != 0
  plotSetup(xtic,ytic,xlabel=xlab,yvals=rep(' ',length(ytic)),
            ylabel='Density',lcVer='grey')
  
  
  if(length(y1) > 1){
    postInt <- cbind(postInt, matrix(0,nrow(postInt),(length(y1)-1)) )
  }
  
  for(j in 1:length(specOrd)){

    w0 <- which(rownames(yy) == specOrd[j])
    w1 <- numeric(0)
    
    for(jk in 1:length(y1)){
      ww <- which(rownames(y1[[jk]]) == specOrd[j])
      if(length(ww) == 0)ww <- 0
      w1 <- c(w1,ww)
    }
    
    if(length(c(w0,w1)) == 0)next
    
    if(length(w0) > 0){
      
      dx <- diff(xx[w0,])[1]
      cj <- cumsum(dx*yy[w0,])
      if(cj[1] > 0)cj[1] <- 0
      if(cj[length(cj)] < 1)cj[length(cj)] <- 1    #cdf
      rj <- xx[w0, findInterval(c(trimx/2,1 - trimx/2),cj) ]   #determine x limits
      yj <- jj[j] + yy[w0,]
      yj[yj > ymax] <- ymax
      
      postInt[j,1:2] <- xx[w0,findInterval( c(tail1,tail2),cj)]
      if(postInt[j,1] > 0) postInt[j,3] <- 1
      if(postInt[j,2] < 0) postInt[j,3] <- -1
    }
    
    if(length(y1) > 1){
      
      for(jk in 2:length(y1)){

        if(w1[jk] == 0)next
        
        xj <- x1[[jk]][w1[jk],]
        dx <- diff(xj)[1]
        cjj <- cumsum(dx*y1[[jk]][w1[jk],])
        if(cjj[1] > 0)cjj[1] <- 0
        if(cjj[length(cjj)] < 1)cjj[length(cjj)] <- 1    #cdf
        rjj <- xj[ findInterval(c(tail1,tail2),cjj) ] 
        if(rjj[1] > 0) postInt[j,2 + jk] <- 1
        if(rjj[2] < 0) postInt[j,2 + jk] <- -1
      }               
    }                    
    
    anyNeg <- anyPos <- F
    if(max(postInt[j,-c(1:2)]) > 0)anyPos <- T
    if(min(postInt[j,-c(1:2)]) < 0)anyNeg <- T

    jmax <- totalHt
    if(j > 1)jmax <- jj[j-1]
    if(anyPos)rect(max(c(0,xtic[1])),jj[j],max(xtic),jmax,col='wheat',border='wheat')
    if(anyNeg)rect(xtic[1],jj[j],min(c(0,max(xtic))),jmax,col='wheat',border='wheat')
    if(anyNeg & anyPos)rect(xtic[1],jj[j],max(xtic),jmax,col='wheat',border='wheat')
  }

  abline(v=0,col='grey',lwd=3)
###########################################

  ldens <- 20

  for(j in 1:length(specOrd)){

    w0 <- which(rownames(yy) == specOrd[j])
    w1 <- numeric(0)
    
    for(jk in 1:length(y1)){
      ww <- which(rownames(y1[[jk]]) == specOrd[j])
      if(length(ww) == 0)ww <- 0
      w1 <- c(w1,ww)
    }

    if(length(c(w0,w1)) == 0)next
  
    xl <- numeric(0)

    if(length(w0) > 0){

       xl <- range(xx[w0,])
       yj <- yy[w0,]
       yj[1] <- yj[length(yj)] <- 0
       if(X5){
         yk <- 5*yj
         yk[yk > ymax] <- ymax
         yk <- jj[j] + yk
       }
       yj[yj > yrw] <- yrw
       yj <- jj[j] + yj

       lines(xx[w0,],yj,lwd=2,col=colSpecs[j])

       if(FILL){
         polygon(c(xx[w0,],xx[w0,1]),c(yj,yj[1]),
                 border=repSpecs[1],col=colSpecs[j])
       }
       if(X5)polygon(c(xx[w0,],xx[w0,1]),c(yk,yk[1]),border=colSpecs[j],
              col=colSpecs[j],density=ldens)
    }
    
    if(length(w1) > 0){
      
      for(jk in 1:length(y1)){
        
        if(w1[jk] == 0)next
        x2 <- range(x1[[jk]][w1[jk],])
        
    #    if(diff(x2) < .1){
          
          yj <- y1[[jk]][w1[jk],]
          xj <- x1[[jk]][w1[jk],]
          yj[yj > yrw] <- yrw
          yj <- jj[j] + yj
        
     #     lines(xj,yj,lwd=5,col='white')
          lines(xj,yj,lwd=2,col=repSpecs[jk])
          if(FILL){
            polygon( c(xj[1],xj,xj[length(xj)], xj[1]) , c(jj[j],yj,jj[j],jj[j]),
                     border=repSpecs[jk],col=repSpecs[jk])
          }
          xl <- range(c(xl,x2))
  #      }
      }
    }
  #  lines(xl,c(jj[j],jj[j]),col=colSpecs[j],lty=2,lwd=2)

    yt <- htFraction*max(yy[w0,])
    xt <- xx[w0,which.max(yy[w0,]) - 12]
    xxl <- max(xl)
    
    leftDist <- min(xl) - min(xtic)
    rightDist <- max(xtic) - max(xl)
    
    LEFT <- T
    
    if(leftDist < 0)LEFT <- F
    if(leftDist > 0 & leftDist > rightDist)LEFT <- T
    if(rightDist > 0 & rightDist > leftDist)LEFT <- F
    

    if(!ALTLABEL){
      if(LEFT){
        xxl <- min(xl)
        postext <- 2
      }
      if(!LEFT){
        postext <- 4
        xxl <- max(xl)
      }
    }

    if(xxl > max(xtic) & !ALTLABEL){
        xxl <- min(xl)
        postext <- 2
    }

   if(postext == 4)xxl <- max(xl)
   if(postext == 2)xxl <- min(xl)

   text(xxl,jj[j],specOrd[j],col=colSpecs[j],pos=postext,cex=textSize)

   if(ALTLABEL & postext == 4){
     postext <- 2
   } else {postext <- 4}

  }

  ys2 <- yvalw[length(yvalw)]
  
  if(!is.null(MAIN))title(MAIN)
  
  
  tk    <- postInt[,-c(1:2)]
  
  if(is.matrix(tk)){
    sump  <- which(tk == 1,arr.ind=T)[,2]
    sumn  <- which(tk == -1,arr.ind=T)[,2]
    sump  <- round(byIndex(sump*0+1,sump,sum,coerce=T)/nrow(tk),3)*100
    sumn  <- round(byIndex(sumn*0+1,sumn,sum,coerce=T)/nrow(tk),3)*100
  }
  if(!is.matrix(tk)){
    sump  <- which(tk == 1,arr.ind=T)
    sumn  <- which(tk == -1,arr.ind=T)
    sump  <- round( length(sump)/length(tk)*100,1)
    sumn  <- round( length(sumn)/length(tk)*100,1)
  }
  
  sump <- paste(sump,'%',sep='')
  sumn <- paste(sumn,'%',sep='')
  
  list(specOrd=specOrd, postInt = postInt, repSpecs = repSpecs, percPos = sump, percNeg = sumn)

}

getScoreNorm <- function(x,mu,xvar){  #Gneiting and Raftery's proper scoring rule

  #outcome x, prediction mean variance (mu, xvar)

  - ( (x - mu)^2)/xvar - log(xvar)

}

getScoreBinary <- function(x,p){    #binary outcome x (0,1), probability p

  (log(p))^x*(log(1 - p))^(1 - x)

}

advection <- function(u,vel,dt,dx){
	
	#finite difference left to right with velocity vel
	#Crank-Nicholson, http://farside.ph.utexas.edu/teaching/329/lectures/node92.html
	#u - current concentration
	
	n <- length(u)
	const <- vel*dt/dx
	w <- rep(0,n)
	
	b <- rep(1,n)
   a <- b*(.25*const)
   c <- -a
   
   for(j in 2:(n-1))w[j] <- u[j] - .25*const*(u[j+1] - u[j-1])
   
   w <- triDiagonal(b,a,c,w)
   w[1] <- w[n] <- 0
   w[w < 0] <- 0
   w
}
   
triDiagonal <- function(a,b,c,y){
	
  #  solves Ax = y for tridiagonal matrix A
  #  a       main diagonal 
  #  b       upper diagonal 
  #  c       lower diagonal
  #  y       right-hand side vector

  n <- length(a)

  #   factorization
  for (i in 1:(n-1)){
     b[i]   <- b[i] / a[i]
     a[i+1] <- a[i+1] - c[i] * b[i]
  }

  #   forward substitution
  y[1] <- y[1]/a[1]
  for(i in (2:n))y[i] <- ( y[i] - c[i-1] * y[i-1] ) / a[i]

  #   back substitution
  for(i in (n-1):1)y[i] <- y[i] - y[i+1] * b[i]
  
  y
}

mcmcSetup <- function(chains=NULL,sums=NULL,chainVals = NULL, 
                      sumVals = NULL,ng=1000){

  chainList <- chainSum <- numeric(0)

  if(!is.null(chainVals))for(k in 1:length(chains))assign(chains[k],chainVals[[k]])
  if(!is.null(sumVals))  for(k in 1:length(sums))  assign(sums[k],sumVals[[k]])

  if(!is.null(chains)){

      ne <- length(chains)
      for(k in 1:ne){

         kv   <- get(chains[k])
         vn   <- paste(chains[k],'Chain',sep='')
         kmat <- matrix(NA,ng,length(kv))

         gnames <- c(1:length(kv))

         if(!is.matrix(kv) & !is.null(names(kv)))gnames = names(kv)

         if(is.matrix(kv)){
            if( is.null(rownames(kv)) ) rownames(kv) <- c(1:nrow(kv))
            if( is.null(colnames(kv)) ) colnames(kv) <- c(1:ncol(kv))
            gnames <- outer(rownames(kv),colnames(kv),paste,sep='_')
         }
         colnames(kmat) <- gnames
         assign(vn,kmat)
         chainList <- append(chainList,list(get(vn)))
         names(chainList)[k] <- vn
     }

   }
         
   if(!is.null(sums)){
      
     ns <- length(sums)
     for(k in 1:ns){

        kv <- get(sums[k])
        v1 <- paste(sums[k],'Sum',sep='')
        v2 <- paste(sums[k],'Sum2',sep='')
        kg <- rep(0,length(kv))
        assign(v1,kg)
        assign(v2,kg)
        chainSum <- append(chainSum,list(get(v1)))
        names(chainSum)[length(chainSum)] <- v1
        chainSum <- append(chainSum,list(get(v2)))
        names(chainSum)[length(chainSum)] <- v2
     }
   }
   list(chainList = chainList, chainSum = chainSum)
}

loadChains <- function(g,burnin=1,chains=chains,sums=NULL, pars, sumVals,
                     chainList, chainSum){


  for(k in 1:length(chains)){
     chainList[[k]][g,] <- pars[[ which(names(pars) == chains[k]) ]]
  }

   if(g >  burnin & !is.null(sums)){
     for(k in 1:length(sums)){

       w1 <- which(names(chainSum) == paste(sums[k],'Sum',sep='') )
       w2 <- which(names(chainSum) == paste(sums[k],'Sum2',sep='') )

       ws <- which(names(sumVals) == sums[k])

       chainSum[[w1]] <- chainSum[[w1]] + sumVals[[ws]]
       chainSum[[w2]] <- chainSum[[w2]] + sumVals[[ws]]^2
     }
   }
   list(chainList = chainList, chainSum = chainSum)
}

processMCMC <- function( chains=NULL, sums = NULL, burnin=1, ng=nrow(chainList[[1]]), PLOTS=F,
                         outfolder=character(0),outfile='mcmc',chainList = chainList, chainSum = chainSum,
                         chainlength= 50000){

  CPLOT <- T
  file  <- outfile

  if(length(outfolder) > 0){
    if(!file.exists(outfolder))dir.create(outfolder)
    file <- paste(outfolder,'/',outfile,sep='')
  }

  postSummary <- postMuSe <- chainOut <- parMeans <- parSds <- numeric(0)

  kk <- 0

  if(!is.null(chains)){

    for(k in 1:length(chains)){

      if(PLOTS)CPLOT <- T
      x  <- col2Mat( chainList[[k]][burnin:ng,] )
      vk <- apply(x,2,var)
      kx <- which(vk > 0)
      if(length(kx) == 0)next

      kk <- kk + 1
 #     x  <- col2Mat(x[,kx])

      nx <- ncol(x)

      if(nx > 12)CPLOT <- F

      tmp <- processPars(x,rep(0,nx),CPLOT=CPLOT,DPLOT=F,
                  burnin=burnin)$summary[,c('estimate','se','0.025','0.975')]

      tmp <- row2Mat(tmp)
      rownames(tmp) <- paste(chains[k],rownames(tmp),sep='_')
      postSummary <- rowBind(postSummary,tmp,rownames(tmp))

      parMeans <- append(parMeans,list(tmp[,1]))
      parSds   <- append(parSds,list(tmp[,2]))
      names(parMeans)[kk] <- names(parSds)[kk] <- chains[k]

      if(nx > 12){

          i1 <- 1
          i2 <- 12
          while(i1 < ncol(x)){
             ii <- i1:i2
             processPars(x[,ii],rep(0,length(ii)),CPLOT=PLOTS)
             if(PLOTS)dev.copy2pdf(file=paste(file,i1,'_',i2,chains[k],'.pdf',sep=''))
             i1 <- i1 + 12
             i2 <- i2 + 12
             if(i2 > ncol(x))i2 <- ncol(x)
          }
       }

      if(nrow(x) > chainlength){
        thin <- round( seq(1,nrow(x),length=chainlength), 0)
        x <- as.matrix(x[thin,])
      }
      chainOut <- append(chainOut,list(x))

      if(PLOTS & CPLOT){
        dev.copy2pdf(file=paste(file,chains[k],'.pdf',sep=''))
      }
    }
  }

  if(!is.null(sums)){

    ntot <- ng - burnin

    for(k in 1:length(sums)){

       w1 <- which(names(chainSum) == paste(sums[k],'Sum',sep='') )
       w2 <- which(names(chainSum) == paste(sums[k],'Sum2',sep='') )

       mu <- chainSum[[w1]]/ntot
       se <- sqrt(chainSum[[w2]]/ntot - mu^2)

        postMuSe <- append(postMuSe,list(cbind(mu,se)))
     }
    names(postMuSe) <- sums
  }

  list(postSummary = postSummary, postMuSe = postMuSe, parMeans = parMeans, 
       parSds = parSds, chainOut = chainOut)
}


priorInteractionIndex <- function(vnames, iSymbol=':'){  #has all variable names, including interactions
  
  if('years' %in% vnames)vnames <- vnames[vnames != 'years']
  vn  <- vnames
  ivn <- grep(iSymbol,vnames)
  
  m <- numeric(0)
  
  if(length(ivn) == 0)return( list(isIntPrior = m, noIntPrior = m, intPrior = m) )
  
  aprior <- allPrior[allPrior %in% vn]
  
  vn <- vn[-ivn]
  ivn <- vnames[ivn]
  
  isIntPrior <- matrix(0,length(aprior),length(ivn))    #indicator: prior on main effects in interaction
  colnames(isIntPrior) <- ivn
  rownames(isIntPrior) <- aprior
  
  tmp <- matrix( unlist(strsplit(ivn,iSymbol)), nrow=2)
  
  for(m in 1:ncol(tmp)){
    mm <- match(tmp[,m],aprior)
    wm <- which(is.finite(mm))
    isIntPrior[mm[wm],m] <- 1
  }
  
  cc <- length(vn) + which(colSums(isIntPrior) == 0)    # neither has prior
  
  noIntPrior <- c(1:length(vn),cc)                    # main effects, interactions without priors
  intPrior   <- c(1:length(vnames))                      # interactions with priors
  if(length(noIntPrior) > 0)intPrior   <- intPrior[-noIntPrior]
  
  list(isIntPrior = isIntPrior, noIntPrior = noIntPrior, intPrior = intPrior)
}



updateIntBeta <- function(xx,yy,bb,lo,hi,sigma,isIntPrior,intPrior,noIntPrior,iSymbol='X'){  
  
  #allPrior in main program
  
  if(is.null(lo))lo <- bb*0 - 100
  if(is.null(hi))hi <- bb*0 + 100

  if(!is.matrix(lo))lo <- as.matrix(lo)
  if(!is.matrix(hi))hi <- as.matrix(hi)
  
  wi <- grep(iSymbol,rownames(bb))
  
  aprior <- rownames(isIntPrior)
  
  #main effects conditional on interaction terms
  
  yz <- yy - as.matrix(xx[,intPrior])%*%bb[intPrior,]   #from interactions
  
  bb[noIntPrior,]  <- bUpdateMVN_Rcpp(xx[,noIntPrior],yz,b=bb[noIntPrior,],
                                      lo=lo[noIntPrior,],hi=hi[noIntPrior,],sigma)
  
  nk <- ncol(bb)
  minB <- -bb[aprior,]      #min bound
  if(!is.matrix(minB))minB <- matrix(minB,nrow=length(aprior),ncol=nk)
  
  ss <- sample(length(wi))
  for(k in ss){

    wk <- wi[k]
    yz <- yy - xx[,-wk]%*%bb[-wk,]
    
    mink <- minB*isIntPrior[,k]
    mink[isIntPrior[,k] == 0] <- -500
    mm   <- apply(mink,2,which.max)
                    
    lok <- mink[cbind(mm,1:nk)]
    
    bb[wk,] <- bUpdateMVN_Rcpp(xx[,wk],yz,b=bb[wk,],lo=lok,hi=hi[wk,],sigma)
    
    minB <- minB - matrix(bb[wk,],nrow(minB),nk,byrow=T)*isIntPrior[,k]
  }
  bb
}

updateBetaMVN <- function(xx,yy,bb,lo=NULL,hi=NULL,sigma,delta,priorO=diag(1,nrow(sigma)),
                          priorDmu=rep(.1,nrow(sigma)),priorDvar=rep(1,nrow(sigma)),INT=F,
                          isIntPrior=NULL,intPrior=NULL,
                          noIntPrior=NULL,varLo=rep(1e-8,ncol(yy)),varHi=rep(10000,ncol(yy)),
                          SIGMA=T){
  
  #bb must have rownames

  ss <- dd <- sigma*0
  
  if(INT)bb <- updateIntBeta(xx,yy,bb=bb,lo=lo,hi=hi,sigma,
                               isIntPrior=isIntPrior,intPrior=intPrior,
                               noIntPrior=noIntPrior)
  
  if(!INT)bb  <- bUpdateMVN_Rcpp(xx,yy,b=bb,lo=lo,hi=hi,sigma)
  
  if(SIGMA){
 #   tmp   <- updateScInvWish(yy,xx%*%bb,delta=delta,
 #                            priorDmu=priorDmu,priorDvar=priorDvar,
 #                            varLo=rep(1e-8,ncol(yy)),varHi=rep(10000,ncol(yy)) )
    
    tmp   <- updateScInvWish(yy,xx%*%bb,delta=delta,
                             priorDmu=priorDmu,priorDvar=priorDvar,
                             varLo=varLo,varHi=varHi )
    
    ss <- tmp$sigma
    dd <- tmp$delta
  }
  
  list(beta = bb, sigma = ss, delta = dd)
  
} 






qprob <- function(p){  #evaluate likelihood for pathogen model: multinomial

    theta <- p[1]
    phi   <- p[2]
    s0    <- p[3]
    s1    <- p[4]
    q00 <- (1 - theta)*(1 - s0) + theta*(1 - phi)*(1 - s1)
    q01 <- (1 - theta)*s0 + theta*(1 - phi)*s1
    q10 <- theta*phi*(1 - s1)
    q11 <- theta*phi*s1

    c(q00,q01,q10,q11)
}

like_multinom <- function(pars){
   q <- qprob(pars)
  -sum(nvec*log(q))
}

outDataFrame <- function(...){  
  # a list of variable names coerced to dataframe (e.g. outDataFrame('a','b'))

  x  <- list(...)

  nx <- length(x)
  nc <- 0
  nr <- 0
  xn <- character(0)
  

  for(j in 1:nx){
    xj <- get(x[[j]])
    if(!is.null(dim(xj)) )xj <- as.matrix(xj)
    if( is.null(dim(xj)) )xj <- matrix(xj,ncol=1)
    nj <- nrow(xj)
    if(is.null(rownames(xj))){
      if(ncol(xj) == 1)colnames(xj) <- x[[j]]
      if(ncol(xj) > 1 )colnames(xj) <- paste(x[[j]],c(1:nj),sep='-')
   }
   xj <- t(xj)

    nr <- max(nr,ncol(xj))
    nc <- nc + nrow(xj)
    xn <- c(xn,rownames(xj))
  }

  out <- data.frame(matrix(NA,nr,nc))
  colnames(out) <- xn

  j1 <- 1

  for(j in 1:nx){
     xj <- get(x[[j]])
     if(is.matrix(xj)) xj <- t(xj)
     if(!is.matrix(xj))xj <- as.matrix(xj,ncol=1)
     j2 <- j1 + ncol(xj) - 1
     out[1:nrow(xj),j1:j2] <- xj
     j1 <- j2 + 1
  }
  format(out,justify='right')
}
    
write2file <- function(x,file,append=F,row.names=T){  
  
  #x is the name of an object (in quotes) or the object itself
  
  z <- x
  if(is.character(x))z <- get(x)
  rnames <- F
  
  if(append){
    xx <- as.data.frame( matrix(c(' ',x),2,1) )
    write.table(xx,file,row.names=row.names,col.names=F,append=append,quote=F)
  }

 # if(is.null(rownames(z)))rnames <- F

  if(!is.null(rownames(z)) & row.names){
    rnames <- rownames(z)
    rownames(z) <- NULL
    nr <- nchar(rnames)
    nc <- max(nr)
    spaces <- rnames
    nl <- length(spaces)
    while(nl > 0){
      spaces[nchar(spaces) < nc] <- paste(spaces[nchar(spaces) < nc],' ',sep='')
      nl <- length( which(nchar(spaces) < nc) )
    }
    rnames <- spaces
  }
   

 # if(!is.data.frame(z))z <- data.frame(z)

  z <- format(z,justify='right')
  write.table(z,file,append=append,quote=F,row.names=rnames)
}

updateSSRW <- function(states,y,zb=y*0,missing,tg,sg){        
  #state-space random walk 
  #update continuous states, random walk
  #missing times, obs y, obs error tg, process error sg
  #zb = z%*%beta
  
  for(t in 1:nt){
    
    VI <- 0
    v  <- 0
    
    if(!t %in% missing){          #observations
      v  <- y[t]/tg
      VI <- 1/tg
    }
     
    if(t < nt){              #t+1 term excluded for last 
      v  <- v + (states[t+1] - zb[t])/sg 
      VI <- VI + 1/sg
    }
    
    if(t > 1){                #t-1 term excluded for 1st 
      v  <- v + (states[t-1] + zb[t-1])/sg
      VI <- VI + 1/sg
    }
    
    V     <- 1/VI
    states[t] <- rnorm(1,V*v,sqrt(V))
  }
  states
}



initialStatesSS <- function(y){  # initialize states in SS model
  
  if(!is.matrix(y))y <- matrix(y,ncol=1)
  r <- ncol(y)
  
  n    <- length(y)
  time <- c(1:n)
  wm   <- which(is.na(y))
  notMiss <- c(1:n)
  
  x <- y
  
  if(length(wm) > 0){
    notMiss <- notMiss[-wm]
    x[wm]   <- predict(smooth.spline(time[-wm],y[-wm]),time)$y[wm]
  }
  
  list(x = as.vector(x), notMiss = notMiss, miss = wm)
}


updateSSbeta <- function(states,z,sg,priorB,priorIVB,addStates=F){  
  #covariates z
  #if addStates = T:  states[t+1] <- states[t] + zbeta[t]
  
  nt <- length(states)
  xx <- states[-1]
  zz <- z
  if(nrow(zz) == nt)zz <- zz[-nt,]     # last value already clipped?
  
  V  <- solve( crossprod(z[-nt,])/sg + priorIVB )
  
  if(addStates)xx <- xx - states[-nt]

  v  <- 1/sg*crossprod(zz,xx)
  t( rmvnorm(1,V%*%v,V) )
}

updateSSB <- function(states,sg,priorB,priorIVB){  # SS model
	
  nt <- length(states)
  
	V <- 1/( (nt - 1)/sg + priorIVB)
	v <- (states[nt] - states[1])/sg + priorIVB%*%priorB
	rnorm(1,V*v,sqrt(V))
}


updateSSparsMet <- function(b,priorB,priorVB,
                            lo=rep(-Inf,length(b)),hi=rep(Inf,length(b)),
                            x,sigma,dt=1,propVar){
  
  k <- length(b)
  if(k == 1)pb <- tnorm(1,lo,hi,b,propVar)
  if(k > 1) pb <- tnorm.mvt(b,b,propVar,lo,hi) 
  
  xnow <- xg[-nt] + fx(xg[-nt],b)*dt     #predicted x
  xnew <- xg[-nt] + fx(xg[-nt],pb)*dt
	
  pnow <- sum(dnorm(xg[-1],xnow,sqrt(sigma*dt),log=T)) +
          dmvnorm(t(b),priorB,priorVB,log=T)
  pnew <- sum(dnorm(xg[-1],xnew,sqrt(sigma*dt),log=T)) +
	      dmvnorm(t(pb),priorB,priorVB,log=T)
  atmp <- acceptMH(pnow,pnew,b,pb,BLOCK=T)
  b   <- atmp$x

  list(b = b, aa = atmp$accept)

}
updateSSstatesMet <- function(x,y,b,sigma,tau,lo=-Inf,hi=Inf,notMiss=c(1:length(x)),
                              dt=1,propSd=.1){   #propose/accept x as a block
  nt   <- length(x)
  xp   <- tnorm(nt,lo,hi,x,rexp(nt,1/propSd) )       #proposed x
  
  xnew <- xp[-nt] + fx(xp[-nt],b)*dt      #mean proposed
  xnow <- x[-nt]  + fx(x[-nt],b)*dt       #mean current

  pnow <- sum(dnorm(x[-1],xnow,sqrt(sigma*dt),log=T)) +
	       sum(dnorm(y[notMiss],x[notMiss],sqrt(tau),log=T))
	       
  pnew <- sum(dnorm(xp[-1],xnew,sqrt(sigma*dt),log=T)) +
	       sum(dnorm(y[notMiss],xp[notMiss],sqrt(tau),log=T))
	       
  tmp <- acceptMH(pnow,pnew,x[-1],xp[-1],BLOCK=T)
  x[-1]   <- tmp$x

  list(x = x, aa = tmp$accept)
}

diagXcovar <- function(diagmat,covmat){

  #  diagmat %*% sigma %*% t(diagmat) 
  #  diagmat - rows are diagonal of a k by k diagonal matrix
  #  covmat  - covariance matrix (symm, pos-def)
  #  returns nn by k^2 matrix, each row are elements of a k by k matrix
  
  nn <- nrow(diagmat)
  k  <- nrow(covmat)
  lm <- matrix(0,nn,k^2)

  km <- matrix(1:(k^2),k,k)

  ii <- km[!upper.tri(km)]
  jj <- t(km)[!upper.tri(km)]

  kk <- 0
  for(i in 1:k){
    for(j in i:k){
      kk <- kk + 1
      lm[,ii[kk]] <- lm[,jj[kk]] <- diagmat[,i]*diagmat[,j]*covmat[i,j]
    }
  }
  lm
}

histWeight <- function(xx,wt,bins,removeZeros=T,useMids=F){

  # xx   - data to be binned
  # wt   - wt assigned to each x
  # bins - partition for x

  #  example: wt is 1/(plotarea), then results are no/area

  xx <- as.vector(xx)
  wt <- as.vector(wt)

  nb <- length(bins)
  db   <- diff(bins)

  if(useMids){

    mids <- bins[-nb] + db/2
    bins <- c(mids,bins[nb]+db[nb-1]/2)
  }
  
  if(!removeZeros)xx[xx < bins[1]] <- bins[1]
  
  x1 <- findInterval(xx,bins) 
  if(removeZeros)x1[xx == 0] <- NA


  freq <- by(as.vector(wt),as.vector(x1),sum,na.rm=T) 

  wb <- as.numeric( names(freq) )
  whist <- bins*0
  whist[wb] <- unlist(freq)
  
  db <- c(db,db[length(db)])
  dens <- whist/sum(whist)

  list(wtHist = whist, bins = bins, dens=dens)
}
  

histByYr <- function(xx,yr,xbin=1,breaks=NULL,dens=F,dt=1,weight=xx*0 + 1,deathyr=NULL){

  # xx has observations as rows, years as columns
  # yr is the vector of years corresponding to columns
  # change - if true, return change in xx distribution
  # dens   - if true, return density
  # deathyr    - length nn vector of death years

  ny <- length(yr)
  nn <- nrow(xx)
  
  if(is.null(deathyr))deathyr <- rep(NA,nn)
  
  dt <- 1

  maxx   <- round( max(xx,na.rm=T) ,0) + xbin
  if(is.null(breaks))breaks <- seq(0,maxx,by=xbin)
  mx     <- length(breaks) 

  times  <- c(1:round(ny/dt,0) )
  nt     <- length(times)
  
  gMat <- matrix(0,mx,nt)

  nh <- length(breaks)
  
  xx[xx == 0] <- NA
  
  ymat <- matrix(yr,nrow(xx),length(yr),byrow=T)
  incr <- t( apply(xx,1,diff)/apply(ymat,1,diff) )
  if(nrow(incr) == 1 & nrow(xx) > 1)incr <- t(incr)
  
  liveMat <- ymat*(xx*0 + 1)
  wlive <- which(is.finite(incr),arr.ind=T)
  wlNo  <- which(is.finite(xx),arr.ind=T)
  
  incr[incr < 0] <- 0
  
  inow  <- findInterval(xx[wlive],breaks)
  iyr   <- match(ymat[wlive],yr)
  
  inowNo  <- findInterval(xx[wlNo],breaks)
  iyrNo   <- match(ymat[wlNo],yr)
  
  groMu <- round(byRcpp(incr[wlive],inow,iyr,gMat*0,gMat*0, fun='mean'),3)
  groVr <- signif(byFunction(incr[wlive],inow,iyr,gMat*0,var),3)
  groNo <- byRcpp(xx[wlNo]*0+1,inowNo,iyrNo,gMat*0,gMat*0, fun='sum')
  
  wt <- byRcpp(weight[wlNo],inowNo,iyrNo,gMat*0,gMat*0, fun='mean')
  
  dist  <- groNo*wt
  ddist <- groNo[,-1] - groNo[,-nt]
  
  deadTable <- groNo*0
  
  wd        <- which(is.finite(deathyr))
  
  if(length(wd) > 0){
    
    wdead     <- cbind( wd, match(deathyr[wd],yr) - 1 )
    
    diamd     <- xx[wdead]
    wdd       <- which(!is.finite(diamd))
    if(length(wdd) > 0){
      dd <- xx[wdead[wdd,1],]
      if(!is.matrix(dd))dd <- matrix(dd,nrow=1)
      dd <- apply(dd,1,max,na.rm=T)

      diamd[wdd] <- dd
      missyr       <- findInterval( deathyr[wd[wdd]] ,yr) + 1 
      missyr[missyr > nt] <- nt
      missyr[missyr == 0] <- 1
      wdead[wdd,2] <- missyr
      wi <- which(!is.finite(diamd))
      if(length(wi) > 0){
        diamd <- diamd[-wi]
        wdead <- wdead[-wi,]
        wd    <- wd[-wi]
        if(!is.matrix(wdead))wdead <- matrix(wdead,1)
      }
    }
    
    idead     <- findInterval(diamd,breaks)
    
    if(length(wd) == 1){
      cc <- cbind(idead,wdead[2])
      deadTable[cc] <- 1
    }
    if(length(wd) > 1)deadTable <- byRcpp(rep(1,length(wd)),idead,wdead[,2],gMat*0,gMat*0, fun='sum')
  }

  dsurv     <- groNo - deadTable
  dsurv[dsurv < 0] <- 0              #check this
  
  if(!is.matrix(ddist))ddist <- matrix(ddist,ncol=1)


  groNo[is.na(groNo)] <- 0
  w0 <- which(groNo == 0,arr.ind=T)
  groMu[w0] <- groVr[w0] <- NA
  
  dsurv[,nt] <- dist[,nt]
  
  colnames(dist) <- yr
  colnames(ddist) <- yr[-ny]

  list(dist = dist, change = ddist, size = breaks+.5*dt,
       surv = dsurv, groMu  = groMu, groVr = groVr, groNo = groNo)
}


getMonthNames  <- function(){

  mm <- c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')
  invisible(mm)
}


subETOPO5 <- function(lon, lat){

  #lon is negative degrees W, positive degrees E
  #lat is positive N
  
  require( GEOmap )
  require( geomapdata )
  
  data(ETOPO5)
  
  mapx <- range(lon)
  mapy <- range(lat)
  
  
  zz <- getETOPO(ETOPO5, glat = mapy, glon = mapx ) #latitudes are upside down
  ilon <- seq(mapx[1], mapx[2], length = nrow(zz))
  ilat <- seq(mapy[1], mapy[2], length = ncol(zz))
  
 # lgrid <-as.matrix( expand.grid(ilon,ilat) )
  
 # ww <- which(zz > 0, arr.ind=T) 
 # plot(ilon[ww[,1]],ilat[ww[,2]], cex=.1)
  
  list(x = ilon, y = ilat, z = zz )
  
}
  
  

subETOPO5old <- function(lon, lat){
  
  #lon is negative degrees W, positive degrees E
  #lat is positive N
  
  require(geomapdata)
    
    data(ETOPO5)

  ilon <- nrow(ETOPO5)/360
  ilat <- ncol(ETOPO5)/180
  
  x <- lon
  y <- lat
  
  
  xx <- 360 + lon
  
#  lon[lon <= 0] <- -lon[lon <= 0]  + 180  #(0, 360) scale
  lat[lat > 0] <- 90 - lat[lat > 0]
  lat[lat < 0] <- -lat[lat < 0] + 90

#  rlon <- lon*ilon
  rlon <- xx*ilon
  rlat <- lat*ilat

  z <- ETOPO5[rlon[1]:rlon[2],]
  z <- z[,rlat[1]:rlat[2]]

  x <- seq(x[1],x[2],length=nrow(z))
  y <- seq(y[1],y[2],length=ncol(z))

  list(x = x, y = y, z = z)

}


long2UTMZone <- function(lon){     #assumes west of prime meridian, (-180, 0)
  floor((lon + 180)/6 %% 60) + 1
}




UTM2latlon <- function( xy, zone = NULL, lon = NULL, southern = FALSE ){
  
  #xy has 2 columns or 3 columns, x, y, zone
  
  if(ncol(xy) == 2){
    if( is.null(zone) & is.null(lon) )stop('zone or lon must be specified')
    if( is.null(zone) )zone <- long2UTMZone(lon)
    xy <- cbind(xy, zone)
  }
  
  if(length(southern) == 1)southern <- rep(southern, nrow(xy))
  
  require(rgdal)
  
  colnames(xy) <- c("x","y","zone")
  
  zones <- sort(unique(xy[,'zone']))
  
  ll <- xy[drop=F,,1:2]
  
  for(j in zones){
    
    wj  <- which(xy[,'zone'] == j)
    xyj <- as.data.frame( xy[drop=F, wj,] )
    
    ws <- which( southern[wj] )
    wn <- which( !southern[wj] )
    
    if(length(ws) > 0){
      xys <- as.data.frame( xyj[drop=F,ws,] )
      coordinates(xys) <- c("x", "y")
      cstring <- paste('+proj=utm +zone=', j,' +ellps=WGS84 +south', sep='')
      upoints <- suppressWarnings( SpatialPoints(coords=xys,
                                                 proj4string=CRS(cstring)) )
      ut <- suppressWarnings( spTransform(upoints, 
                                          CRS("+proj=longlat +ellps=WGS84 +south")) ) # ll to utm
      ut <- as.data.frame(ut)
      if(nrow(ut) == 1){
        ut <- unlist(ut)
      }else{
        ut <- as.matrix(ut)
      }
      ll[ wj[ws],1:2] <- ut
    }
    if(length(wn) > 0){
      xys <- as.data.frame( xyj[drop=F,wn,] )
      coordinates(xys) <- c("x", "y")
      
      cstring  <- paste('+proj=utm +zone=', j,' +ellps=WGS84 +north', sep='')
      upoints <- suppressWarnings( SpatialPoints(coords=xys,
                                                 proj4string=CRS(cstring)) )
      ut <- suppressWarnings( spTransform(upoints, 
                                          CRS("+proj=longlat +ellps=WGS84 +north")) )
      ut <- as.data.frame(ut)
      if(nrow(ut) == 1){
        ut <- unlist(ut)
      }else{
        ut <- as.matrix(ut)
      }
      ll[ wj[wn],1:2] <- ut
    }
  }
  colnames(ll) <- c('lon','lat')
  ll
}

lonLat2UTM <- function( lonLat, zone=NULL ){
  
  #xy has 2 columns, long, lat
  #long is longitude, used to find UTM zone
  
 # require(PBSmapping)
  require(rgdal)
  
  if( length(lonLat) == 2 & !is.matrix(lonLat) & !is.data.frame(lonLat) )
    lonLat <- matrix(lonLat, 1 )
  
  lonLat <- as.data.frame(lonLat)
  
  colnames(lonLat) <- c('x','y')
  
  if(is.null(rownames(lonLat)))rownames(lonLat) <- c(1:nrow(lonLat))
  
  wf <- which(is.finite(lonLat[,1]) & is.finite(lonLat[,2]))
  xy <- lonLat[drop=F, wf,]
  
  if(length(wf) < nrow(lonLat)){
    wna <- which( !c(1:nrow(xy)) %in% wf )
    warning( paste('missing values in lonLat:', paste0(wna, collapse=', ')) )
  }
  
  utm <- xy*0
  utm$zone <- NA
  
  zone <- long2UTMZone(lonLat[wf,'x'])
  
  zz <- sort(unique(zone))
  nz <- length(zz)
  
  for(j in zz){
    
    wj  <- which(zone == j)
    xyj <- as.data.frame( xy[drop=F, wj,] )
    utm$zone[wj] <- j
    
    ws <- which(xyj[,2] < 0)
    wn <- which(xyj[,2] >= 0)
    
    if(length(ws) > 0){
      xys <- as.data.frame( xyj[drop=F,ws,] )
      coordinates(xys) <- c("x", "y")
      cstring <- paste('+proj=utm +zone=', j,' +ellps=WGS84 +south', sep='')
      
      llpoints <- suppressWarnings( SpatialPoints(coords=xys,
                                proj4string=CRS("+proj=longlat +ellps=WGS84 +south")) )
      ut <- suppressWarnings( spTransform(llpoints, CRS(cstring)) ) # ll to utm
      ut <- as.data.frame(ut)
      utm[wj[ws],1:2] <- ut
    }
    if(length(wn) > 0){
      xys <- as.data.frame( xyj[drop=F,wn,] )
      coordinates(xys) <- c("x", "y")
      cstring <- paste('+proj=utm +zone=', j,' +ellps=WGS84 +north', sep='')
      
      llpoints <- suppressWarnings( SpatialPoints(coords=xys,
                                proj4string=CRS("+proj=longlat +ellps=WGS84 +north")) )
      ut <- suppressWarnings( spTransform(llpoints, CRS(cstring)) ) # ll to utm
      ut <- as.data.frame(ut)
      utm[wj[wn],1:2] <- ut
    }
  }

  colnames(utm) <- c('UTMx','UTMy','UTMzone')
  unew <- lonLat*0 + NA
  unew$UTMzone <- NA
  unew[rownames(utm),] <- utm
  
  unew
}

UTM2latlonOld <- function(xy, zone = NULL, long=NULL, southern=F){
  
  require(PBSmapping)
  
  #xy has 2 columns, 'X', 'Y' for UTMx, UTMy
  #long is longitude, used to find UTM zone
  
  if(is.null(zone) & is.null(long))stop('reguire either zone or long')
  
  xy <- as.matrix(xy)
  colnames(xy) <- c('X','Y')
  
  if(is.null(zone))zone <- long2UTMZone(long)
  
  if(length(zone) == 1)zone <- rep(zone,nrow(xy))
  
  zi <- unique(zone)
  utm <- xy*0
  
  for(i in zi){
    wi  <- which(zone == i & is.finite(rowSums(xy)))
    xyi <-  matrix( xy[wi,], ncol=2 )
    colnames(xyi) <- c('X','Y')
    attr(xyi,'zone') <- i
    attr(xyi,'projection') <- 'UTM'
    utm[wi,] <- as.matrix( convUL(xyi,km=F,southern=southern) )
  }
  utm
}


columnSplit <- function(vec, sep='_', ASFACTOR = F, ASNUMERIC=F,
                        LASTONLY=F){
  
  vec <- as.character(vec)
  nc  <- length( strsplit(vec[1], sep, fixed=T)[[1]] )
  
  mat <- matrix( unlist( strsplit(vec, sep) ), ncol=nc, byrow=T )
  if(LASTONLY & ncol(mat) > 2){
    rnn <- mat[,1]
    for(k in 2:(ncol(mat)-1)){
      rnn <- columnPaste(rnn,mat[,k])
    }
    mat <- cbind(rnn,mat[,ncol(mat)])
  }
  if(ASNUMERIC){
    mat <- matrix( as.numeric(mat), ncol=nc )
  }
  if(ASFACTOR){
    mat <- data.frame(mat)
  }
  if(LASTONLY)mat <- mat[,2]
  mat
}

columnPaste <- function(c1, c2, sep='-', NOSPACE = FALSE){
  
  c1    <- as.character(c1)
  c2    <- as.character(c2)
  if(NOSPACE){
    c1   <- .replaceString(c1, ' ', '')
    c2   <- .replaceString(c2, ' ', '')
  }
  c12   <- apply( cbind(c1, c2) , 1, paste0, collapse=sep)
  
  c12
}

distributionMapLittle <- function(spec, 
                                  dataLoc="/Volumes/research/clark/clark.unix/fia/little/" ){
  #returns shape file
  
  require(maptools)
  require(foreign)

  treeCodes <- read.csv('/Volumes/research/clark/clark.unix/traitTable.csv',
                          stringsAsFactors = F)
  fiaCode <- as.character(treeCodes[ match(spec,treeCodes[,'code4']),'codeUSDA'] )
  
  fold <- tolower(spec)
  ff <- list.files( paste(dataLoc,fiaCode,sep=''), full.names=T )
#  ff <- paste(dataLoc,fiaCode,ff,sep='/')
  
 # lf <- list.files(ff)
  
  file <- ff[grep('dbf',ff)]

  frame <- read.dbf( file )
 # attr  <- read.dta( ff[grep('dta',ff)] )
   shap  <- readShapePoly( ff[grep('shp',ff)] )
  shap
}

dayLength <- function(JD,lat){
  
  nl  <- length(lat)
  nj  <- length(JD)
  latNew <- lat/180
  
  JDmat <- JD
  
  if(nl > 1 & nj > 1){
    latNew <- matrix(latNew,nl,nj)
    JDmat  <- matrix(JD,nl,nj,byrow=T)
  }
  
  p <- asin( .39795*cos(.2163108 + 2*atan(.9671396*tan(.0086*(JDmat - 186)))) )
  tmp <- 24 - (24/pi)*acos( (sin(.8333*pi/180) + sin(latNew*pi)*sin(p)) /
                       (cos(latNew*pi)*cos(p)))
  tmp
}
monthlyPET <- function(yi,mi,tempMatrix,precMatrix,lat){
  
  # m  = nyr*12 (serial months)
  # mi = month index (1:12)
  # yi = yr index    (e.g., 1990, 1991)
  # returns n X month*nyr matrix and n X nyr matrix 
  # tempMatrix - n X m location by month
  # precMatrix - n X m location by month
  # yi     - length-m year index
  # mi     - length-m month index
  
  if(!is.matrix(tempMatrix))tempMatrix <- matrix(tempMatrix,1)
  if(!is.matrix(precMatrix))precMatrix <- matrix(precMatrix,1)
  
  n <- nrow(tempMatrix)
  m <- ncol(tempMatrix)
  yrvec <- min(yi):max(yi)
  yseq  <- yrvec - yrvec[1] + 1
  nyr   <- length(yrvec)
  ymat  <- matrix(yseq,n,nyr,byrow=T)
  
  yindex <- match(yi,yrvec)
  
  JD <- c(1:365)
  
  dl <- dayLength(JD,lat)
  mo <- daysSinceDate(1,1,yrvec[1],c(2:12),1,yrvec[1],nineteen=F)
  di <- findInterval(JD,mo) + 1
  
  yvec <- as.vector( matrix(yindex,n,m,byrow=T) )
  ivec <- as.vector( matrix(c(1:n),n,m) )
  
  # heat index by year
  tmp1   <- tempMatrix
  tmp1[tmp1 < 0] <- 0
 # HI    <- byRcpp(as.vector((tmp1/5)^1.514),ivec,yvec,ymat*0,ymat*0, fun='sum')
  
  
  index <- list(year = yvec, month = ivec)
  HI    <- tapply(as.vector((tmp1/5)^1.514), index, sum)
  
  
  HI    <- HI[yvec]
  alpha <- (6.75*1e-7)*HI^3 - (7.71*1e-5)*HI^2 + (1.792*1e-2)*HI + .49239
  
  tmp <- aggregateSequence(di,dl,action='mean')  #mean monthly daylength
  dlMon <- tmp[[1]]
  dlLen <- tmp[[2]]
  dlLen <- dlLen[,mi]
  
  PET     <- 16*dlLen/12*(10*tmp1/HI)^alpha
  list(PET = PET, daylength = dlLen)
}

monthlyPHr <- function(yi,mi,tempMatrix=tempj,precMatrix=precj,lat=lonlatj[2],tempThresh=4){
  
  tmp <- monthlyPET(yi,mi,tempMatrix,precMatrix,lat)
  PET <- tmp$PET
  daylen <- tmp$daylength
  
  pExcess <- precMatrix - PET
  
  dayFraction <- daylen
  #  pHours  <- dayFraction*pExcess
  #  pHours[pHours < 0] <- 0
  
  posE <- pExcess        #surplus
  posE[posE > 0] <- 1
  posE[posE <= 0] <- 0
  posD <- 1 - posE       #deficit
  
  DD    <- tempMatrix
  DD[DD < tempThresh] <- 0   #threshold for DD

  degreeHrs      <- DD*dayFraction
  degreeHrExcess <- degreeHrs*pExcess
  degreeHrPos    <- degreeHrs*posE        #degreeMonthHrs
  degreeHrNeg    <- degreeHrs*posD
  
  tmp <- aggregateSequence(yi,degreeHrExcess,action='sum')  #total annual daylength
  degHrExYr <- tmp$data
  
  tmp <- aggregateSequence(yi,degreeHrPos ,action='sum')  #total annual daylength
  degHrPosYr <- tmp$data
  
  tmp <- aggregateSequence(yi,degreeHrNeg ,action='sum')  #total annual daylength
  degHrNegYr <- tmp$data
  
  tmp <- aggregateSequence(mi,degreeHrExcess ,action='sum')  #total monthly daylength
  degHrExMo <- t(tmp$data)
  
  tmp <- aggregateSequence(mi, degreeHrPos ,action='sum')  #total monthly daylength
  degHrPosMo <- tmp$data
  
  list(degHrExYr = degHrExYr, degHrPosYr = degHrPosYr, degHrNegYr = degHrNegYr,
       degHrExMo = degHrExMo, degHrPosMo = degHrPosMo,
       degreeHrExcess = degreeHrExcess, degreeHrPos = degreeHrPos, degreeHrNeg = degreeHrNeg)
}

##############################################################################3
plotstart <- function(plotfile,REMOTE=F){
  
  # initiates plot
  # postscript plotfile ends in .ps
  # pdf plotfile ends in .pdf
  
  
  cc  <-NULL
  PDF <- F
  if(length(grep('.ps',plotfile)) > 0) cc <- 1
  if(length(grep('.pdf',plotfile)) > 0)cc <- 2
  
  if(is.null(cc)){
    warning('plotfile must end in .ps or .pdf')
    return()
  }
  if(cc == 2)PDF <- T
  
  if(!REMOTE)graphics.off()
  if(REMOTE & !PDF)postscript(file=plotfile,width=6, height=9,horizontal=FALSE)
}

####################################################

plotend <- function(plotfile,REMOTE=F){
  
  # terminates plot
  # postscript plotfile ends in .ps
  # pdf plotfile ends in .pdf
  
  cc  <- NULL
  PDF <- F
  if(length(grep('.ps',plotfile)) > 0) cc <- 1
  if(length(grep('.pdf',plotfile)) > 0)cc <- 2
  
  if(is.null(cc)){
    warning('plotfile must end in .ps or .pdf')
    return()
  }
  if(cc == 2)PDF <- T
  
  
  if (REMOTE)dev.off()
  if(!REMOTE)dev.print(device=postscript,file=plotfile,width=7, 
                       height=10, horizontal=FALSE)
  if(REMOTE & PDF)dev.copy2pdf(file=plotfile)

  
}

speciesCountSummary <- function(y,aggregate=F,TYPE='continuous'){
  
  #y is a sample by species matrix of counts
  
  S <- ncol(y)
  n <- nrow(y)
  
  ii <- list(spec = as.vector(matrix(c(1:S),n,S,byrow=T)), ss = as.vector(y))
  classBySpec <- imat <- byIndex(as.vector(y)*0+1,ii,sum)
  
  classAll <- table(y)
  
  lowestClass <- matrix(1:ncol(classBySpec),S,ncol(classBySpec),byrow=T)
  
  tmp <- classBySpec
  tmp[tmp > 0] <- 1
  tmp <- tmp*lowestClass  
  tmp[tmp == 0] <- NA
  lowestClass <- apply(tmp,1,min,na.rm=T)
  highestClass <- apply(tmp,1,max,na.rm=T)
  
  ii <- list(spec = as.vector(matrix(c(1:S),n,S,byrow=T)), ss = as.vector(y))
  imat <- byIndex(as.vector(y)*0+1,ii,sum)
  imat[imat > 1] <- 1
  isum <- rowSums(imat)
  novary <- which(isum < 2)

  if(length(novary) > 0 & !TYPE %in% c('continuous', 'composition') ){
    warning( paste('novariation in species: ',names(isum)[novary],
                                      sep='') )
  }
                                                                          
  #aggregate spp that don't vary
  
  if(aggregate){
    novary <- 100
    while(length(novary) > 0){
      ii <- list(spec = as.vector(matrix(c(1:S),n,S,byrow=T)), ss = as.vector(y))
      classBySpec <- imat <- byIndex(as.vector(y)*0+1,ii,sum)
      imat[imat > 1] <- 1
      isum <- rowSums(imat)
      novary <- which(isum < 2)
      if(length(novary) == 1){
        y <- y[,-novary]
      }
      if(length(novary) > 1){
        other <- rowSums(y[,novary])
        y <- cbind(y[,-novary],other)
      }
      S <- ncol(y)             
      snames <- colnames(y)
    }
  }
  
  list(y = y, classBySpec = classBySpec, classAll = classAll, 
       lowestClass = lowestClass, highestClass = highestClass)
}



multivarChainNames <- function(rowNames,colNames){
  as.vector( t(outer(colNames,rowNames,paste,sep='_')) )
}

add2matrix <- function(values,xmat=NULL){
  
  #xmat   - n X ? matrix with one row, columns are integer values
  #values - length-n vector be added/slotted in to xvec
  
  if(is.null(xmat)){
    n    <- length(values)
    cc   <- sort(unique(values))
    xmat <- matrix(0,n,length(cc),dimnames = list(1:n,cc))
    xmat[ cbind( c(1:n),match(values,cc)) ] <- 1
    return(xmat)
  }
  
  n <- nrow(xmat)
  if(length(values) != n)stop('vector length must equal rows in xmat')
  
  all <- sort( unique( c(values,as.numeric(colnames(xmat))) ))
  nc       <- length(all)
  
  xnew <- matrix(0,n,nc,dimnames = list(1:n,all))
  xnew[,colnames(xmat)] <- xmat
  
  xnew[ cbind(c(1:n),match(values,all)) ] <- xnew[ cbind(c(1:n),match(values,all)) ] + 1
  xnew
}



add2vector <- function(values,xvec=NULL){
  
  #xvec - matrix with one row, columns are integer values
  #values - to be added/slotted in to xvec
  
  if(is.null(xvec)){
    xvec <- table(values)
    cc   <- names(xvec)
    xvec <- matrix(xvec,1,dimnames = list(' ',cc))
    return(xvec)
  }
  x1 <- table(values)
  
  allNames <- sort( unique(as.numeric( c(names(x1),colnames(xvec)))) )
  nc       <- length(allNames)
  
  xnew <- matrix(allNames*0,1,dimnames = list(' ',allNames))
  ww   <- which(names(x1) %in% allNames)
  if(length(ww) > 0)xnew[,names(x1)[ww]] <- x1[ww]
  ww   <- which(colnames(xvec) %in% allNames)
  if(length(ww) > 0)xnew[,colnames(xvec)] <- xnew[,colnames(xvec)] + xvec[,ww]
  
  xnew
}

predictY2X_nonLinear <- function(xx,yy,bb,ss,priorIV = diag(1e-10,ncol(xx)), 
                                 priorX=matrix(0,ncol(xx)),
                                 predCols=c(2:ncol(xx)),isInt=NULL,intMat=NULL,
                                 isFactor=NULL,factorList=NULL){
  
  #inverse prediction for multivariate nonlinear in x and factors, metropolis
  
  iFcol  <- NULL
  priorX <- priorX[predCols]
  if(!is.matrix(priorX))priorX <- matrix(priorX)
  
  nn <- nrow(yy)
  intercept <- xx[,1]
  
  xnew <- xx
  xnew[,predCols] <- myrmvnorm(nn,xx[,predCols],diag(.01,length(predCols)))
  
  if(!is.null(isFactor)){          # all factors, main effects
    xnew[,isFactor] <- 0
    for(k in 1:length(factorList)){
      nf <- length(factorList[[k]]) + 1
      wk <- sample(1:nf,nn,replace=T)
      tm <- cbind( 1:nn, wk )
      xnew[,c('intercept',factorList[[k]])][tm] <- 1
    }
    iFcol <- match(isFactor,colnames(xx))
  }
  
  if(length(intMat) > 0){                     #are some of the nlin terms interactions?
    pindex <- unique(as.vector(intMat[,-1]))
    if(length(iFcol) > 0)pindex <- pindex[!pindex %in% iFcol]  # those that are not factors
    if(length(pindex) > 0) xnew[,pindex] <- myrmvnorm(nn,xx[,pindex],diag(.01,length(pindex)))
    xnew[,intMat[,1]] <- xnew[,intMat[,2]]*xnew[,intMat[,3]]
  }
  
  pnow <- dmvnormZeroMean(yy - xx%*%bb,smat=ss)
  pnew <- dmvnormZeroMean(yy - xnew%*%bb,smat=ss)
  
  a  <- exp(pnew - pnow)
  z  <- runif(nn,0,1)
  wa <- which(z < a)
  xx[wa,] <- xnew[wa,]
  
  list(x = xx, accept = length(wa))
}

predictY2X_linear <- function(xx,yy,bb,ss,priorIV = diag(1e-10,ncol(xx)), 
                              priorX=matrix(0,ncol(xx)),predCols=c(2:ncol(xx))){
  
  #inverse prediction for multivariate linear in x
  
  priorX <- priorX[predCols]
  if(!is.matrix(priorX))priorX <- matrix(priorX)
  
  nn <- nrow(yy)
  notPred <- c(1:ncol(xx))[-predCols]
  
  bn <- matrix(bb[notPred,],length(notPred))
  bp <- matrix(bb[predCols,],length(predCols))
  
  yi <- yy - xx[,notPred]%*%bn
  pp <- length(predCols)
  
  bs <- bp%*%solve(ss)
  
  V <- solve( bs%*%t(bp) + priorIV[predCols,predCols] )
  v <- yi%*%t(bs) + matrix( priorIV[predCols,predCols] %*% priorX,nn,pp,byrow=T)
  mu <- v%*%V
  if(ncol(mu) > 1) xn <- myrmvnorm(nn,mu,V)
  if(ncol(mu) == 1)xn <- rnorm(nn,mu,sqrt(V))
  xx[,predCols] <- xn
  xx
}
conditionalMVNRcpp <- function(xx, mu, sigma, cdex, p=ncol(mu)){  
  
  # xx,mu are matrices
  
  if(!length(xx) > nrow(sigma))return( conditionalMVN(xx,mu,sigma,cdex) )
  
  gdex <- (1:p)[-cdex] - 1
  cdex <- cdex - 1
  conditionalMVNInline_cpp(cdex, gdex, xx, mu, sigma)
}

tnorm.mvtVectorize <- function(column,w,mu,smat,lo,hi){
  
  n   <- nrow(w)
  
  tmp <- conditionalMVNRcpp(w,mu,smat,column)
  mu  <- tmp$mu
  vr  <- max(tmp$vr,1e-10)
  tmp <- tnorm(n,lo[,column],hi[,column],mu,sqrt(vr))
  tmp
}

tnorm.mvtVec <- Vectorize(FUN=tnorm.mvtVectorize, vectorize.args="column")


imputX_MVN <- function(xx,yy,beta,wmiss,sinv,xprior=0,xbound=NULL,priorWT=1){
  
  # priorWT is inverse of variance
  
  wcol <- unique(wmiss[,2])
  S    <- nrow(sinv)
  Q    <- nrow(beta)
  
  if(is.null(xbound))xbound <- apply(x,2,range,na.rm=T)
  
  for(j in wcol){
    
    wj <- wmiss[wmiss[,2] == j,]
    if(!is.matrix(wj))wj <- matrix(wj,1,2)
    wr <- wj[,1]
    xp <- xprior[wmiss[,2] == j]
    
    bj <- matrix(beta[j,],1)
    bn <- matrix(beta[-j,],Q - 1)
    
    xn <- matrix(xx[wr,-j],length(wr))
    z <- yy[wr,] - xn%*%bn
    datwt <- bj%*%sinv%*%t(bj)
    V     <- 1/( datwt + priorWT*datwt )
    v     <- z %*%sinv%*%t(bj) + xp*priorWT
    
    xx[wj] <- tnorm(length(wr),xbound[1,j],xbound[2,j],v%*%V,sqrt(V))
  }
  xx
}


shadeInterval <- function(xvalues,loHi,col='grey',PLOT=T,add=T,xlab=' ',ylab=' '){
  
  #draw shaded interval
  
  tmp <- smooth.na(xvalues,loHi)
  xvalues <- tmp[,1]
  loHi <- tmp[,-1]
  
  xbound <- c(xvalues,rev(xvalues))
  ybound <- c(loHi[,1],rev(loHi[,2]))
  if(!add)plot(xvalues,loHi[,1]*0,cex=.01,ylim=c(range(loHi,na.rm=T)),
               xlab=xlab,ylab=ylab)
  if(PLOT)polygon(xbound,ybound,border=NA,col=col)
  
  invisible(cbind(xbound,ybound))
  
}



getW2Y <- function(ww,pg){ 1 - (1 - pg)^(ww/pg) }
getY2W <- function(ww,pg){ pg*log(1 - ww)/log(1 - pg) }


Y2W <- function(yy,pg){    # y and w with same number of columns (expand y)
  
  S <- ncol(yy)
  wy <- yy
  W  <- 1 - rowSums(yy[,-S])
  
  wh <- which(W > pg)
  
  if(length(wh) > 0){
    expand  <- pg*log(1 - W[wh])/log(1 - pg)/W[wh]
    wy[wh,] <- wy[wh,]*expand
  }
  wy
}


oceanValues <- function(minz=0,minc=0,xx,yy,zz,lonLat,climVec){
  
  # adds low values to climVec so ocean is not shaded on map
  # minz - minimum elevation to define ocean
  # minc - minimum value for climate variable
  
  ww <- which(zz < minz,arr.ind=T)  #ocean
  
  if(length(ww) == 0) return( list(clim = climVec, lonLat = lonLat) )
  
  if(nrow(ww) > 5000)ww <- ww[sample(nrow(ww),5000),]
  
  zl     <- cbind(xx[ww[,1]],yy[ww[,2]])
  colnames(zl) <- colnames(lonLat)
  lonLat <- rbind(lonLat,zl)
  
  clim   <- rep(minc,nrow(ww))
  clim   <- c(climVec,clim)
  
  list(clim = clim, lonLat = lonLat, ocean = zl)
}

multivarEmat <- function(bchains,covx,snames,orderB){
  
  #SB - columns in bchains, excludes other
  
  Q <- ncol(covx)
  SM <- length(snames)
  SB <- length(orderB)
  
  onames <- snames[orderB]
  o1     <- paste(onames,'_',sep='')
  o2     <- paste('_',onames,sep='')
  
  bcols <- ccols <- numeric(0)
  ccols <- character(0)
  for(j in 1:SB){
    wj    <- grep(o1[j],colnames(bchains))
    if(length(wj) == 0)wj    <- grep(o2[2],colnames(bchains))
    if(length(wj) == 0)next
    bcols <- c(bcols,min( wj ) )
    ccols <- c(ccols,colnames(bchains)[wj])
  }
  
  cnames <- colnames(bchains[,bcols])
  cnames <- .replaceString(cnames,now='intercept_',new='')
  cnames <- .replaceString(cnames,now='_intercept',new='')
  
  nsim  <- 2000
  ng    <- nrow(bchains)  
  index <- c(round(ng/10,0):ng)
  
  jj <- sample(index,nsim,replace=T)
  rj <- matrix(NA,nsim,SB*SB)
  
  for(j in 1:nsim){
    
    bb     <- matrix( bchains[jj[j],ccols],Q,SB)  #includes omitSpec, but not other
    ss     <- t(bb)%*%covx%*%bb
    rr     <- cov2cor(ss)
    rj[j,] <- as.vector( rr )
  }
  
  tmp    <- apply(rj,2,quantile,c(.5,.025,.975))
  bm     <- matrix(tmp[1,],SB,SB)
  rownames(bm) <- colnames(bm) <- cnames
  loMar  <- which(tmp[2,] < 0 & tmp[3,] < 0)
  hiMar  <- which(tmp[2,] > 0 & tmp[3,] > 0)
  inMar  <- which(tmp[2,] < 0 & tmp[3,] > 0)
  
  ss <- matrix(0,SB,SB)
  ss[inMar] <- 1
  whichZero <- which(ss == 1,arr.ind=T)
  
  list(bm = bm, whichZero = whichZero)
}

gjam_trueVest <- function(chains,true,typeCode,allTypes,xlim=NULL,ylim=NULL,
                          label=NULL,outfile=NULL,colors=NULL,add=F,legend=T){
  
  true   <- as.vector(true)
  ntypes <- length(allTypes)
  
  if(is.null(ylim))ylim <- range(chains,na.rm=T)
  if(is.null(xlim))xlim <- range(true,na.rm=T)
  
  if(!is.matrix(chains)){
    chains <- matrix(chains,ncol=1)
    bCoeffTable <- c(mean(chains),sd(chains),quantile(chains,c(.025,.975)),true)
    bCoeffTable <- matrix(bCoeffTable,1)
  } else {
    bCoeffTable <- processPars(chains,xtrue=true )
  }
  
  if(is.null(colors)){
    colors <- 1
    if(ntypes > 1)colors <- typeCode
  }
  if(length(colors) == 1) colors <- rep(colors,ntypes)
  
  predVsObs(true,p=chains,xlab='true',xlim=xlim,ylim=ylim,ylab='estimated',
            colors=colors,add=add)
  
  if(ntypes > 1 & legend)legend('topleft',allTypes,text.col=colors,bty='n')
  if(!is.null(label))plotLabel(label,above=T)
  
  if(!is.null(outfile))dev.copy2pdf(file=outfile) 
  
  bCoeffTable
}


omitChainCol <- function(cmat,omitCols){
  
  #omitCols - characterVector
  
  keep <- c(1:ncol(cmat))
  ocol <- numeric(0)
  for(j in 1:length(omitCols)){
    ocol <- c(ocol,grep(omitCols[j],colnames(cmat)))
  }
  if(length(ocol) > 0)keep <- keep[-ocol]
  list(keep = keep, omit = ocol)
}


gjam_plotPars <- function(type='CA',y1,yp,censm=NULL){
  
  if(!is.matrix(y1))y1 <- matrix(y1)
  if(!is.matrix(yp))yp <- matrix(yp)
  
  n       <- nrow(y1)
  nk      <- ncol(y1)
  nbin    <- NULL
  nPerBin <- n*nk/15
  breaks  <- NULL
  xlimit  <- range(y1,na.rm=T)
  ylimit  <- quantile(yp,c(.01,.99),na.rm=T)
  vlines  <- NULL
  wide    <- NULL
  MEDIAN  <- T
  LOG     <- F
  yss     <- sd(as.vector(y1))
  
  if(type == 'PA'){
    breaks  <- c(-.5,.5,1.5)
    wide    <- c(.04,.04)
    nPerBin <- NULL
    ylimit  <- xlimit <- c(0,1)
  } 
  if(type == 'OC'){
    breaks  <- seq(-.5,max(y1) + .5,by=1)
    wide    <- 1/max(y1)
    nPerBin <- NULL
    ylimit  <- range(yp,na.rm=T)
    xlimit  <- c( min(floor(y1)), max(ceiling(y1)) )
  } 
  if(type == 'DA')MEDIAN <- F
  if(type %in% c('DA','CA')){
    if(yss > 5){
      xlimit <- range(y1)
      LOG <- T
    }
  }
  if(type %in% c('FC','CC')){
    xlimit <- range(y1)
    MEDIAN <- F
    nPerBin <- round( n*nk/15,0 )
    if(type  == 'CC')LOG <- T
  } 
  if( !is.null(censm) ){
    
    cc  <- censm$partition
    vlines  <- numeric(0)
    breaks  <- NULL
    nPerBin <- n*nk/30
    xlimit  <- range(y1)
    ylimit  <- quantile(yp,c(.01,.99),na.rm=T)
    
    if(ncol(cc) > 1){
      cm     <- unique( as.vector(cc[-1,]) )
      vlines <- cm[is.finite(cm)]
      breaks <- vlines
      nbin   <- nPerBin <- NULL
      uncens <- cbind(cc[3,-ncol(cc)],cc[2,-1])
      wu     <- which( uncens[,1] != uncens[,2] )
      for(m in wu){
        sm <- seq(uncens[m,1],uncens[m,2],length=round(10/length(wu),0))
        if(type == 'DA') sm <- c(uncens[m,1]:uncens[m,2])
        breaks <- c(breaks,sm)
      }
      breaks <- c(breaks,max(y1) + 1)
      breaks <- sort( unique(breaks) )
    }
  }
  
  if(LOG){
    xlimit[1] <- ylimit[1] <- 1
    w0     <- which(y1 == 0)
    y1[w0] <- ylimit[1]
    w0     <- which(yp == 0)
    yp[w0] <- ylimit[1]
    nPerBin <- nPerBin/4
    ylimit[1] <- xlimit[1] <- 1
  }
  
  list( y1 = y1, yp = yp, nbin=nbin, nPerBin=nPerBin, vlines=vlines,
        xlimit=xlimit,ylimit=ylimit,breaks=breaks,wide=wide,LOG=LOG,
        POINTS=F,MEDIAN=MEDIAN )
}

gjam_baselineHist <- function(y1,ylo = 0,nclass=20){
  
  # add histogram to base of current plot
  
  hh <- hist(y1,nclass=nclass,plot=F)
  dx <- diff(hh$breaks)[1]
  yo <- par('usr')[3]
  dy <- diff(par()$usr[3:4])
  hx <- hh$mids
  hy <- hh$density
  hy <- .3*hy*dy/max(hy) + yo
  nh <- length(hy)
  px <- hx[nh] + diff(hy)[nh-1]
  hx <- c(par()$xaxp[1],hx,rep(px,2))
  hy <- c(hy[1],hy,hy[nh],0)
  
  lines(hx,hy+yo,type='s',lwd=4,col='black')
  lines(hx,hy+yo,type='s',lwd=2,col='white')
}

invMatZero <- function(sgibbs,nsim=2000,knames,index=NULL,SS = length(knames)){   # return conditional independence
  
  if(is.null(index))index <- c(1:nrow(sgibbs))
  
  S1 <- sqrt(ncol(sgibbs))
  
  if(is.null(colnames(sgibbs))){
    cnames <- paste('S',c(1:S1),sep='')
    cnames <- outer( cnames,cnames,paste,sep='_')
    colnames(sgibbs) <- cnames
    knames <- cnames
  }
  
#  dnames <- as.vector( outer( knames,knames,paste,sep='_') )
 # gnames <- matrix( unlist(strsplit(colnames(sgibbs),'_')),ncol=2,byrow=T)
 # keepS  <- which(dnames %in% colnames(sgibbs))
  
 # mnames <- gnames[1:S1,1]
 # keep   <- which(mnames %in% knames)
  
 # SS <- length(keep)
  
  jj   <- sample(index,nsim,replace=T)
  sj   <- rj <- matrix(NA,nsim,S1*S1)
  
  for(j in 1:nsim){
    ss     <- matrix(sgibbs[jj[j],],S1,S1) 
    rr     <- cov2cor(ss)
    rj[j,] <- as.vector( rr )
    sj[j,] <- as.vector( solve( ss ) )
  }
  
  tmp    <- apply(rj,2,quantile,c(.5,.025,.975))
  loMar  <- which(tmp[2,] < 0 & tmp[3,] < 0)
  hiMar  <- which(tmp[2,] > 0 & tmp[3,] > 0)
  inMar  <- which(tmp[2,] < 0 & tmp[3,] > 0)
  
  tmp    <- apply(sj,2,quantile,c(.5,.025,.975))
  loCon  <- which(tmp[2,] < 0 & tmp[3,] < 0)
  hiCon  <- which(tmp[2,] > 0 & tmp[3,] > 0)
  inCon  <- which(tmp[2,] < 0 & tmp[3,] > 0)
  
  ss <- matrix(0,S1,S1)
  ss[inCon] <- 1
#  ss[upper.tri(ss)] <- 0
  inConMat <- which(ss == 1,arr.ind=T)
  ss <- matrix(0,S1,S1)
  ss[inMar] <- 1
 # ss[upper.tri(ss)] <- 0
  inMarMat <- which(ss == 1,arr.ind=T)
  
  list( inMar = inMar, inCon = inCon, inMarMat = inMarMat, inConMat = inConMat )
}
  

myMBA.surf <- function (xyz, no.X, no.Y, n = 1, m = 1, h = 8, extend = FALSE, 
                    sp = FALSE, ...) {
  
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  if (!extend) {
    hpts <- chull(xyz[, c(1, 2)])
    hpts <- c(hpts, hpts[1])
  }
  else {
    hpts <- NULL
  }
  x.min.dom <- min(xyz[, 1])
  x.max.dom <- max(xyz[, 1])
  y.min.dom <- min(xyz[, 2])
  y.max.dom <- max(xyz[, 2])
  if ("b.box" %in% elip.args) {
    b.box <- list(...)$b.box
    if (b.box[1] < x.min.dom) {
      x.min.dom <- b.box[1]
    }
    else {
      warning("b.box[1] set to min x value")
    }
    if (b.box[2] > x.max.dom) {
      x.max.dom <- b.box[2]
    }
    else {
      warning("b.box[2] set to max x value")
    }
    if (b.box[3] < y.min.dom) {
      y.min.dom <- b.box[3]
    }
    else {
      warning("b.box[3] set to min y value")
    }
    if (b.box[4] > y.max.dom) {
      y.max.dom <- b.box[4]
    }
    else {
      warning("b.box[4] set to max y value")
    }
  }
  xyz <- as.matrix(xyz)
  storage.mode(xyz) <- "double"
  storage.mode(x.min.dom) <- "double"
  storage.mode(x.max.dom) <- "double"
  storage.mode(y.min.dom) <- "double"
  storage.mode(y.max.dom) <- "double"
  out <- .Call("MBASurf", xyz, as.integer(no.X), as.integer(no.Y), 
               as.integer(m), as.integer(n), as.integer(h), as.integer(extend), 
               as.integer(hpts), 
               x.min.dom, x.max.dom, y.min.dom, y.max.dom)
  if (sp) {
    xy <- expand.grid(out[["x"]], out[["y"]])
    grid <- data.frame(z = matrix(out[["z"]], 
                                  as.integer(length(out[["x"]]) * length(out[["y"]])), 1), 
                       x = xy[, 1], y = xy[, 2])
    coordinates(grid) = ~x + y
    gridded(grid) <- TRUE
  } else {
    grid <- out[c("x", "y", "z")]
  }
  out <- list()
  out$xyz.est <- grid
    out$no.X <- no.X
    out$no.Y <- no.Y
    out$n <- n
    out$m <- m
    out$h <- h
    out$extend <- extend
    out$sp <- sp
    out$b.box <- c(x.min.dom, x.max.dom, y.min.dom, y.max.dom)
    out
}

getCoverageNorm <- function(meanVec,sdVec,values,fraction){
  
  #coverage for values contained in interval fraction
  #  for mean and sd
  
  q  <- -qnorm( (1 - fraction)/2 )
  lo <- meanVec - q*sdVec
  hi <- meanVec + q*sdVec
  ii <- which(values >= lo & values <= hi)
  length(ii)/length(lo)
}


fitSummary <- function(yy,mu,stdev,coverage=.95){
  
  wf <- which(stdev > 0 & is.finite(stdev))
  
  ww <- which(stdev == 0 | !is.finite(stdev))
  stdev[ww] <- min(stdev,na.rm=T)/100
  
  rmse  <-  sqrt( mean((yy - mu)^2,na.rm=T) )
  score <-  mean( getScoreNorm(yy[wf],mu[wf],stdev[wf]^2),na.rm=T )
  cover <- getCoverageNorm(mu,stdev,yy,coverage)
  mvar  <- quantile(stdev^2,.5)
  
  fitTable <- matrix(c(rmse,score,mvar,cover),ncol=1)
  rownames(fitTable) <- c('rmse','score','median var',
                          paste(100*coverage,'% coverage',sep=''))
  signif(fitTable,3)
}

updateWishartNoPrior <- function(xx,yy,df,beta=NULL,IXX=NULL,WX=NULL,WIX=NULL,TRYPRIOR=F){
  
  # more stable without prior
  # TRYPRIOR includes non-informative prior if cholesky fails 
  
  index <- 0
  
  if(is.null(IXX)) IXX <- solve( crossprod(xx) )
  if(is.null(WX))  WX  <- crossprod(xx,yy)
  if(is.null(WIX)) WIX <- IXX%*%WX
  
  SS   <- crossprod(yy) - t(WX)%*%WIX
  testv <- try(chol(SS),T)
  
  if( inherits(testv,'try-error') ){

    message('warning: updateWishartNoPrior')
 #   df <- nrow(SS) + nrow(xx)
    SS <- crossprod(yy - xx%*%beta) +  diag(diag(SS)*.001)*nrow(SS)
    testv <- try(chol(SS),T)
  }
  
  SI <- chol2inv(testv)
  
  testChol <- try(chol(SI),T)
    
  if( inherits(testChol,'try-error') ){
    message('warning: prior used in updateWishartNoPrior')
    if(TRYPRIOR){
      index  <- 1
      SI     <- SI + diag(diag(SI)*.01)
      df     <- nrow(SI) + nrow(xx)
      testChol <- try(chol(SI),T)
    }
  }

  z     <- matrix(rnorm(df*nrow(SS)),df,nrow(SS))%*%testChol
  sinv  <- crossprod(z)
  
  testSolve <- try( solve(sinv),T )
  if( !inherits(testSolve,'try-error') )sigma <- testSolve
  
  if( inherits(testSolve,'try-error') ){
    message('warning: prior used in updateWishartNoPrior')
    if(TRYPRIOR){
      sinv   <- sinv + diag(diag(sinv)*.001)
      df     <- nrow(sinv) + nrow(xx)
      testSolve <- try(chol(sinv),T)
      sigma   <- chol2inv(testSolve)
    }
  }
  
  list( sigma = sigma, sinv = sinv, indicator = index )
}

pointsInBox <- function(xy,xbox,ybox){
  
  which(xy[,1] >= xbox[1] & xy[,1] <= xbox[2] &
        xy[,2] >= ybox[1] & xy[,2] <= ybox[2])
}

interactionsFromGibbs <- function(mainx,bchain,schain,specs,xmnames=names(xx),
                                  xx=colMeans(x),omitY=NULL){
  
  # returns main effects and interactions for variable named main
  # xx are values of covariates to condition on
  # mainx is the name of a main effect
  
 # xmnames <- colnames(x)
  
  if(length(omitY) > 0){
    
    wob <- grep(omitY,colnames(bchain))
    wos <- grep(omitY,colnames(schain))
    
    bchain[,wob] <- 0
    schain[,wob] <- 0
    
    specs <- specs[!specs %in% omitY]
  }
  
  ww   <- grep(':',xmnames)
  int  <- unique( unlist( strsplit(xmnames[ww],':') ) ) 
  int  <- int[int != mainx]
  
  xj <- paste(mainx,specs,sep='_')
  wj <- which(colnames(bchain) %in%  xj)
  if(length(wj) == 0){
    xj <- paste(specs,mainx,sep='_')
    wj <- which(colnames(bchain) %in%  xj)
  }
  
  maine <- bchain[,xj]
  inter  <- maine*0
  
  m1 <- paste(mainx,':',sep='')
  m2 <- paste(':',mainx,sep='')
  i1 <- grep( m1,xmnames )
  i2 <- grep( m2,xmnames )
  if( length(i1) > 0 ){
    ww <- match(unlist( strsplit(xmnames[i1],m1) ),xmnames)
    ox <- xmnames[ww[is.finite(ww)]]
    for(kk in 1:length(i1)){
      xi    <- paste(xmnames[i1[kk]],specs,sep='_')
      wi <- which(colnames(bchain) %in%  xi)
      if(length(wi) == 0){
        xi    <- paste(specs,xmnames[i1[kk]],sep='_')
        wi <- which(colnames(bchain) %in%  xi)
      }
      inter <- inter + bchain[,xi]*xx[ox[kk]]
    }
  }
  if( length(i2) > 0 ){
    ww <- match(unlist( strsplit(xmnames[i2],m2) ),xmnames)
    ox <- xmnames[ww[is.finite(ww)]]
    for(kk in 1:length(i2)){
      xi    <- paste(xmnames[i2[kk]],specs,sep='_')
      wi <- which(colnames(bchain) %in%  xi)
      if(length(wi) == 0){
        xi    <- paste(specs,xmnames[i2[kk]],sep='_')
        wi <- which(colnames(bchain) %in%  xi)
      }
      inter <- inter + bchain[,xi]*xx[ox[kk]]
    }
  }
  list(main = maine, inter = inter)
}


fitLogit <- function(x,y){
  
  p <- ncol(x)
  par0 <- rep(0,p)
  lo   <- par0*0 - 10
  hi   <- par0*0 + 10
  
  out <- nlminb(par0,likeBernLogit,lower=lo,upper=hi)
  b   <- out$par
  L   <- out$objective
  
  list(beta = b, likelihood = L)
}


yaxisHorizLabs <- function( labels, at=c(1:length(labels)), xshift=.05 ){
  
  #add horizontal y axis labels to existing plot
  
  text(par('usr')[3] - xshift*par('usr')[4] - par('usr')[3],y=at,labels,xpd=T)
}


getKernel    <- function(a,d)     exp(-(d/a)^2)           #kernel for distance dij                
getKernSigma <- function(k,s,t) crossprod(t(k))*s + diag(t,nrow(k)) #covariance

updateKernAS <- function(ag,sg,bg,tg,dij,x,
                         amax=100,smax=10){  #update alpha and sigma
  
  asnow <- c(ag,sg)
  asnew <- tnorm.mvt(asnow,asnow,vas,c(0,0),c(amax,smax))
  
  know <- getKernel(ag,dij)	
  knew <- getKernel(asnew[1],dij)
  snow <- getKernSigma(know,sg,tg)
  snew <- getKernSigma(knew,asnew[2],tg)
  
  pnow <- dmvnorm(t(y),know%*%x%*%bg,snow,log=T)
  pnew <- dmvnorm(t(y),knew%*%x%*%bg,snew,log=T)
  
  a <- exp(pnew - pnow)
  z <- runif(1,0,1)
  if(z < a)asnow <- asnew
  asnow
  
}

updateKernB <- function(k,sg,tg,x){
  
  s  <- getKernSigma(k,sg,tg)
  v1 <- t(x)%*%t(k)%*%s
  V  <- solve(v1%*%k%*%x)
  v  <- v1%*%y
  t(myrmvnorm(1,V%*%v,V))
}

updateKernT <- function(bg,sg,k,x,t1=1,t2=1){
  
  e  <- matrix(rnorm(J,0,sqrt(sg)),J,1)
  rs <- y - k%*%(x%*%bg + e)
  u1 <- t1 + n/2
  u2 <- t2 + .5*crossprod(rs)
  1/rgamma(1,u1,u2)
}

plotAspectEffect <- function( betaSlope, slopeRange=c(.1,.3),
                              aspect=seq(-pi,pi,length=100), minEffect=0, maxNumber=20,
                              ylim=NULL, textSize=1, xlab='Aspect', ylab='Effect'){
  
  # slopeRange - two values giving envelop to plot
  # betaSlope  - 3 u coefficients
  # min2Plot   - plot only those having at least this aspect effect
  
  tmp    <- predictSlopeAspect(betaSlope, slopeValue = slopeRange[2],aspect=aspect)
  ub     <- tmp$ubeta
  aspect <- tmp$aspect
  tmp    <- predictSlopeAspect(betaSlope, slopeValue = slopeRange[1], aspect=aspect)
  uc     <- tmp$ubeta
  
  mm  <- apply(ub,2,max)
  
  if(minEffect > 0){
    wm  <- which(mm > minEffect)
    mm  <- mm[wm]
    ub  <- ub[,wm]
    uc  <- uc[,wm]
  }
  if(ncol(ub) > maxNumber){
    o <- order(mm,decreasing=T)
    wm <- o[1:maxNumber]
    ub <- ub[,wm]
    uc <- uc[,wm]
  }
  
  bname <- colnames(ub)
  if(is.null(bname))bname <- paste('S',wm,sep='_')
  
  if(is.null(ylim))ylim <- c(0,max(ub)*1.4)
  
  nr <- ncol(ub)
  
  palette(terrain.colors(ncol(ub)+1))
  
  for(s in 1:nr){
    if(s == 1){
      plot(aspect,ub[,s],type='l',ylim=ylim,xaxt='n',xlab=xlab,ylab=ylab)
      axis(1,at=c(-pi,-pi/2,0,pi/2,pi),labels=c('S','W','N','E','S') )
    }
    polygon(c(aspect,rev(aspect)),c(ub[,s],rev(uc[,s])),col=s,border=s)
  }
  for(s in 1:nr){
    lines(aspect,ub[,s],lwd=2,col=s)
    lines(aspect,uc[,s],lwd=2,col=s)
  }
  for(s in 1:nr){
    wx <- which.max(ub[,s])
    text(aspect[wx],ub[wx,s],bname[s],srt=50,pos=4,col='black',cex=textSize)
  }
}

u2slopeAspect <- function(umat){
  
  # umat - slope, sin(slope)sin(aspect), sin(slope)cos(aspect)
  
  aspect <- atan2(umat[,2],umat[,3])*180/pi
  slope  <- asin(umat[,1])*180/pi
  
  cbind(slope, aspect)
}
  
predictSlopeAspect <- function(betaSlope,slopeValue,aspect=seq(-pi,pi,length=100)){
  
  # betaSlope - matrix of coefficients for u1, u2, u3
  
  nas     <- length(aspect)
  S       <- ncol(betaSlope)
  
  ub <- getSlopeAspect(slopeValue, aspect)
  
  if(ncol(betaSlope) == 3) betaSlope <- t(betaSlope)
  
  ub <- ub%*%betaSlope
  rb <- apply(ub,2,range)
  ub <- ub - matrix(rb[1,],nas,S,byrow=T)
  
  list(ubeta = ub, aspect = aspect)
}


getSlopeAspect <- function(slope,aspect, CLOCKWISE = F){
  #slope   - one value in radians
  #aspect  - vector of values
  #assumes aspect is counterclockwise from North
  
  if(CLOCKWISE)aspect <- 2*pi - aspect
  
  u1 <- sin(slope)
  u2 <- sin(slope)*sin(aspect)
  u3 <- sin(slope)*cos(aspect)
  
  cbind(u1,u2,u3)
}

traitLabel <- function(tname){
  
  tname <- .replaceString(tname,now='soilFactor',new='')
  tname[tname == 'gmPerSeed'] <- 'Seed mass'
  tname[tname == 'gmPerCm']   <- 'Wood dens'
  tname[tname == 'woodSG']    <- 'Wood dens (green)'
  tname[tname == 'maxHt']     <- 'Max ht'
  tname[tname == 'leafN']     <- 'leaf [N]'
  tname[tname == 'leafP']     <- 'leaf [P]'
  tname[tname == "other"]  <- 'Deciduous'
  tname[tname == "broaddeciduous"]  <- 'Deciduous'
  tname[tname == "broadevergreen"]  <- 'BL evergrn'
  tname[tname == "needleevergreen"] <- 'NL evergrn'
  tname[tname == "dioecious"] <- 'Dioecious'
  tname[tname == "u1"] <- 'Slope'
  tname[tname == "u2"] <- 'Aspect 1'
  tname[tname == "u3"] <- 'Aspect 2'
  tname[tname == "ringPorous"] <- 'RP xylem'
  tname[tname == "temp"] <- 'Winter temperature'
  tname[tname == "stdage"] <- 'Stand age'
  for(j in length(tname)){
    tname[j] <- paste(toupper(substring(tname[j], 1, 1)), substring(tname[j], 2),sep = "", collapse = " ")
  }
  tname
}

rangeFromList <- function(xlist){
  
  # get range from a list
  
  xrange <- c(Inf,-Inf)
  for(j in 1:length(xlist)){
    rj <- range( xlist[[j]],na.rm=T )
    if(rj[1] < xrange[1])xrange[1] <- rj[1]
    if(rj[2] > xrange[2])xrange[2] <- rj[2]
  }
  xrange
}

ellipsoidVolume <- function(sigma,alpha=.05){
  
  d     <- nrow(sigma)
  chs   <- qchisq(p = 1 - alpha,d)^(d/2)
  up    <- 2*(pi^(d/2))*sqrt(det(sigma))*chs
  down  <- d*gamma(d/2)
  
  out   <- c(chs,det(sigma),up/down)
  names(out) <- c('chs','det','vol')
  out
}

myMerge <- function(data1,data2,vname){
  
  # merge on character variable vname
  # both must have column names
  
  cols <- union(colnames(data1),colnames(data2))
  qq   <- length(cols)
  
  m1 <- as.character( data1[,vname] )
  m2 <- as.character( data2[,vname] )
  
  m3 <- sort(unique( c(m1,m2) ))
  nn <- length(m3)
  
  w1 <- match(m1,m3)
  w2 <- match(m2,m3)
  c1 <- match(colnames(data1),cols)
  c2 <- match(colnames(data2),cols)
  
  mat <- as.data.frame(matrix(NA,nn,qq))
  colnames(mat) <- cols
  
  for(j in c1){
    tj <- data1[,cols[j]]
    if(is.factor(tj))tj <- as.character(tj)
    mat[w1,j] <- tj
  }
  
  for(j in c2){
    tj <- data2[,cols[j]]
    if(is.factor(tj))tj <- as.character(tj)
    mat[w2,j] <- tj
  }
  mat
}


kmeansEqualGroups <- function( cvar, minSize, maxSize, 
                               fillNA = T, progress=T){
  
  # k means clustering to approx equal size groups
  
  cvar <- as.matrix( cvar )
  if( fillNA ){
    ww <- which( is.na( cvar ), arr.ind = T )
    if( length( ww ) > 0 ){
      cmeans <- colMeans( cvar, na.rm = T )
      cvar[ ww ] <- cmeans[ ww[,2] ]
    }
  }
  n    <- nrow(cvar)
  nall <- maxSize*1.5
  midSize <- mean(c(minSize,maxSize))
  
  if(is.null(rownames(cvar)))rownames(cvar) <- c(1:n)
  
  cid <- rep(1,n)
  cvalues <- 0
  keepRows <- c(1:n)
  nk <- n
  nc <- round(n/midSize,0)
  
  ntot <- n/midSize
  
 # while( nk > nall ){
    
    while( nk > ntot ){
    
    tmp <- kmeans( cvar[keepRows,], nc, iter.max=20 )
    ci <- tmp$cluster 
    
    ig   <- which(tmp$size > minSize & tmp$size < maxSize)
    wg   <- length(ig)
    
    if( wg > 0 ){
      
      nb <- 1 + max(cvalues)
      newVals <- c( nb:(nb + wg - 1) )
      
      mm <- match(ci,ig)
      wm <- which(is.finite(mm))
      
      cid[ keepRows[wm] ] <- newVals[ mm[wm] ]
      
      cvalues <- c(cvalues,newVals)
      keepRows <- keepRows[!keepRows %in% keepRows[wm] ]
      nc <- ceiling( length(keepRows)/midSize )
      
    } else {
      
      maxSize <- maxSize + 1
      minSize <- minSize - 1
    }
    
    nk <- length(keepRows)
    if(progress)print(nk)
  }
  
  allc <- sort(unique(cid))
  
  out <- matrix(NA,length(allc),ncol(cvar))
  rownames(out) <- allc
  colnames(out) <- colnames(cvar)
  
  osd <- out
  
  for(k in 1:ncol(cvar)){
    out[,k] <- tapply(cvar[,k],cid,mean)
    osd[,k] <- tapply(cvar[,k],cid,sd)
  }
  
  list(cluster = cid, clusterMeans = out, clusterSds = osd)
}

colorNeg2Pos <- function(zlevs){
  
  nz <- length(zlevs)
  
  pos <- neg <- F
  
  mid <- findInterval(0,zlevs)
  if(mid < length(zlevs))neg <- T
  
  
  wp <- which(zlevs > 0)
  wn <- which(zlevs <= 0)
  
  np <- length(wp)
  nn <- length(wn)
  
  nm <- max(c(np,nn))
  
  clo    <- colorRampPalette(c('darkblue','blue','green','lightgreen','white'))
  chi    <- colorRampPalette(c('white','yellow','orange','red','brown'))
  
  clof <- clo(nm)
  chif <- chi(nm)
  
  cseq <- numeric(0)
  if(length(wn) > 0)cseq <- clof[ rev(nm - wn) + 1 ]
  if(length(wp) > 0)cseq <- c(cseq,chif[ wp - min(wp) + 1 ] )
  
  cseq
}

speciesNAmerMap <- function(lon,lat,yy,ngrid=100,trim=c(.001,.999),
                            mfrow=c(1,1),
                            POINTS=F,add=F,zlevs=NULL,STATES=F,
                            newPage=T, cRamp=NULL){
  INT <- F
  
  maplon <- xm <- range(lon)
  maplat <- ym <- range(lat)
  zm <- NULL
  
  dscale <- min( diff(xm),diff(ym) )
  
  mdata <- 'world'
  if(dscale < 10){
    INT <- T
    mdata <- 'state'
  }
  if(!STATES)mdata <- 'world'
  
  if(!is.matrix(yy)){
    yy <- matrix(yy,ncol=1)
    colnames(yy) <- ' '
  }
  scaleCoords <- rbind(maplon[1] + diff(maplon)*c(.9,.97),
                       maplat[1] + diff(maplat)*c(.05,.25))

  n <- nrow(yy)
  S <- ncol(yy)
  
  ncol   <- 50
  
  if(is.null(cRamp))cRamp <- c('darkblue','green','orange',
                               'red','brown','black')
  
  colF   <- colorRampPalette(cRamp)
  
  if(!add & newPage)par(mfrow=mfrow, bty='n', mar=c(.1,.1,.1,.1) )
  
  for(j in 1:S){
    
    pj <- colnames(yy)[j]
    ssj <- yy[,j]              #predicted traits
    
 #   data(countyMapEnv)
 #   data(stateMapEnv)
    
    maps::map(boundary=T,col='grey',lwd=2,xlim=range(xm),ylim=range(ym))
    
    sc <- zlevs
    
    if(is.null(zlevs)){
      sc <- seq(min(ssj,na.rm=T),max(ssj,na.rm=T),length=ncol)
      zlevs <- round(seq(ss[1],ss[2],by=diff(ss)/50),4)
    }
    ss <- range(sc)
    zint <- findInterval(sc,zlevs)
                
  #  colq <- colorNeg2Pos(zlevs)
    
  #  colF   <- colorRampPalette(c('white','orange','brown',
  #                               'black','black'))
    colq   <- colF(length(zlevs))
    ncc  <- length(colq)
  
    colz  <- colq[zint]
    
    values2contour(xx=lon,yy=lat,
                   z=ssj,nx=ngrid,ny=ngrid,col=colz,
                   zlevs=zlevs,lwd=4,add=T,fill=T)
    
    zcol <- colq[findInterval(ssj,sc)]
    if(POINTS)symbols(lon,lat,circles=ssj/max(ssj)/2,inches=F,fg=zcol,bg=zcol,add=T)
    
    mapMask(lon,lat,dx=.5*dscale/ngrid,dy=.5*dscale/ngrid,
            whiteOutDist=2*dscale/ngrid,col='white')
    
    map(mdata,add=T,interior=INT,col='white',lwd=5)
    
#    data(countyMapEnv)
#    data(stateMapEnv)
    
    map(mdata,add=T,interior=INT,col='black',lwd=1)
    
    if(max(ss) > 1000)ss <- round(ss,-2)
    if(max(ss) > 100)ss <- round(ss,-1)
    if(max(ss) > 10)ss <- round(ss,0)
    if(max(ss) > 1)ss <- round(ss,1)
    if(max(ss) < 1)ss <- signif(ss,2)
    
    colorLegend(scaleCoords[1,],scaleCoords[2,],
                scale=c(1:ncc),cols=colq,text.col='black',
                labside='left',endLabels=ss)
  }
}

