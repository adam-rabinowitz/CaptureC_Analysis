require('data.table')
require('parallel')
require('mgcv')

###############################################################################
## Calculate intrachromosomal background for a single bait
###############################################################################
bait.interchr.background <- function(
  bait.data, chr.sizes
) {
  # Extract bait chromosome
  bait.chr <- unique(bait.data$baitChr)
  if (length(bait.chr) != 1) {
    stop('non-unique bait chromosome found')}
  # Calculate total background count and size
  background.size <- 0 
  background.counts <- c(0, 0)
  for (chrom in names(chr.sizes)) {
    if (chrom == bait.chr) {next}
    background.size <- background.size + chr.sizes[[chrom]]
    background.counts <- background.counts + colSums(
      bait.data[bait.data$fragChr == chrom, c('N1', 'N2')])}
  # Calculate, check and return background
  background <- background.counts / background.size
  if (min(background) <= 0) {
    stop('no background detetcted')}
  return(background)
}

###############################################################################
## Calculate binned distance-frequency profile
###############################################################################
bait.distance.frequency <- function(
  bait.data, min.dist, bin.size, chr.sizes
) {
  # Extract bait chromosome data
  bait.chr <- unique(bait.data$baitChr)
  if (length(bait.chr) != 1) {
    stop('non-unique bait chromosome found')}
  bait.chr.length <- chr.sizes[[bait.chr]]
  # Create distance breaks
  max.dist <- max(unlist(chr.sizes)) + 1
  length.out <- ((max.dist - min.dist) / bin.size) + 1 
  dist.breaks <- seq(from=min.dist, by=bin.size, length.out=length.out)
  if (tail(dist.breaks, 1) < max.dist) {
    stop('error in generating breaks')}
  # Find distance groups for counts
  dist.bins <- cut(abs(bait.data$fragDist), breaks=dist.breaks, right=F)
  dist.df <- data.frame(
    'binStart' = head(dist.breaks, -1),
    'binEnd' = tail(dist.breaks, -1) -1,
    'N1' = tapply(bait.data$N1, dist.bins, sum),
    'N2' = tapply(bait.data$N2, dist.bins, sum))
  # Set intrachromosomal NA counts to zero
  intra.filter <- matrix(
    rep(dist.df$binStart < bait.chr.length, times=4),
    ncol=4)
  dist.df[is.na(dist.df) & intra.filter] <- 0
  # Calculate and add background
  background <- bait.interchr.background(bait.data, chr.sizes)
  bin.background <- background * bin.size
  dist.df$N1 <- dist.df$N1 + bin.background['N1']
  dist.df$N2 <- dist.df$N2 + bin.background['N2']
  # Convert counts to ratio
  dist.df$N1 <- dist.df$N1 / sum(dist.df$N1, na.rm=T)
  dist.df$N2 <- dist.df$N2 / sum(dist.df$N2, na.rm=T)
  # Generate matrix and return
  dist.matrix <- as.matrix(dist.df)
  return(dist.matrix)
}

###############################################################################
## Function to generate monotonic spline
###############################################################################
distance.decay.monospline <- function(x, y, limits, k) {
  # Calculate knots
  knots <- data.frame(
    'x' = 10^(
      seq(from=log10(limits[1]), to=log10(limits[2]), length.out=k)))
  # perform spline and save terms
  sfit <- mgcv::gam(y~s(x,bs="cr"))
  # Generate smoothing terms
  dat <- data.frame('x'=x, 'y'=y)
  knots <- data.frame('x'=knots)
  smooth.terms <- mgcv::smoothCon(
    s(x, k=k, bs="cr"),
    dat,
    knots=knots)[[1]]
  # Set constraints
  constraints <- mgcv::mono.con(
    smooth.terms$xp,
    up=FALSE) # Set to true for increase 
  # Set up list of variable for spline
  M <- list(
    X=smooth.terms$X,
    C=matrix(0,0,0),
    sp=sfit$sp,
    p=smooth.terms$xp * -1, # keep positive for increase
    y=y,
    w=y*0+1,
    Ain=constraints$A,
    bin=constraints$b,
    S=smooth.terms$S,
    off=0)
  # Perform smooth spline and return data
  penalties <- mgcv::pcls(M)
  mono.fit <- list(
    'smooth' = smooth.terms,
    'penalties' = penalties,
    'knots' = knots$x)
  return(mono.fit)
}

##############################################################################
## Create fit for series of baits
##############################################################################
create.decay.fit <- function(
  bait.list, min.dist, bin.size, k, chr.sizes, cores
) {
  # Calculate frequency for bins
  frequency.list <- mclapply(
    bait.list,
    bait.distance.frequency,
    min.dist=min.dist,
    bin.size=bin.size,
    chr.sizes=chr.sizes,
    mc.cores=cores)
  # Merge frequency data and convert to dataframe
  frequency.data <- apply(simplify2array(frequency.list), c(1,2), mean, na.rm=T)
  frequency.data <- as.data.frame(frequency.data)
  # Extract data for fit
  binCentre <- (frequency.data$binStart + frequency.data$binEnd) / 2
  xlimits <- c(min(frequency.data$binStart), max(frequency.data$binEnd))
  # Perform fits and return data
  fits <- mclapply(
    frequency.data[,c('N1', 'N2')],
    distance.decay.monospline,
    x=binCentre,
    limits=xlimits,
    k=k,
    mc.cores=min(cores, 2))
  # name and return fits
  return(fits)
}

##############################################################################
## Calculates fits for given dataset
##############################################################################
distance.decay.fit <- function(
  path, min.dist, bin.size, k, chr.sizes, cores
) {
  # Read in data and check dataset
  frag.data <- fread(path, showProgress=F)
  dataset <- unique(frag.data$dataset)
  if (length(dataset) != 1) {
    stop('multiple datasets found')}
  # Extract proabilites for each bait
  bait.list <- split(frag.data, frag.data$baitID)
  fits <- create.decay.fit(
    bait.list=bait.list,
    min.dist=min.dist,
    bin.size-bin.size,
    k=k,
    chr.sizes=chr.sizes,
    cores=cores
  )
  # name and return fits
  names(fits) <- paste0(dataset, c('.Rep1', '.Rep2'))
  return(fits)
}




























