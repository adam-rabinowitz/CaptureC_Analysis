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
  # Extract replicates
  reps <- colnames(bait.data)[grepl('^Rep\\d+$', colnames(bait.data))]
  # Calculate total background count and size
  background.size <- 0 
  background.counts <- c(0, 0)
  for (chrom in names(chr.sizes)) {
    if (chrom == bait.chr) {next}
    background.size <- background.size + chr.sizes[[chrom]]
    background.counts <- background.counts + colSums(
      bait.data[bait.data$fragChr == chrom, reps])}
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
    stop('non-unique bait chromosome found')
  }
  bait.chr.length <- chr.sizes[[bait.chr]]
  # Calculate and add background
  background <- bait.interchr.background(bait.data, chr.sizes)
  bin.background <- background * bin.size
  # Create distance breaks
  max.dist <- max(unlist(chr.sizes)) + 1
  length.out <- ((max.dist - min.dist) / bin.size) + 1
  dist.breaks <- seq(from=min.dist, by=bin.size, length.out=length.out)
  if (tail(dist.breaks, 1) < max.dist) {
    stop('error in generating breaks')
  }
  dist.bins <- cut(abs(bait.data$fragDist), breaks=dist.breaks, right=F)
  # Create distance bins and create output dat frame
  dist.df <- data.frame(
    'binStart' = head(dist.breaks, -1),
    'binEnd' = tail(dist.breaks, -1) -1
  )
  # Add data for each replicate to dataframe
  reps <- colnames(bait.data)[grepl('^Rep\\d+$', colnames(bait.data))]
  for (rep in reps) {
    # Count reads in each distance bin
    rep.counts <- tapply(bait.data[,rep], dist.bins, sum)
    # Replace NA with 0 for distances less than chromosome length
    rep.counts[is.na(rep.counts) & dist.df$binStart < bait.chr.length] <- 0
    # Add background and calculate bin ratio
    rep.counts <- rep.counts + bin.background[rep]
    rep.ratio <- rep.counts / sum(rep.counts, na.rm=T)
    # Store data
    dist.df[,rep] <- rep.ratio
  }
  # Convert to matrix and return
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
    mc.cores=cores
  )
  # Merge frequency data and convert to dataframe
  frequency.data <- apply(simplify2array(frequency.list), c(1,2), mean, na.rm=T)
  frequency.data <- as.data.frame(frequency.data)
  # Extract data for fit
  binCentre <- (frequency.data$binStart + frequency.data$binEnd) / 2
  xlimits <- c(min(frequency.data$binStart), max(frequency.data$binEnd))
  # Perform fits and return data
  reps <- colnames(frequency.data)[grepl('^Rep\\d+$', colnames(frequency.data))]
  fits <- mclapply(
    frequency.data[,reps],
    distance.decay.monospline,
    x=binCentre,
    limits=xlimits,
    k=k,
    mc.cores=min(cores, length(reps))
  )
  # Return fits
  return(fits)
}

##############################################################################
## Calculates fits for given dataset
##############################################################################
distance.decay.fit <- function(
  path, min.dist, bin.size, k, chr.sizes, cores
) {
  # Read in data and extract values
  frag.data <- fread(path, showProgress=F)
  dataset <- unique(frag.data$dataset)
  if (length(dataset) != 1) {
    stop('multiple datasets found')}
  reps <- colnames(frag.data)[grepl('^Rep\\d+$', colnames(frag.data))]
  # Split data by bait and convert to dataframe
  bait.list <- split(frag.data, frag.data$baitID)
  bait.list <- lapply(bait.list, as.data.frame)
  # Extract proabilites for each bait
  fits <- create.decay.fit(
    bait.list=bait.list,
    min.dist=min.dist,
    bin.size=bin.size,
    k=k,
    chr.sizes=chr.sizes,
    cores=cores
  )
  # name and return fits
  names(fits) <- paste(dataset, reps, sep='.')
  return(fits)
}

###############################################################################
## create distance decay fit for all datasets
###############################################################################
distance.decay.fit.all <- function(
  paths, min.dist, bin.size, k, chr.sizes, cores
) {
  # Extract all fits
  replicate.fits <- lapply(
    paths,
    distance.decay.fit,
    min.dist=params$mindist,
    bin.size=params$binsize,
    k=params$k,
    chr.sizes=chr.sizes,
    cores=cores
  )
  # flatten fits and return
  replicate.fits <- do.call(c, replicate.fits)
  return(replicate.fits)
}





















