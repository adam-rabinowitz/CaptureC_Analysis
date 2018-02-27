require('data.table')

##############################################################################
## Function to replace NA values in one vector with values from other vector
##############################################################################
fillNA <- function(v1, v2) {
  indices <- is.na(v1)
  v1[indices] <- v2[indices]
  return(v1)
}

##############################################################################
## Function to replace NA values in one vector with zero
##############################################################################
zeroNA <- function(v) {
  indices <- is.na(v)
  v[indices] <- 0
  return(v)
}

###############################################################################
## Read in count data
###############################################################################
read.split.replicates <- function(path, min.dist, intra.only) {
  # Check arguments
  if (class(intra.only) != 'logical') {
    stop('intra.only argument must be logical')}
  if (class(min.dist) != 'numeric') {
    stop('min.dist argument must be numeric')}
  if (min.dist < 0) {
    stop('min.dist cannot be negatove')}
  # Read in data and remove proximal reads
  data <- data.table::fread(path, showProgress=F)
  data <- data[abs(data$fragDist) >= min.dist,]
  # Remove intrachromosomal counts if requested
  if (intra.only) {
    data <- data[data$baitChr == data$fragChr,]
  }
  # Split data by bait and return
  data.list <- base::split(data, data$baitID)
  return(data.list)
}

##############################################################################
## Function to merge view points
##############################################################################
merge.bait.data <- function(x1, x2) {
  # Extract data and merge
  cond1 <- x1$dataset[1]
  cond2 <- x2$dataset[1]
  merged <- merge(x=x1, y=x2, by='fragID', all=T)
  # Create output
  merged.df <- data.frame(
    'cond1' = cond1,
    'cond2' = cond2,
    'baitName' = fillNA(merged$baitName.x, merged$baitName.y),
    'baitID' = fillNA(merged$baitID.x, merged$baitID.y),
    'baitChr' = fillNA(merged$baitChr.x, merged$baitChr.y),
    'baitStart' = fillNA(merged$baitStart.x, merged$baitStart.y),
    'baitEnd' = fillNA(merged$baitEnd.x, merged$baitEnd.y),
    'fragID' = merged$fragID,
    'fragChr' = fillNA(merged$fragChr.x, merged$fragChr.y),
    'fragStart' = fillNA(merged$fragStart.x, merged$fragStart.y),
    'fragEnd' = fillNA(merged$fragEnd.x, merged$fragEnd.y),
    'fragDist' = fillNA(merged$fragDist.x, merged$fragDist.y),
    'cond1.1' = zeroNA(merged$N1.x),
    'cond1.2' = zeroNA(merged$N2.x),
    'cond2.1' = zeroNA(merged$N1.y),
    'cond2.2' = zeroNA(merged$N2.y),
    stringsAsFactors=F)
  # Order and return output
  merged.df <- merged.df[order(merged.df$baitID, merged.df$fragID),]
  return(merged.df)
}

###############################################################################
## Function to merge all baits in dataset
###############################################################################
merge.datasets <- function(d1, d2, cores) {
  # Find common baits
  baits <- union(names(d1), names(d2))
  # Merge data for each bait
  merged.data <- mclapply(
    baits,
    function(bait) {
      merge.bait.data(d1[[bait]], d2[[bait]])},
    mc.cores=cores)
  return(merged.data)
}
  



