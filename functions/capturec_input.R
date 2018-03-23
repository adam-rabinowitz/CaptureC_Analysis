require('data.table')
require('plyr')

# ##############################################################################
# ## Function to replace NA values in one vector with values from other vector
# ##############################################################################
# fillNA <- function(v1, v2) {
#   indices <- is.na(v1)
#   v1[indices] <- v2[indices]
#   return(v1)
# }
# 
# ##############################################################################
# ## Function to replace NA values in one vector with zero
# ##############################################################################
# zeroNA <- function(v) {
#   indices <- is.na(v)
#   v[indices] <- 0
#   return(v)
# }

# ##############################################################################
# ## Function to merge view points
# ##############################################################################
# merge.bait.data <- function(x1, x2) {
#   # Extract data and merge
#   cond1 <- x1$dataset[1]
#   cond2 <- x2$dataset[1]
#   merged <- merge(x=x1, y=x2, by='fragID', all=T)
#   # Create output
#   merged.df <- data.frame(
#     'cond1' = cond1,
#     'cond2' = cond2,
#     'baitName' = fillNA(merged$baitName.x, merged$baitName.y),
#     'baitID' = fillNA(merged$baitID.x, merged$baitID.y),
#     'baitChr' = fillNA(merged$baitChr.x, merged$baitChr.y),
#     'baitStart' = fillNA(merged$baitStart.x, merged$baitStart.y),
#     'baitEnd' = fillNA(merged$baitEnd.x, merged$baitEnd.y),
#     'fragID' = merged$fragID,
#     'fragChr' = fillNA(merged$fragChr.x, merged$fragChr.y),
#     'fragStart' = fillNA(merged$fragStart.x, merged$fragStart.y),
#     'fragEnd' = fillNA(merged$fragEnd.x, merged$fragEnd.y),
#     'fragDist' = fillNA(merged$fragDist.x, merged$fragDist.y),
#     'cond1.1' = zeroNA(merged$N1.x),
#     'cond1.2' = zeroNA(merged$N2.x),
#     'cond2.1' = zeroNA(merged$N1.y),
#     'cond2.2' = zeroNA(merged$N2.y),
#     stringsAsFactors=F)
#   # Order and return output
#   merged.df <- merged.df[order(merged.df$baitID, merged.df$fragID),]
#   return(merged.df)
# }

###############################################################################
## Function to calculate distances
###############################################################################
calculate.distance <- function(
  data
) {
  # Calculate distance
  distance <- pmax(
    (data$baitStart - data$fragEnd) - 1,
    (data$fragStart - data$baitEnd) - 1
  )
  # Adjust sign of distance
  distance <- distance * sign(data$fragStart - data$baitStart)
  distance[data$baitChr != data$fragChr] <- Inf
  return(distance)
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
  # Read in data, recalculate distance and remove proximal reads
  data <- data.table::fread(path, showProgress=F)
  data$fragDist <- calculate.distance(data)
  data <- data[abs(data$fragDist) >= min.dist,]
  # Remove intrachromosomal counts if requested
  if (intra.only) {
    data <- data[data$baitChr == data$fragChr,]
  }
  # Split data by bait and return
  data.list <- base::split(data, data$baitID)
  data.list <- lapply(data.list, as.data.frame)
  return(data.list)
}

##############################################################################
## Function to merge replicates of a single bait
##############################################################################
merge.replicates.bait <- function(x1, x2) {
  # Extract conditions and check identity
  cond <- unique(c(unique(x1$dataset), c(unique(x2$dataset))))
  if (length(cond) != 1) {
    stop('multiple conditions found')
  }
  # Check replicates are different
  rep1 <- colnames(x1)[grepl('^Rep\\d$', colnames(x1))]
  rep2 <- colnames(x2)[grepl('^Rep\\d$', colnames(x2))]
  reps <- c(rep1, rep2)
  if (length(reps) != length(unique(reps))) {
    stop('overlapping replicates found')
  }
  # Merge data
  merged <- merge(x=x1, y=x2, all=T)
  for (rep in reps) {
    rep.counts <- merged[,rep,drop=T]
    rep.counts[is.na(rep.counts)] <- 0
    merged[,rep] <- rep.counts
  }
  # Return merged
  return(merged)
}

##############################################################################
## Function to merge replicates of a single bait
##############################################################################
merge.replicates.all <- function(x1, x2) {
  # Extract names of all baits
  baits <- union(names(x1), names(x2))
  # Transpose baits
  bait.list <- list()
  for (b in baits) {
    bait.list[[b]] <- merge.replicates.bait(x1[[b]], x2[[b]])
  }
  return(bait.list)
}

##############################################################################
## Read counts for files and merge by condition
##############################################################################
extract.condition.counts <- function(input.paths, suffix) {
  # Split paths by condition
  condition.paths <- split(
    input.paths,
    gsub(suffix, '', basename(input.paths))
  )
  # Process paths for each condition
  condition.counts <- lapply(
    condition.paths,
    function(paths) {
      print(paths)
      # Read in count data
      counts <- lapply(
        paths, 
        read.split.replicates,
        min.dist=params$mindist,
        intra.only=T
      )
      # Merge counts if required
      if (length(counts) == 1) {
        counts <- counts[[1]]
      } else if (length(counts) == 2) {
        counts <- merge.replicates.all(
          counts[[1]], counts[[2]]
        )
      } else {
        stop('more than two replicates not considered')
      }
      # Return counts
      return(counts)
    }
  )
}

##############################################################################
## Function to merge view points
##############################################################################
merge.conditions.bait <- function(x1, x2) {
  # Extract conditions and check identity
  cond1 <- unique(x1$dataset)
  cond2 <- unique(x2$dataset)
  if (length(cond1) != 1 | length(cond2) != 1) {
    stop('multiple conditions found')
  }
  if (cond1 == cond2) {
    stop('conditions are not unique')
  }
  # Rename columns
  x1 <- plyr::rename(x1, replace=c('dataset'='cond1'))
  x1.replicates <- colnames(x1)[grepl('Rep\\d', colnames(x1))]
  x1.replicates <- paste(cond1, x1.replicates, sep='.')
  colnames(x1)[grepl('Rep\\d', colnames(x1))] <- x1.replicates
  x2 <- plyr::rename(x2, replace=c('dataset'='cond2'))
  x2.replicates <- colnames(x2)[grepl('Rep\\d', colnames(x2))]
  x2.replicates <- paste(cond2, x2.replicates, sep='.')
  colnames(x2)[grepl('Rep\\d', colnames(x2))] <- x2.replicates
  # Merge data
  merged <- merge(x=x1, y=x2, all=T)
  # Fill in missing data
  merged$cond1 <- cond1
  merged$cond2 <- cond2
  reps <- colnames(merged)[grepl('\\.Rep\\d', colnames(merged))]
  for (rep in reps) {
    rep.counts <- merged[,rep,drop=T]
    rep.counts[is.na(rep.counts)] <- 0
    merged[,rep] <- rep.counts
  }
  # Rorder columns and return
  merged <- cbind(
    merged[,c('cond1', 'cond2')],
    merged[,!grepl('cond\\d', colnames(merged))]
  )
  return(merged)
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
  



