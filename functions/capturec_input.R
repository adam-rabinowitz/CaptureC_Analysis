require('data.table')
require('purrr')

###############################################################################
## function to extract comparisons
###############################################################################
extract.comparisons <- function(
  path
) {
  # Read in file
  comparisons <- read.table(
    params$comparisons, sep='\t', stringsAsFactors=F
  )
  # check comparisons
  if (
    any(
      substring(comparisons[,1], 1, 3) !=
      substring(comparisons[,2], 1, 3))
  ) {
    stop('unmatched pulldowns')
  }
  # Extract names
  cmp.names <- paste(
    substring(comparisons[,1], 1, 3),
    substring(comparisons[,1], 5),
    substring(comparisons[,2], 5),
    sep='.'
  )
  # Split comparsions by name and return data
  comparisons <- split(comparisons, cmp.names)
  comparisons <- lapply(comparisons, as.character)
  return(comparisons)
}

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
merge.replicates.bait <- function(x) {
  # Check input
  if (class(x) != 'list') {
    stop('input must be a list')
  }
  if (length(x) != 2) {
    stop('list must have a length of two')
  }
  if (any(sapply(x, class) != 'data.frame')) {
    stop('input must be a list of data.frame')
  }
  # Extract conditions and check identity
  cond <- unique(c(unique(x[[1]]$dataset), c(unique(x[[2]]$dataset))))
  if (length(cond) != 1) {
    stop('multiple conditions found')
  }
  # Check replicates are different
  rep1 <- colnames(x[[1]])[grepl('^Rep\\d$', colnames(x[[1]]))]
  rep2 <- colnames(x[[2]])[grepl('^Rep\\d$', colnames(x[[2]]))]
  reps <- c(rep1, rep2)
  if (length(reps) != length(unique(reps))) {
    stop('overlapping replicates found')
  }
  # Merge data and check all fragments are unique 
  merged <- merge(x=x[[1]], y=x[[2]], all=T)
  if (length(merged$fragID) != length(unique(merged$fragID))) {
    stop('Non-unique fragments present in output data')
  }
  # Fill in missing data
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
merge.replicates.all <- function(x, cores) {
  # Check input
  if (class(x) != 'list') {
    stop('input must be a list')
  }
  if (length(x) != 2) {
    stop('list must have a length of two')
  }
  # Check names of all baits
  if (!all(names(x[[1]]) == names(x[[2]]))) {
    stop('different baits present')
  }
  # Transpose baits
  tran.baits <- purrr::transpose(x)
  bait.list <- parallel::mclapply(
    tran.baits,
    merge.replicates.bait,
    mc.cores=cores
  )
  # Return data
  return(bait.list)
}

##############################################################################
## Function to merge view points
##############################################################################
merge.conditions.bait <- function(x) {
  # Check input
  if (class(x) != 'list') {
    stop('input must be a list')
  }
  if (length(x) != 2) {
    stop('list must have a length of two')
  }
  if (any(sapply(x, class) != 'data.frame')) {
    stop('input must be a list of data.frame')
  }
  # Extract conditions and check identity
  cond1 <- unique(x[[1]]$dataset)
  cond2 <- unique(x[[2]]$dataset)
  if (length(cond1) != 1 | length(cond2) != 1) {
    stop('multiple conditions found')
  }
  if (cond1 == cond2) {
    stop('conditions are not unique')
  }
  # Rename columns for first condition
  colnames(x[[1]]) <- gsub('^dataset$', 'cond1', colnames(x[[1]]))
  x1.replicates <- colnames(x[[1]])[grepl('Rep\\d', colnames(x[[1]]))]
  x1.replicates <- paste(cond1, x1.replicates, sep='.')
  colnames(x[[1]])[grepl('Rep\\d', colnames(x[[1]]))] <- x1.replicates
  # Rename columns for second condition
  colnames(x[[2]]) <- gsub('^dataset$', 'cond2', colnames(x[[2]]))
  x2.replicates <- colnames(x[[2]])[grepl('Rep\\d', colnames(x[[2]]))]
  x2.replicates <- paste(cond2, x2.replicates, sep='.')
  colnames(x[[2]])[grepl('Rep\\d', colnames(x[[2]]))] <- x2.replicates
  # Merge data and check all entries are unique
  merged <- merge(x=x[[1]], y=x[[2]], all=T)
  if (length(merged$fragID) != length(unique(merged$fragID))) {
    stop('Non-unique fragments present in output data')
  }
  # Fill in missing data
  merged$cond1 <- cond1
  merged$cond2 <- cond2
  reps <- colnames(merged)[grepl('\\.Rep\\d+$', colnames(merged))]
  for (rep in reps) {
    rep.counts <- merged[,rep,drop=T]
    rep.counts[is.na(rep.counts)] <- 0
    merged[,rep] <- rep.counts
  }
  # Rorder columns and return
  merged <- cbind(
    merged[,c('cond1', 'cond2')],
    merged[,!grepl('^cond\\d+$', colnames(merged))]
  )
  return(merged)
}

###############################################################################
## Function to merge all baits in dataset
###############################################################################
merge.conditions.all <- function(x, cores) {
  # Check input
  if (class(x) != 'list') {
    stop('input must be a list')
  }
  if (length(x) != 2) {
    stop('list must have a length of two')
  }
  # Check names of all baits
  if (!all(names(x[[1]]) == names(x[[2]]))) {
    stop('different baits present')
  }
  # Transpose baits
  tran.baits <- purrr::transpose(x)
  bait.list <- parallel::mclapply(
    tran.baits,
    merge.conditions.bait,
    mc.cores=cores
  )
  # Return data
  return(bait.list)
}
  
##############################################################################
## Read counts for files and merge by condition
##############################################################################
extract.condition.counts <- function(
  input.paths, suffix, min.dist, intra.only, cores
) {
  # Split paths by condition
  condition.paths <- split(
    input.paths,
    gsub(suffix, '', basename(input.paths))
  )
  # Process paths for each condition
  condition.counts <- lapply(
    condition.paths,
    function(paths) {
      # Read in count data
      counts <- lapply(
        paths, 
        read.split.replicates,
        min.dist=min.dist,
        intra.only=intra.only
      )
      # Merge counts if required
      if (length(counts) == 1) {
        counts <- counts[[1]]
      } else if (length(counts) == 2) {
        counts <- merge.replicates.all(
          counts, cores=cores
        )
      } else {
        stop('more than two replicates not considered')
      }
      # Return counts
      return(counts)
    }
  )
  # Rename condition counts and return
  return(condition.counts)
}

