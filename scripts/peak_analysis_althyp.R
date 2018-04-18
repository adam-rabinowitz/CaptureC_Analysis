require('data.table')
require('GenomicRanges')

###############################################################################
## Read input files
###############################################################################
read.filter.file <- function(path, padj.thr) {
  # Process files dependent on compression
  if (endsWith(path, '.gz')) {
    command <- paste('zcat', path)
  } else {
    command <- path
  }
  # Read in data, filter and return
  data <- fread(command, showProgress=F, stringsAsFactors=F)
  filt.data <- data[data$padj < padj.thr,]
  filt.data <- as.data.frame(filt.data)
  return(as.data.frame(filt.data))
}

###############################################################################
## Function to calculate distance between neighbouring fragments
###############################################################################
neighbour.distance <- function(start, end) {
  # Calculate upstream distance
  distance <- (tail(start, -1) - head(end, - 1)) -1
  # Check for unexpected overlapping intervals and return data
  if (any(distance < 0)) {
    stop('regions overlap')
  }
  return(distance)
}

###############################################################################
## Function to split results into list elements seperated by 2 * window size
###############################################################################
split.window <- function(x, window) {
  # Calculate and store neighbour distance
  ndistance <- neighbour.distance(x$fragStart, x$fragEnd)
  # Calculate groups seperated by further than the window
  groups <- cumsum(ndistance > window)
  groups <- c(1, groups + 1)
  groups <- factor(groups, levels=1:tail(groups, 1))
  # Split and return input data
  x.list <- split(x, groups)
  return(x.list)
}

###############################################################################
## Function to split peaks overlapping with specified fragments
###############################################################################
split.overlapping.peaks <- function(peak, frag.ids) {
  # Find overlaps
  overlap.no <- length(intersect(peak$fragID, frag.ids))
  # Return inut as list if no overlaps found
  if (overlap.no == 0) {
    output <- list('1' = peak)
    # Return empty list if no non-overlaps found
  } else if (overlap.no == nrow(peak)) {
    output <- list('1' = NULL)
    # Split input by overlaps
  } else {
    # Find overlaps
    overlap <- peak$fragID %in% frag.ids
    overlap.groups <- cumsum(peak$fragID %in% frag.ids)
    # Remove overlaps
    peak.filt <- peak[!overlap,]
    overlap.groups.filt <- overlap.groups[!overlap]
    # Create output
    output <- split(peak.filt, overlap.groups.filt)
    names(output) <- 1:length(output)
  }
  # Return data
  return(output)
}

###############################################################################
## Function to split datasets by window and overlap
###############################################################################
generate.peaks.bait <- function(
  bait.name, change.bait.list, change, window, sig.padj
) {
  # Extract bait data from list
  bait.list <- lapply(change.bait.list, function(z) {z[[bait.name]]})
  # Extract target data
  target.data <- bait.list[[change]]
  if (class(target.data) != 'data.frame') {
    stop('target.data must be a data.frame')
  }
  if (nrow(target.data) == 0) {
    stop('target.data data.frame is empty')
  }  
  # Extract frag ids to filter from target
  filter.data <- bait.list[setdiff(names(change.bait.list), change)]
  filter.data <- filter.data[!sapply(filter.data, is.null)]
  if (length(filter.data) == 0) {
    filter.frag.ids = numeric(0)
  } else {
    filter.frag.ids <- lapply(filter.data, function(z) {z$fragID})
    filter.frag.ids <- unique(do.call(c, filter.frag.ids))
  }
  # Split target by window size
  target.windows <- split.window(
    target.data,
    window=window
  )
  # Further split peaks by overlap
  target.unique <- lapply(
    target.windows,
    split.overlapping.peaks,
    frag.ids=filter.frag.ids
  )
  target.unique <- do.call(c, target.unique)
  target.unique <- target.unique[!sapply(target.unique, is.null)]
  # Filter peaks not containing fragments below threshold
  if (length(target.unique) == 0) {
    output <- target.unique
  } else {
    min.padj <- sapply(target.unique, function(z) {min(z$padj)})
    output <- target.unique[min.padj < sig.padj]
  }
  # Rename output and return
  if (length(output) > 0) {
    names(output) <- 1:length(output)
    output <- lapply(
      output, function(z) {
        z$change = change
        return(z)
      }
    )
  }
  return(output)
}

###############################################################################
## Function to check list overlaps
###############################################################################
check.list.overlaps <- function(
  change.list, sig.padj
) {
  # Check overlaps
  frag.ids <- unlist(lapply(change.list, function(z) {z$fragID}))
  duplicate.ids <- unique(frag.ids[duplicated(frag.ids)])
  # Loop through duplicates
  for (dup in duplicate.ids) {
    dup.list <- lapply(change.list, function(z) {z[z$fragID == dup,]})
    dup.list <- dup.list[sapply(dup.list, nrow) > 0]
    dup.df <- do.call(rbind, dup.list)
    if (sum(dup.df$padj < sig.padj) > 1) {
      message('DANGER! duplicate significant fragment found')
      warning('duplicate significant fragment found')
    } else {
      message('duplicate fragment found')
    }
    print(dup.df)
  }
}

###############################################################################
## Function to find all pssible peaks from two files
###############################################################################
generate.peaks.all <- function(
  diff.data, same.data, max.padj, sig.padj, window
) {
  # Create list containing changes
  change.list <- list(
    'incr' = diff.data[diff.data$lfc > 0 & diff.data$padj < max.padj,],
    'decr' = diff.data[diff.data$lfc < 0 & diff.data$padj < max.padj,],
    'same' = same.data[same.data$padj < max.padj,]
  )
  # Check overlaps
  check.list.overlaps(change.list, sig.padj)
  # Check datasets are from a common comparison
  dataset <- unique(unlist(
    lapply(change.list, function(z) {unique(z$dataset)})))
  if (length(dataset) != 1) {
    stop('files contain multiple datasets')
  }
  cond1 <- unique(unlist(
    lapply(change.list, function(z) {unique(z$cond1)})))
  if (length(cond1) != 1) {
    stop('files contain multiple cond1')
  }
  cond2 <- unique(unlist(
    sapply(change.list, function(z) {unique(z$cond2)})))
  if (length(cond2) != 1) {
    stop('files contain multiple cond2')
  }
  # Split datasets by bait
  change.bait.list <- lapply(
    change.list,
    function(z) {
      split(z, z$baitName)
    }
  )
  # Loop through all baits for all groups
  peaks <- list()
  for (change in names(change.bait.list)) {
    for (bait in names(change.bait.list[[change]])) {
      # Extract bait data
      peak.group <- paste(dataset, cond1, cond2, bait, change, sep='.')
      peak.data <- generate.peaks.bait(
        bait.name=bait, change.bait.list=change.bait.list, change=change,
        window=window, sig.padj=sig.padj
      )
      peaks[[peak.group]] <- peak.data
    }
  }
  # Concatenate and return peaks
  peaks <- do.call(c, peaks)
  return(peaks)
}

###############################################################################
## Extract output data for peaks
###############################################################################
create.peak.df <- function(
  peak, sig.padj
) {
  # Create output df
  peak.df <- data.frame(
    'dataset' = unique(peak$dataset),
    'normalisation' = unique(peak$normalisation),
    'cond1' = unique(peak$cond1),
    'cond2' = unique(peak$cond2),
    'change' = unique(peak$change),
    'baitName' = unique(peak$baitName),
    'baitID' = unique(peak$baitID),
    'baitStart' = unique(peak$baitStart),
    'baitEnd' = unique(peak$baitEnd),
    'peakName' = NA,
    'peakChr' = unique(peak$fragChr),
    'peakStart' = min(peak$fragStart),
    'peakEnd' = min(peak$fragEnd),
    'peakFragStart' = min(peak$fragID),
    'peakFragEnd' = max(peak$fragID),
    'peakFrags' = nrow(peak),
    'peakSigFrags' = sum(peak$padj < sig.padj),
    'peakMeanCount' = mean(peak$repMean),
    'padjMin' = min(peak$padj),
    'lfcMean' = mean(peak$lfc),
    'lfcWeightedMean' = weighted.mean(
      x=peak$lfc,
      w=-log10(peak$padj)
    ),
    'lfcAtMin' = peak$lfc[base::which.min(peak$pvalue)]
  )
  # Check and return data
  if (nrow(peak.df) > 1) {
    stop('non-unique values present')
  }
  return(peak.df)
}

###############################################################################
## Extract all peaks in a file
###############################################################################
generate.peaks.from.file <- function(
  sig.file, diff.hyp, same.hyp, normalisation, max.padj, sig.padj, window
) {
  # Read in significant fragments
  sig.fragments <- read.table(
    sig.file, header=T, sep='\t', stringsAsFactors=F
  )
  # Filter data to remove unwanted hypotheses and normalisation
  hypotheses <- paste(sig.fragments$altHyp, sig.fragments$lfcThr, sep='.')
  sig.fragments <- sig.fragments[
    sig.fragments$normalisation == normalisation &
    (hypotheses == diff.hyp | hypotheses == same.hyp)
  ,]
  # Split dataset by dataset, conditions and bait
  groups <- paste(
    sig.fragments$dataset,
    sig.fragments$cond1,
    sig.fragments$cond2,
    sig.fragments$baitName
  )
  sig.fragments.list <- split(sig.fragments, groups)
  # Loop through groups and find peaks
  peak.list <- list()
  for (dataset in names(sig.fragments.list)) {
    # Extract data for dataset
    dataset.data <- sig.fragments.list[[dataset]]
    dataset.hypotheses <- paste(
      dataset.data$altHyp, dataset.data$lfcThr, sep='.')
    # Seperate differential and same regions
    diff.data <- dataset.data[dataset.hypotheses == diff.hyp,]
    same.data <- dataset.data[dataset.hypotheses == same.hyp,]
    # Generate peaks and store
    peaks <- generate.peaks.all(
      diff.data=diff.data, same.data=same.data, max.padj=max.padj,
      sig.padj=sig.padj, window=window
    )
    peak.list[[dataset]] <- peaks
  }
  # Flatten peak list and rename
  peak.names <- unlist(lapply(peak.list, names))
  peak.list <- do.call(c, peak.list)
  # Format peaks, convert to datframe and return
  peak.df.list <- lapply(peak.list, create.peak.df, sig.padj=sig.padj)
  peak.df <- as.data.frame(data.table::rbindlist(peak.df.list))
  peak.df$peakName <- peak.names
  return(peak.df)
}

###############################################################################
## Function to find overlaps
###############################################################################
annotate.overlaps <- function(peaks) {
  # Create peak genomic range
  peak.list <- split(peaks, peaks$dataset)
  for (dataset.no in 1:length(peak.list)) {
    # Create genomic range
    dataset <- peak.list[[dataset.no]]
    dataset.gr <- GenomicRanges::makeGRangesFromDataFrame(
      dataset,
      seqnames.field='peakChr',
      start.field='peakStart',
      end.field='peakEnd',
      ignore.strand=T,
      keep.extra.columns=T
    )
    # Find overlaps by peak name
    overlaps <- findOverlaps(dataset.gr, dataset.gr)
    overlaps <- overlaps[queryHits(overlaps) != subjectHits(overlaps)]
    query.names <- mcols(dataset.gr)$peakName[queryHits(overlaps)]
    subject.names <- mcols(dataset.gr)$peakName[subjectHits(overlaps)]
    overlap.list <- split(subject.names, query.names)
    # Add overlap peaks to datasets
    overlappingPeakNo <- sapply(
      overlap.list, length
    )
    dataset$overlappingPeakNo <- overlappingPeakNo[dataset$peakName]
    dataset$overlappingPeakNo[is.na(dataset$overlappingPeakNo)] <- 0
    overlappingPeakNames <- sapply(
      overlap.list, paste, collapse=';'
    )[dataset$peakName]
    dataset$overlappingPeakNames <- overlappingPeaks[dataset$peakName]
    # Add overlapping comparisons to dataset
    overlappingComparisons <- sapply(
      overlap.list,
      function(z) {
        comparisons <- sapply(
          strsplit(z, '\\.'),
          function(z) {
            paste(z[1:3], collapse='.')
          }
        )
        comparisons <- unique(comparisons)
        return(comparisons)
      }
    )
    dataset$overlappingComparisonNo <- sapply(
      overlappingComparisons, length
    )[dataset$peakName]
    dataset$overlappingComparisonNo[is.na(dataset$overlappingComparisonNo)] <- 0
    dataset$overlappingComparisonNames <- sapply(
      overlappingComparisons, paste, collapse=';'
    )[dataset$peakName]
    dataset$overlappingComparisonNames
    # Add overlapping baits
    overlappingBaits <- sapply(
      overlap.list,
      function(z) {
        baits <- sapply(
          strsplit(z, '\\.'),
          function(z) {
            z[4]
          }
        )
        baits <- unique(baits)
        return(baits)
      }
    )
    dataset$overlappingBaitNo <- sapply(
      overlappingBaits, length
    )[dataset$peakName]
    dataset$overlappingBaitNo[is.na(dataset$overlappingBaitNo)] <- 0
    dataset$overlappingBaitNames <- sapply(
      overlappingBaits, paste, collapse=';'
    )[dataset$peakName]
    # Store data
    peak.list[[dataset.no]] <- dataset
  }
  # Join peak list and return
  peaks <- do.call(rbind, peak.list)
  return(peaks)
}

###############################################################################
## Perform analysis
###############################################################################
# Set parameters
params <- list(
  sig.file='/g/furlong/project/37_Capture-C/data/diffinter/five_hypotheses_merged/fragments/significant_fragments_fdr_0.1.txt',
  peak.file='/g/furlong/project/37_Capture-C/data/diffinter/five_hypotheses_merged/peaks/significant_peaks_fdr_0.1_0.05.txt',
  diff.hyp='greaterAbs.0',
  same.hyp='lessAbs.1',
  normalisation='distance_size',
  max.padj=0.1,
  sig.padj=0.05,
  window=1000
)
# Generate peaks
peaks <- generate.peaks.from.file(
  sig.file=params$sig.file,
  diff.hyp=params$diff.hyp,
  same.hyp=params$same.hyp,
  normalisation=params$normalisation,
  max.padj=params$max.padj,
  sig.padj=params$sig.padj,
  window=params$window
)
# Annotate peaks with overlaps and write to file
peaks <- annotate.overlaps(peaks)
# Write parameters to file
header <- character(0)
for (p in names(params)) {
  line = paste0('# ', p, ': ', params[p])
  header <- c(header, line)
}
header <- paste0(header, collapse='\n')
write(header, params$peak.file)
# Write peaks to file
write.table(
  peaks, params$peak.file, col.names=T, row.names=F, quote=F, sep='\t',
  append=T
)