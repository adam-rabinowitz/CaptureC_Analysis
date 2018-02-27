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
  # Calculate groups split by 2 * window
  groups <- cumsum(ndistance > window)
  groups <- c(1, groups + 1)
  groups <- factor(groups, levels=1:tail(groups, 1))
  # Split and return input data
  x.list <- split(x, groups)
  return(x.list)
}

###############################################################################
## Function to split peaks
###############################################################################
running.lfc <- function(lfc, padj, fragments, weights) {
  # Check all fragments have weights
  if (length(weights) != fragments) {
    stop('length of weights does not equal fragments')}
  # Check no. frags is odd
  if ((fragments %% 2) != 1) {
    stop('fragments must be odd')}
  # Extract data
  output <- rep(NA, length(lfc))
  overhang <- (fragments - 1) / 2
  lfc <- c(rep(0, overhang), lfc, rep(0, overhang))
  padj <- c(rep(1, overhang), padj, rep(1, overhang))
  log10.padj <- -1 * (log10(padj))
  # Perform running lfc calculation
  last.start <- (length(lfc) - fragments) + 1
  for (i in 1:last.start) {
    indices <- i:(i + (fragments - 1))
    adj.weights <- log10.padj[indices] * weights
    mean.lfc <- weighted.mean(lfc[indices], adj.weights)
    output[i] <- mean.lfc
  }  
  # Return data
  return(output)
}

###############################################################################
## function to trim peaks
###############################################################################
trim.peaks <- function(x, max.padj, sign) {
  # Extract which values to trim
  if (sign == 'pos') {
    ftrim <- cumsum(x$lfc > 0 & x$padj <= max.padj) == 0
    rtrim <- rev(cumsum(rev(x$lfc) > 0 & rev(x$padj) <= max.padj)) == 0
  } else if (sign == 'neg') {
    ftrim <- cumsum(x$lfc < 0 & x$padj <= max.padj) == 0
    rtrim <- rev(cumsum(rev(x$lfc) < 0 & rev(x$padj) <= max.padj)) == 0
  } else {
    stop('uknown sign argument')
  }
  # Trim values and return
  x.filt <- x[!(ftrim | rtrim),]
  if (nrow(x.filt) == 0) {
    x.filt <- NA
  }
  return(x.filt)
}

###############################################################################
## Split on lfc changes
###############################################################################
split.peaks <- function(x, sig.padj, max.padj, sign) {
  # Calculate breaks
  if (sign == 'pos') {
    breaks <- sign(x$lfc) < 0 & x$padj <= max.padj
  } else if (sign == 'neg') {
    breaks <- sign(x$lfc) > 0 & x$padj <= max.padj
  } else {
    stop('unknown sign argument')
  }
  # Return data if no breaks found
  if (!any(breaks)) {
    return(list('1'=x))
  }
  # Remove break values
  print('splitting a peak')
  x.filt <- x[!breaks,]
  x.list <- split(x.filt, cumsum(breaks)[!breaks])
  # Trim peaks
  x.trim <- lapply(x.list, trim.peaks, max.padj=max.padj, sign=sign)
  x.trim <- x.trim[sapply(x.trim, is.data.frame)]
  # Remove non-significant peaks
  x.filter <- x.trim[sapply(x.trim, function(z) {any(z$padj <= sig.padj)})]
  names(x.filter) <- 1:length(x.filter)
  return(x.filter)
}

###############################################################################
## Find significant peaks
###############################################################################
find.peaks <- function(
  x, sig.padj, max.padj, window
) {
  # Rename row and filter data
  x.filt <- x[!is.na(x$padj) & x$padj <= max.padj,]
  if (nrow(x.filt) == 0) {
    return(list())}
  # Split data into positive and negative lfc
  sign.factors <- factor(
    sign(x.filt$lfc),
    levels=c(-1, 1),
    labels=c('neg', 'pos'))
  x.filt.list <- split(x.filt, sign.factors)
  x.filt.list <- x.filt.list[sapply(x.filt.list, nrow) > 0]
  # split into windows
  windows <- lapply(x.filt.list, split.window, window=window)
  # Find peaks for positive and negative lfc sequentially
  output.peaks <- list()
  for (sign in c('pos', 'neg')) {
    # Find signficant peaks for strand
    sign.windows <- windows[[sign]]
    sign.peaks <- sign.windows[
      sapply(sign.windows, function(z) {min(z$padj)}) <= sig.padj]
    # Generate empty list in absence of signifcant peaks
    if (length(sign.peaks) == 0) {
      sign.split.peaks <- list()
    # Else process signficant peaks
    } else {
      # Extract peaks for original data
      sign.full.peaks <- lapply(
        sign.peaks,
        function(z) {
          indices <- (
            x$fragID >= min(z$fragID) &
            x$fragID <= max(z$fragID))
          x[indices,]})
      # Split peaks
      sign.split.peaks <- lapply(
        sign.full.peaks,
        split.peaks,
        sig.padj=sig.padj,
        max.padj=max.padj,
        sign=sign)
      sign.split.peaks <- do.call(c, sign.split.peaks)
    }
    # Store data
    output.peaks[[sign]] <- sign.split.peaks
  }
  # Flatten, sort and return data
  output.peaks <- do.call(c, output.peaks)
  peak.order <- order(sapply(output.peaks, function(z) {min(z$fragStart)}))
  output.peaks <- output.peaks[peak.order]
  return(output.peaks)
}

###############################################################################
## Extract peaks from results
###############################################################################
find.results.peaks <- function(
  results, sig.padj, max.padj, window
) {
  # Extrcat peaks from results and flatten
  peaks <- lapply(
    results,
    find.peaks,
    sig.padj=sig.padj,
    max.padj=max.padj,
    window=window)
  peaks <- do.call(c, peaks)
  # Extract summary statistics from peaks
  peak.summary <- lapply(
    peaks,
    function(z) {
      # Calculate peak metrics
      min.index <- which.min(z$padj)
      dist.sign <- sign(z$fragDist)
      if (all(dist.sign == 1)) {
        peak.dist <- min(z$fragDist)
      } else if (all(dist.sign == -1)) {
        peak.dist <- max(z$fragDist)
      } else {
        print(z)
        stop('unexpected distance sign')
      }
      # Create peak dataframe
      peak.df <- data.frame(
        'peakID' = NA,
        'cond1' = z$cond1[1],
        'cond2' = z$cond2[1],
        'baitID' = z$baitID[1],
        'baitChr' = z$baitChr[1],
        'baitStart' = z$baitStart[1],
        'baitEnd' = z$baitEnd[1],
        'peakChr' = z$fragChr[1],
        'peakStart' = min(z$fragStart),
        'peakEnd' = max(z$fragEnd),
        'peakDist' = min(z$fragDist),
        'peakFragments' = nrow(z),
        'peakSigFragments' = sum(z$padj <= sig.padj),
        'peakWeightedLFC' = weighted.mean(
          z$lfc,
          -1 * log10(z$padj)),
        'peakMinStart' = z$fragStart[min.index],
        'peakMinEnd' = z$fragEnd[min.index],
        'peakMinLFC' = z$lfc[min.index],
        'peakMinPadj' = z$padj[min.index]
      )
      # Check data and return
      if (peak.df$peakFragments < 1) {
        print(peak.df)
        stop('peaks has no signficant fragments!')}
      if (peak.df$peakMinPadj > sig.padj) {
        print(peak.df)
        stop('peak has high padj!')}
      return(peak.df)})
  # Merge peaks and add peak ID
  all.peaks <- do.call(rbind, peak.summary)
  all.peaks <- all.peaks[order(
    all.peaks$baitChr, all.peaks$baitStart, all.peaks$peakChr,
    all.peaks$peakStart),]
  peak.no <- unlist(
    sapply(rle(as.character(all.peaks$baitID))$lengths, seq))
  # Create amd store peak ID
  all.peaks$peakID <- paste(
    substring(all.peaks$cond1, 1, 3), substring(all.peaks$cond1, 5),
    substring(all.peaks$cond2, 5), all.peaks$baitID, peak.no, sep='.')
  # Return data
  return(all.peaks)
}