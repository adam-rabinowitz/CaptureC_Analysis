require('DESeq2')

###############################################################################
## Generate sample data
###############################################################################
gen.sample.data <- function(x) {
  # Check consistency of conditions
  if (any(x$cond1 != x$cond1[1])) {
    stop('inconsistent condition 1')
  }
  if (any(x$cond2 != x$cond2[1])) {
    stop('inconsistent condition 1')
  }
  # Create data frame
  sample.data <- data.frame(
    'sample' = rep(
      c(as.character(x$cond1[1]), as.character(x$cond2[2])),
      each=2),
    'replicate' = c('Rep1', 'Rep2', 'Rep1', 'Rep2'))
  # modify row names and return data
  row.names(sample.data) <- paste(
    sample.data$sample,
    sample.data$replicate,
    sep='.')
  return(sample.data)
}

###############################################################################
## Function to extract sample counts from data
###############################################################################
gen.frag.counts <- function(x) {
  # Extract counts
  counts <- x[,c('cond1.1', 'cond1.2', 'cond2.1', 'cond2.2')]
  colnames(counts) <- paste0(
    rep(c(as.character(x$cond1[1]), as.character(x$cond2[1])),
        each=2),
    c('.Rep1', '.Rep2', '.Rep1', '.Rep2'))
  # Convert to integer and return
  counts <- as.matrix(counts)
  counts <- apply(counts, 2, as.integer)
  return(counts)
}

###############################################################################
## Function to extract sample counts from data
###############################################################################
gen.frag.distances <- function(x) {
  # Exract distance and return
  dist <- abs(x$fragDist)
  return(dist)
}

###############################################################################
## Predict distance decay values
###############################################################################
predict.monospline <- function(mono.fit, x) {
  # Generate prediction
  predict <- mgcv::Predict.matrix(
    mono.fit$smooth,
    data.frame(x=x))%*%mono.fit$penalties
  return(predict[,1])
}

###############################################################################
## Extract normalisation factors
###############################################################################
calculate.normalisation <- function(mono.fits, x) {
  # Calculate fits
  fitted <- sapply(
    mono.fits,
    predict.monospline,
    x=x)
  colnames(fitted) <- names(mono.fits)
  # Calculate normalisation factors
  geo.mean <- exp(rowMeans(log(fitted)))
  if (any(is.na(geo.mean))) {
    stop('failure to calulcate normalisation')
  }
  factors <- fitted / geo.mean
  # Create and return output
  output <- list(
    'fitted'=fitted,
    'factors'=factors)
  return(output)
}

###############################################################################
## Format output data
###############################################################################
format.deseq.output <- function(x, results) {
  output <- data.frame(
    'cond1' = x$cond1,
    'cond2' = x$cond2,
    'baitID' = x$baitID,
    'baitChr' = x$baitChr,
    'baitStart' = x$baitStart,
    'baitEnd' = x$baitEnd,
    'fragID' = x$fragID,
    'fragChr' = x$fragChr,
    'fragStart' = x$fragStart,
    'fragEnd' = x$fragEnd,
    'fragDist' = x$fragDist,
    'baseMean' = results$baseMean,
    'lfc' = results$log2FoldChange,
    'pvalue' = results$pvalue,
    'padj' = results$padj
  )
  return(output)
}

###############################################################################
## Generate count filter
###############################################################################
perform.deseq.analysis <- function(x, fits, min.sum) {
  # Extract counts
  counts <- gen.frag.counts(x)
  passed <- rowSums(counts) >= min.sum
  counts <- counts[passed,]
  # Extract additional data
  sample.data <- gen.sample.data(x)
  distances <- gen.frag.distances(x)
  distances <- distances[passed]
  # Create DESeq2 object and add metadata
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData=counts,
    colData=sample.data,
    design=~sample)
  S4Vectors::metadata(dds)$distances <- distances
  # Perform distance decay normalisation if fits supplied
  if (!is.null(fits)) {
    normMatrix <- calculate.normalisation(
      fits[colnames(dds)], distances)$factors
  } else {
    normMatrix <- matrix(1, nrow=nrow(counts), ncol=ncol(counts))
    colnames(normMatrix) <- colnames(counts)
  }
  # Perform DESeq2 analysis
  dds <- estimateSizeFactors(dds, normMatrix=normMatrix)
  dds <- DESeq2::DESeq(dds, fitType='local', betaPrior=T)
  results <- DESeq2::results(dds)
  # Format output
  output <- format.deseq.output(
    x=x[passed,],
    results=results)
  return(output)
}