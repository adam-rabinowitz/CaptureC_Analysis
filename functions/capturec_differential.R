require('DESeq2')

###############################################################################
## Generate sample data
###############################################################################
gen.sample.data <- function(x) {
  # Check consistency of conditions
  cond1 <- x$cond1[1]
  if (any(x$cond1 != cond1)) {
    stop('inconsistent condition 1')
  }
  cond2 <- x$cond2[2]
  if (any(x$cond2 != cond2)) {
    stop('inconsistent condition 1')
  }
  # Create data frame
  sample.data <- data.frame(
    row.names = colnames(x)[grepl('.Rep\\d+$', colnames(x))])
  # Add sample names
  sample.data$sample <- factor(
    gsub('-', '_', gsub('.Rep\\d+$', '', row.names(sample.data))),
    levels = gsub('-', '_', c(cond1, cond2))
  )
  # Add replicate information and return
  sample.data$replicate <- gsub(
    '^.*?(Rep\\d+)$', '\\1', row.names(sample.data))
  return(sample.data)
}

###############################################################################
## Function to extract sample counts from data
###############################################################################
gen.frag.counts <- function(x) {
  # Extract replicates and counts
  counts <- x[,grepl('.Rep\\d+$', colnames(x))]
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
    x=x
  )
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
    'baitName' = x$baitName,
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
## Extract results from dds object using hypothesis
###############################################################################
extract.deseq.results <- function(
  alt.hypothesis, dds, x
) {
  # Extract results for supplied hypothesis
  dds.results <- DESeq2::results(
    dds,
    lfcThreshold=alt.hypothesis$lfcThreshold,
    altHypothesis=alt.hypothesis$altHypothesis  
  )
  # Format and return output
  output <- format.deseq.output(x=x, results=dds.results)
  return(output)
}

###############################################################################
## Perform deseq analysis for a single bait
###############################################################################
deseq.analysis.bait <- function(
  bait.data, fits, min.mean, alt.hypotheses
) {
  # Extract counts
  counts <- gen.frag.counts(bait.data)
  passed <- apply(counts, 1, mean) >= min.mean
  counts <- counts[passed,]
  # Extract additional data
  sample.data <- gen.sample.data(bait.data)
  distances <- gen.frag.distances(bait.data)
  distances <- distances[passed]
  # Create DESeq2 object and add metadata
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData=counts,
    colData=sample.data,
    design=~sample
  )
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
  dds <- DESeq2::estimateSizeFactors(dds, normMatrix=normMatrix)
  dds <- DESeq2::DESeq(dds, fitType='local', betaPrior=F)
  results.list <- lapply(
    alt.hypotheses,
    extract.deseq.results,
    dds=dds,
    x=bait.data[passed,]
  )
  # Format output
  return(results.list)
}

###############################################################################
## Perform deseq analysis for all baits
###############################################################################
deseq.analysis.all <- function(
  bait.list, fits, min.mean, alt.hypotheses, cores
) {
  # Perform differential analysis
  bait.results <- parallel::mclapply(
    bait.list,
    deseq.analysis.bait,
    fits=fits,
    min.mean=min.mean,
    alt.hypotheses=alt.hypotheses,
    mc.cores=cores
  )
  # Transpose data and join
  hypothesis.results <- purrr::transpose(bait.results)
  hypothesis.results <- lapply(
    hypothesis.results,
    data.table::rbindlist
  )
  hypothesis.results <- lapply(
    hypothesis.results,
    base::as.data.frame
  )
  # Adjust pvalue and return
  hypothesis.results <- lapply(
    hypothesis.results,
    function(hr) {
      hr$padj <- p.adjust(hr$pvalue, method='fdr')
      return(hr)
    }
  )
  return(hypothesis.results)
}


