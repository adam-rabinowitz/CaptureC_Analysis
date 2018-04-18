require('DESeq2')
require('parallel')

###############################################################################
## Generate sample data
###############################################################################
gen.sample.data <- function(x) {
  # Create data frame
  sample.data <- data.frame(
    row.names = colnames(x)[grepl('.Rep\\d+$', colnames(x))]
  )
  # Add dataset
  sample.data$dataset <- gsub(
    '^(CRM|TSS)_(\\d+)-(\\d+h)_([[:alnum:]]+)\\.(Rep\\d+)$',
    '\\1',
    row.names(sample.data)
  )
  # Add time
  sample.data$time <- gsub(
    '^(CRM|TSS)_(\\d+)-(\\d+h)_([[:alnum:]]+)\\.(Rep\\d+)$',
    paste('\\2', '\\3', sep='.'),
    row.names(sample.data)
  )
  # Add tissue
  sample.data$tissue <- gsub(
    '^(CRM|TSS)_(\\d+)-(\\d+h)_([[:alnum:]]+)\\.(Rep\\d+)$',
    '\\4',
    row.names(sample.data)
  )
  # Add condition
  sample.data$condition <- paste(
    sample.data$time, sample.data$tissue, sep='_'
  )
  sample.data$condition <- factor(
    sample.data$condition,
    levels=unique(sample.data$condition)
  )
  # Add replicate
  sample.data$replicate <- gsub(
    '^(CRM|TSS)_(\\d+)-(\\d+h)_([[:alnum:]]+)\\.(Rep\\d+)$',
    '\\5',
    row.names(sample.data)
  )
  # Return data
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
  rownames(counts) <- x$fragID
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
## Create DESEqDataSet object for a single bait
###############################################################################
create.dds.bait <- function(bait.data, min.mean) {
  # Extract counts
  counts <- gen.frag.counts(bait.data)
  passed <- apply(counts, 1, mean) >= min.mean
  counts <- counts[passed,]
  # Extract additional data
  sample.data <- gen.sample.data(bait.data)
  distances <- gen.frag.distances(bait.data)
  distances <- distances[passed]
  # Create DESeq2 object
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData=counts,
    colData=sample.data,
    design=~condition
  )
  # Add metadata
  row.names(bait.data) <- bait.data$fragID
  S4Vectors::metadata(dds)$baitdata <- bait.data[
    passed, !grepl('\\.Rep\\d+$', colnames(bait.data))]
  # Return dds object
  return(dds)
}

###############################################################################
## Create DESEqDataSet object for a single bait
###############################################################################
create.dds.all <- function(bait.list, min.mean, cores=cores) {
  # Create dds from list of bait data and return
  dds.list <- parallel::mclapply(
    bait.list,
    create.dds.bait,
    min.mean=min.mean,
    mc.cores=cores
  )
  return(dds.list)
}

###############################################################################
## Perform deseq analysis for a single DESeqDataSet
###############################################################################
deseq.analysis.bait <- function(
  dds, fits
) {
  # Extract distances
  distances <- abs(metadata(dds)$baitdata$fragDist)
  # Perform distance decay normalisation if fits supplied
  if (!is.null(fits)) {
    normMatrix <- calculate.normalisation(
      fits[colnames(dds)], distances)$factors
    dds <- DESeq2::estimateSizeFactors(dds, normMatrix=normMatrix)
  }
  # Perform DESeq2 analysis and return dds object
  dds <- DESeq2::DESeq(dds, fitType='local', betaPrior=F)
  return(dds)
}

###############################################################################
## Perform deseq analysis for all baits
###############################################################################
deseq.analysis.all <- function(
  dds.list, fits, cores
) {
  # Perform differential analysis and return results
  dds.list <- parallel::mclapply(
    dds.list,
    deseq.analysis.bait,
    fits=fits,
    mc.cores=cores
  )
  return(dds.list)
}

###############################################################################
## Extract results from single dds object using single hypothesis and contrast
###############################################################################
deseq.results.bait <- function(
  dds, dataset, normalisation, alt.hypothesis, contrast, alpha
) {
  # Edit contrasts and extract count for replicates
  alt.contrast <- gsub('-', '.', contrast)
  replicates <- colnames(dds)[colData(dds)$condition %in% alt.contrast]
  replicate.counts <- counts(dds, normalized=T)[,replicates]
  replicate.means <- apply(replicate.counts, 1, mean)
  # Extract results for supplied hypothesis
  dds.results <- DESeq2::results(
    dds,
    contrast = c('condition', alt.contrast[1], alt.contrast[2]),
    lfcThreshold=alt.hypothesis$lfcThreshold,
    altHypothesis=alt.hypothesis$altHypothesis,
    independentFiltering=T,
    filter=replicate.means, ############# Required to remove low counts 
    alpha=alpha
  )
  # Extract bait data and check identity
  bait.data <- metadata(dds)$baitdata
  if (!identical(rownames(bait.data), row.names(dds.results))) {
    stop('result order differs to bait')
  }
  # Create output and return
  output <- data.frame(
    'dataset' = dataset,
    'normalisation' = normalisation,
    'cond1' = contrast[1],
    'cond2' = contrast[2],
    'altHyp' = paste(
      alt.hypothesis$altHypothesis,
      alt.hypothesis$lfcThreshold,
      sep='_'
    ),
    'baitName' = bait.data$baitName,
    'baitID' = bait.data$baitID,
    'baitChr' = bait.data$baitChr,
    'baitStart' = bait.data$baitStart,
    'baitEnd' = bait.data$baitEnd,
    'fragID' = bait.data$fragID,
    'fragChr' = bait.data$fragChr,
    'fragStart' = bait.data$fragStart,
    'fragEnd' = bait.data$fragEnd,
    'fragDist' = bait.data$fragDist,
    'baseMean' = dds.results$baseMean,
    'repMean' = replicate.means,
    'lfc' = dds.results$log2FoldChange,
    'pvalue' = dds.results$pvalue,
    'padj' = dds.results$padj
  )
  return(output)
}

###############################################################################
## Extract results from many dds objects using many hypotheses and contrasts
###############################################################################
deseq.results.all <- function(
  dds.list, contrasts, alt.hypotheses, alpha, dataset, normalisation, outdir,
  suffix, cores
) {
  # Loop through contrasts and hypotheses to generate results
  combinations <- expand.grid(names(contrasts), names(alt.hypotheses))
  for (i in 1:nrow(combinations)) {
    # Extract names
    contrast.name <- as.character(combinations[i, 1])
    alt.hypothesis.name <- as.character(combinations[i, 2])
    out.file <- paste(
      dataset, contrast.name, alt.hypothesis.name, suffix, sep='.'
    )
    out.path <- file.path(outdir, out.file)
    # Extract results
    results <- parallel::mclapply(
      dds.list,
      deseq.results.bait,
      dataset=dataset,
      normalisation=normalisation,
      alt.hypothesis=alt.hypotheses[[alt.hypothesis.name]],
      contrast=contrasts[[contrast.name]],
      alpha=alpha,
      mc.cores=cores
    )
    results <- data.table::rbindlist(results)
    # Write to file
    data.table::fwrite(
      results, out.path, sep='\t', quote=F, col.names=T, row.names=F,
      showProgress=F
    )
    # Gzip file
    command = paste('gzip', out.path)
    system(command, wait=F)
    # Clean up after loop
    rm(results)
    gc(verbose=F)
  }
}

