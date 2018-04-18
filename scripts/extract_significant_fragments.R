rm(list=ls())
require(purrr)
require(parallel)
require(data.table)

###############################################################################
## Function to extract significant interactions
###############################################################################
find.significant.fragments <- function(
  path, fdr.threshold
) {
  # Read in data
  command <- paste('zcat', path)
  data <- data.table::fread(command, showProgress=F, sep='\t')
  # Find and store significant fragments
  sig <- (!is.na(data$padj)) & (data$padj < 0.1)
  sig.data <- data[sig,]
  # Format and return output
  counts <- c(sum(sig), length(sig))
  output <- list(
    'sig.frags' = sig.data,
    'counts' = counts
  )
  return(output)
}

###############################################################################
## Function to process input files in parallel
###############################################################################
find.significant.fragments.parralel <- function(
  paths, fdr.threshold, cores
) {
  # Read in data
  frag.data <- mclapply(
    paths,
    find.significant.fragments,
    fdr.threshold=fdr.threshold,
    mc.cores=cores
  )
  # Transpose and reorganise data
  frag.data <- purrr::transpose(frag.data)
  frag.data$sig.frags <- data.table::rbindlist(frag.data$sig.frags)
  counts <- do.call(rbind, frag.data$counts)
  frag.data$counts <- data.frame(
    'path' = paths,
    'significant' = counts[,1],
    'total' = counts[,2]
  )
  return(frag.data)
}

###############################################################################
## Perform analysis
###############################################################################
# Set parametrs
params <- list(
  deseqDir='/g/furlong/project/37_Capture-C/data/diffinter/five_hypotheses_merged/fragments',
  deseqPattern='deseq_results.txt.gz$',
  sigFragFile='/g/furlong/project/37_Capture-C/data/diffinter/five_hypotheses_merged/fragments/significant_fragments_fdr_0.1.txt',
  fragCountFile='/g/furlong/project/37_Capture-C/data/diffinter/five_hypotheses_merged/fragments/fragment_counts_fdr_0.1.txt',
  fdrThreshold=0.1,
  cores=4
)
# Extract signficant fragments
deseq.files <- list.files(
  params$deseqDir, params$deseqPattern, full.names=T
)
fragment.data <- find.significant.fragments.parralel(
  deseq.files, params$fdrThreshold, cores=params$cores
)
# Write data to file
write.table(
  fragment.data$sig.frags, file=params$sigFragFile, sep='\t', row.names=F,
  col.names=T, quote=F
)
write.table(
  fragment.data$counts, file=params$fragCountFile, sep='\t', row.names=F,
  col.names=T, quote=F
)








