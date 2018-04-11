rm(list=ls())
require(data.table)
# Set parametrs
params <- list(
  deseqDir='/g/furlong/project/37_Capture-C/data/diffinter/five_hypotheses_merged/fragments/',
  deseqPattern='deseq_results.txt.gz',
  sigFragFile='/g/furlong/project/37_Capture-C/data/diffinter/five_hypotheses_merged/fragments/significant_fragments.txt',
  outdir='all_replicate_results/fragments',
  chicagoCRMdir='/g/furlong/project/37_Capture-C/analysis/TS_Capture/CRM_all',
  chicagoTSSdir='/g/furlong/project/37_Capture-C/analysis/TS_Capture/TSS_all',
  chicago.pattern='_(Rep\\d+){2,}.Rds',
  fdrThreshold=0.1
)
###############################################################################
## Function to extract significant interactions
###############################################################################
find.significant.fragments <- function(
  file.paths, fdr.threshold
) {
  # Create output variable
  sig.fragments <- data.frame()
  for (path in file.paths) {
    cat(paste(path, '\n'))
    command <- paste('zcat', path)
    data <- fread(command, showProgress=F)
    sig <- (!is.na(data$padj)) & (data$padj < 0.1)
    cat(paste('no. sig:', sum(sig), '\n'))
    sig.data <- data[sig,]
    sig.fragments <- rbind(sig.fragments, sig.data)
  }
  return(sig.fragments)
}


# Read in significant fragments
if (file.exists(params$sigFragFile)) {
  sig.fragments <- read.table(
    params$sigFragFile, header=T, sep='\t'
  )
} else {
  deseq.files <- list.files(
    params$deseqDir, params$deseqPattern, full.names=T
  )
  sig.fragments <- find.significant.fragments(
    deseq.files, params$fdrThreshold  
  )
  write.table(
    sig.fragments, file=params$sigFragFile, sep='\t', row.names=F,
    col.names=T, quote=F
  )
}
deseq.files <- list.files(
  params$deseqDir, params$deseqPattern
)
sig.fragments <- data.frame(
  for (file)
)