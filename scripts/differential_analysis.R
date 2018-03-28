# Script to generate distance-decay plots for replicates of all samples
rm(list=ls())
require('parallel')
source('~/github/CaptureC_Analysis/functions/capturec_normalisation.R')
source('~/github/CaptureC_Analysis/functions/capturec_input.R')
source('~/github/CaptureC_Analysis/functions/capturec_differential.R')
# Set and print parameters
params <- list(
  'indir' = '/g/furlong/project/37_Capture-C/data/diffinter/input',
  'pattern' = '_Rep\\d+Rep\\d+.counts.txt$',
  'outdir' = '/g/furlong/project/37_Capture-C/data/diffinter/all_replicate_results/fragments',
  'comparisons' = '/g/furlong/project/37_Capture-C/data/diffinter/comparisons.txt',
  'chrfile' = '/g/furlong/project/37_Capture-C/data/diffinter/chr_sizes.txt',
  'mindist' = 2000,
  'binsize' = 1000,
  'k' = 20,
  'minmean' = 1,
  'cores' = 8
)
cat('Parameters:\n')
for (p in names(params)) {
  cat(paste0('\t', p, ' : ', params[[p]], '\n'))}
# Set and print alternative hypotheses for deseq
alt.hypotheses <- list(
  'above0.0' = list(
    'lfcThreshold' = 0,
    'altHypothesis' = 'greaterAbs'
  ),
  'above1.0' = list(
    'lfcThreshold' = 1,
    'altHypothesis' = 'greaterAbs'
  ),
  'below0.5' = list(
    'lfcThreshold' = 0.5,
    'altHypothesis' = 'lessAbs'
  )
)
cat('\nAltHypotheses:\n')
for (a in names(alt.hypotheses)) {
  cat(
    paste0(
      '\t', a, ':\n\t\tlfcThreshold : ', alt.hypotheses[[a]]$lfcThreshold,
      '\n\t\taltHypothesis : ', alt.hypotheses[[a]]$altHypothesis, '\n'
    )
  )
}
# Extract input file paths and print
input.paths <- list.files(
  params$indir,
  pattern=params$pattern,
  full.names=T
)
cat('\nInput Paths:\n')
for (ip in input.paths) {
  cat(paste0('\t - ', ip, '\n'))}
# Extract chromosome sizes and print
chr.sizes <- read.table(
  params$chrfile, sep='\t', col.names=c('name', 'length'))
chr.sizes <- split(chr.sizes$length, chr.sizes$name)
cat('\nChromosome Sizes:\n')
for (chr in names(chr.sizes)) {
  cat(paste0('\t', chr, ' : ', chr.sizes[chr], '\n'))}
# Extract comparisons
comparisons <- extract.comparisons(params$comparisons)
cat('\nComparisons:\n')
for (cn in names(comparisons)) {
  cat(
    paste0(
      '\t', cn, ':\n', '\t\t - ', comparisons[[cn]][1], '\n\t\t - ',
      comparisons[[cn]][2], '\n'
    )
  )
}
# Create fits
distance.fits <- distance.decay.fit.all(
  input.paths,
  min.dist=params$mindist,
  bin.size=params$binsize,
  k=params$k,
  chr.sizes=chr.sizes
)
# Read in counts for each condition
condition.counts <- extract.condition.counts(
  input.paths, suffix=params$pattern, min.dist=params$mindist, intra.only=T,
  cores=params$cores
)
# Loop through comparisons
cat('\nStarted:\n')
for (cmp.name in names(comparisons)) {
  # Extract data
  cat(paste0('\t- ', cmp.name, '\n'))
  cond1 <- comparisons[[cmp.name]][1]
  cond2 <- comparisons[[cmp.name]][2]
  # Merge data
  merged.data <- merge.conditions.all(
    condition.counts[c(cond1, cond2)],
    cores=params$cores
  )
  # Perform deseq analysis
  norm.results <- deseq.analysis.all(
    bait.list=merged.data,
    fits=distance.fits,
    min.mean=1,
    alt.hypotheses=alt.hypotheses,
    cores=1
  )
  # Save file
  for (ah in names(alt.hypotheses)) {
    # Create output path
    norm.path <- file.path(
      params$outdir,
      paste(cmp.name, ah, 'deseq2_norm_results.txt', sep='.')
    )
    # Write data to file
    write.table(
      norm.results[[ah]], norm.path, row.names=F, col.names=T, sep='\t', quote=F 
    )
  }
  # Merge data and recalculate p-value
  raw.results <- deseq.analysis.all(
    bait.list=merged.data,
    fits=NULL,
    min.mean=params$minmean,
    alt.hypotheses=alt.hypotheses,
    cores=1
  )
  # Save raw results for each hypothesis
  for (ah in names(alt.hypotheses)) {
    # Create output path
    raw.path <- file.path(
      params$outdir,
      paste(cmp.name, ah, 'deseq2_raw_results.txt', sep='.')
    )
    # Write data to file
    write.table(
      raw.results[[ah]], raw.path, row.names=F, col.names=T, sep='\t', quote=F 
    )
  }
}

