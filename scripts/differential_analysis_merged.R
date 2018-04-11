# Script to generate distance-decay plots for replicates of all samples
rm(list=ls())
source('/g/furlong/adamr/github/CaptureC_Analysis/functions/capturec_normalisation.R')
source('/g/furlong/adamr/github/CaptureC_Analysis/functions/capturec_input.R')
source('/g/furlong/adamr/github/CaptureC_Analysis/functions/capturec_differential.R')
# Set and print parameters
params <- list(
  'indir' = '/g/furlong/project/37_Capture-C/data/diffinter/input',
  'pattern' = '_(Rep\\d+){2,}.counts.txt$',
  'outdir' = '/g/furlong/project/37_Capture-C/data/diffinter/five_hypotheses_merged/fragments',
  'normsuffix' = 'norm_deseq_results.txt',
  'rawsuffix' = 'raw_deseq_results.txt',
  'chrfile' = '/g/furlong/project/37_Capture-C/data/diffinter/chr_sizes.txt',
  'mindist' = 2000,
  'binsize' = 1000,
  'k' = 20,
  'minmean' = 0,
  'cores' = 8,
  'fdrcutoff' = 0.1
)
cat('Parameters:\n')
for (p in names(params)) {
  cat(paste0('\t', p, ' : ', params[[p]], '\n'))}
rm(p)
# Set and print alternative hypotheses for deseq
alt.hypotheses <- list(
  'above0.0' = list(
    'lfcThreshold' = 0,
    'altHypothesis' = 'greaterAbs'
  ),
  'above0.5' = list(
    'lfcThreshold' = 0.5,
    'altHypothesis' = 'greaterAbs'
  ),
  'above1.0' = list(
    'lfcThreshold' = 1,
    'altHypothesis' = 'greaterAbs'
  ),
  'below0.5' = list(
    'lfcThreshold' = 0.5,
    'altHypothesis' = 'lessAbs'
  ),
  'below1.0' = list(
    'lfcThreshold' = 1,
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
rm(a)
# Extract input file paths and print
input.paths <- list(
  'CRM' = list.files(
    params$indir,
    pattern=paste0('CRM.*?', params$pattern),
    full.names=T
  ),
  'TSS' = list.files(
    params$indir,
    pattern=paste0('TSS.*?', params$pattern),
    full.names=T
  )
)
cat('\nInput Paths:\n')
for (dataset in names(input.paths)) {
  cat(paste0('\t', dataset, ':\n'))
  for (path in input.paths[dataset]) {
    cat(paste0('\t\t - ', path, '\n'))
  }
}
rm(dataset, path)
# Set contrasts for analysis and print
contrasts <- list(
  '6-8h_DN.6-8h_elav' = c('6-8h_DN', '6-8h_elav'),
  '6-8h_DN.6-8h_Mef2' = c('6-8h_DN', '6-8h_Mef2'),
  '6-8h_elav.6-8h_Mef2' = c('6-8h_elav', '6-8h_Mef2'),
  '10-12h_DN.10-12h_elav' = c('10-12h_DN', '10-12h_elav'),
  '10-12h_DN.10-12h_Mef2' = c('10-12h_DN', '10-12h_Mef2'),
  '10-12h_elav.10-12h_Mef2' = c('10-12h_elav', '10-12h_Mef2'),
  '6-8h_DN.10-12h_DN' = c('6-8h_DN', '10-12h_DN'),
  '6-8h_elav.10-12h_elav' = c('6-8h_elav', '10-12h_elav'),
  '6-8h_Mef2.10-12h_Mef2' = c('6-8h_Mef2', '10-12h_Mef2')
)
cat('\nContrasts:\n')
for (contrast in names(contrasts)) {
  cat(paste0('\t', contrast, ':\n'))
  for (condition in contrasts[contrast]) {
    cat(paste0('\t\t - ', condition, '\n'))
  }
}
rm(contrast, condition)
# Extract chromosome sizes and print
chr.sizes <- read.table(
  params$chrfile, sep='\t', col.names=c('name', 'length'))
chr.sizes <- split(chr.sizes$length, chr.sizes$name)
cat('\nChromosome Sizes:\n')
for (chr in names(chr.sizes)) {
  cat(paste0('\t', chr, ' : ', chr.sizes[chr], '\n'))}
rm(chr)
# Loop through 
for (dataset in names(input.paths)) {
  # Create fits
  distance.fits <- distance.decay.fit.all(
    input.paths[[dataset]],
    min.dist=params$mindist,
    bin.size=params$binsize,
    k=params$k,
    chr.sizes=chr.sizes,
    cores=params$cores
  )
  # Read in counts for each condition
  merged.counts <- extract.merged.counts(
    input.paths[[dataset]],
    suffix=params$pattern,
    min.dist=params$mindist,
    intra.only=T,
    cores=params$cores
  )
  # Create list of DSESeqDataSets and delete counts
  dds.list <- create.dds.all(
    merged.counts,
    min.mean=params$minmean,
    cores=params$cores
  )
  rm(merged.counts)
  # Perform normalised analysis
  norm.dds <- deseq.analysis.all(
    dds.list=dds.list,
    fits=distance.fits,
    cores=params$cores
  )
  # Save normalised results and delete
  deseq.results.all(
    dds.list=norm.dds,
    contrasts=contrasts,
    alt.hypotheses=alt.hypotheses,
    alpha=params$fdrcutoff,
    dataset=dataset,
    normalisation='distance_size',
    outdir=params$outdir,
    suffix=params$normsuffix,
    cores=params$cores
  )
  rm(norm.dds)
  # Perform raw analysis
  raw.dds <- deseq.analysis.all(
    dds.list=dds.list,
    fits=NULL,
    cores=params$cores
  )
  rm(dds.list)
  # Save raw results and delete data
  deseq.results.all(
    dds.list=raw.dds,
    contrasts=contrasts,
    alt.hypotheses=alt.hypotheses,
    alpha=params$fdrcutoff,
    dataset=dataset,
    normalisation='size',
    outdir=params$outdir,
    suffix=params$rawsuffix,
    cores=params$cores
  )
  rm(raw.dds)
}