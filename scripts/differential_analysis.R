# Script to generate distance-decay plots for replicates of all samples
rm(list=ls())
require('parallel')
source('/g/furlong/adamr/github/CaptureC_Analysis/functions/capturec_normalisation.R')
source('/g/furlong/adamr/github/CaptureC_Analysis/functions/capturec_input.R')
source('/g/furlong/adamr/github/CaptureC_Analysis/functions/capturec_differential.R')
# Set and print parameters
params <- list(
  'indir' = '/g/furlong/project/37_Capture-C/data/diffinter/input/',
  'pattern' = '_Rep\\dRep\\d.counts.txt$',
  'outdir' = '/g/furlong/project/37_Capture-C/data/diffinter/results/',
  'comparisons' = '/g/furlong/project/37_Capture-C/data/diffinter/comparisons.txt',
  'chrfile' = '/g/furlong/project/37_Capture-C/data/diffinter/chr_sizes.txt',
  'mindist' = 2000,
  'binsize' = 1000,
  'k' = 20,
  'minsum' = 4,
  'cores' = 8)
cat('Parameters:\n')
for (p in names(params)) {
  cat(paste0('\t', p, ' : ', params[[p]], '\n'))}
# Extract input file paths and print
input.paths <- list.files(
  params$indir,
  pattern=params$pattern,
  full.names=T)
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
comparisons <- read.table(
  params$comparisons, sep='\t', stringsAsFactors=F)
# check comparisons and get names
if (
  any(
    substring(comparisons[,1], 1, 3) !=
    substring(comparisons[,2], 1, 3))
) {
  stop('unmatched pulldowns')}
cmp.names <- paste(
  substring(comparisons[,1], 1, 3),
  substring(comparisons[,1], 5),
  substring(comparisons[,2], 5),
  sep='.')
# Split comparsions
comparisons <- split(comparisons, cmp.names)
comparisons <- lapply(comparisons, as.character)
cat('\nComparisons:\n')
for (cn in names(comparisons)) {
  cat(paste0('\t', cn, ':\n', '\t\t - ', comparisons[[cn]][1], '\n\t\t - ',
    comparisons[[cn]][2], '\n'))}
# Create fits
frequency.fits <- lapply(
  input.paths,
  distance.decay.fit,
  min.dist=params$mindist,
  bin.size=params$binsize,
  k=params$k,
  chr.sizes=chr.sizes,
  cores=params$cores)
frequency.fits <- do.call(c, frequency.fits)
# Extract intrachromosomal count data
input.counts <- lapply(
  input.paths,
  read.split.replicates,
  min.dist=params$mindist,
  intra.only=T)
names(input.counts) <- gsub('_Rep1Rep2.counts.txt', '', basename(input.paths))
# Loop through comparisons
cat('\nStarted:\n')
for (cmp.name in names(comparisons)) {
  # Extract data
  cat(paste0('\t- ', cmp.name, '\n'))
  cond1 <- comparisons[[cmp.name]][1]
  cond2 <- comparisons[[cmp.name]][2]
  # Merge data
  merged.data <- merge.datasets(
    input.counts[[cond1]],
    input.counts[[cond2]],
    cores=params$cores)
  # Perform differential analysis with fit
  deseq.norm.results <- mclapply(
    merged.data,
    perform.deseq.analysis,
    fits=frequency.fits,
    min.sum=params$minsum,
    mc.cores=params$cores)
  # Merge data and recalculate p-value
  deseq.norm.results <- rbindlist(deseq.norm.results)
  deseq.norm.results$padj <- p.adjust(deseq.norm.results$pvalue, method='fdr')
  # Save file
  norm.path <- file.path(
    params$outdir,
    paste0(cmp.name, '.deseq2_norm_results.txt'))
  write.table(
    deseq.norm.results, norm.path, row.names=F, col.names=T, sep='\t', quote=F)
  # Perform differential analysis without fit
  deseq.raw.results <- mclapply(
    merged.data,
    perform.deseq.analysis,
    fits=NULL,
    min.sum=params$minsum,
    mc.cores=params$cores)
  # Merge data and recalculate p-value
  deseq.raw.results <- rbindlist(deseq.raw.results)
  deseq.raw.results$padj <- p.adjust(deseq.raw.results$pvalue, method='fdr')
  # Save file
  raw.path <- file.path(
    params$outdir,
    paste0(cmp.name, '.deseq2_raw_results.txt'))
  write.table(
    deseq.raw.results, raw.path, row.names=F, col.names=T, sep='\t', quote=F)
}

