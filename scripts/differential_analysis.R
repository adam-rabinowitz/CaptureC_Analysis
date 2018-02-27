# Script to generate distance-decay plots for replicates of all samples
rm(list=ls())
require('parallel')
source('~/github/CaptureC_Analysis/functions/capturec_normalisation.R')
source('~/github/CaptureC_Analysis/functions/capturec_input.R')
source('~/github/CaptureC_Analysis/functions/capturec_differential.R')
# Set and print parameters
params <- list(
  'indir' = '~/differential/newCounts/',
  'pattern' = '_Rep\\dRep\\d.counts.txt$',
  'outdir' = '~/differential/deseq2Results/',
  'comparisons' = '~/differential/comparisons.txt',
  'chrfile' = '~/differential/chr_sizes.txt',
  'mindist' = 2000,
  'binsize' = 1000,
  'k' = 20,
  'cores' = 6)
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
  # Perform differential analysis
  deseq.results <- mclapply(
    merged.data,
    perform.deseq.analysis,
    fits=frequency.fits,
    mc.cores=params$cores)
  # Merge data and recalculate p-value
  deseq.results <- rbindlist(deseq.results)
  deseq.results$padj <- p.adjust(deseq.results$pvalue, method='fdr')
  # Save file
  out.path <- file.path(params$outdir, paste0(cmp.name, '.deseq2_results.txt'))
  write.table(
    deseq.results, out.path, row.names=F, col.names=T, sep='\t', quote=F)
}

