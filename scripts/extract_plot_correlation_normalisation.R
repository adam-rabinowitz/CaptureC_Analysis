# Script to generate distance-decay plots for replicates of all samples
rm(list=ls())
require(ggplot2)
require(reshape2)
source('~/github/CaptureC_Analysis/normalisation/capturec_normalisation.R')
source('~/github/CaptureC_Analysis/normalisation/capturec_input.R')
# Find input files
indir <- '~/differential/newCounts/'
outdir <- '~/differential/factorPlots/'
comp.file <- '~/differential/comparisons.txt'
input.paths <- list.files(
  '~/differential/newCounts',
  pattern='_Rep\\dRep\\d.counts.txt',
  full.names=T)
# Set chromosome sizes
chr.sizes <- list(
  'chr2L'=23513712,
  'chr2R'=25286936,
  'chr3L'=28110227,
  'chr3R'=32079331,
  'chr4'=1348131,
  'chrX'=23542271)
# Create fits
frequency.fits <- lapply(
  input.paths,
  distance.decay.fit,
  min.dist=2000,
  bin.size=1000,
  k=20,
  chr.sizes=chr.sizes,
  cores=6)
frequency.fits <- do.call(c, frequency.fits)
# Read input data
input.data <- lapply(
  input.paths,
  function(z) {
    data <- fread(z)
    data.filt <- data[data$baitChr == data$fragChr,]
    data.list <- split(data.filt, data.filt$baitID)
    return(data.list)})
names(input.data) <- gsub('_Rep1Rep2.counts.txt', '', basename(input.paths))
# Loop through comparisons and merge counts
comparisons <- read.table(comp.file, sep='\t', stringsAsFactors=F)
common.counts <- list()
for (cmp.no in 1:nrow(comparisons)) {
  # Extract data
  cond1 <- comparisons[cmp.no, 1]
  cond2 <- comparisons[cmp.no, 2]
  # Extract comparison name
  if (substring(cond1, 1, 3) != substring(cond2, 1, 3)) {
    stop('datasets from different series')}
  cmp.name <- paste(
    substring(cond1, 1, 3),
    substring(cond1, 5),
    substring(cond2, 5),
    sep='.')
  print(cmp.name)
  # Merged data by baits
  all.baits <- union(names(input.data[[cond1]]), names(input.data[[cond2]]))
  merged.baits <- mclapply(
    all.baits,
    function(bait) {
      merge.bait.data(
        input.data[[cond1]][[bait]],
        input.data[[cond2]][[bait]])},
    mc.cores=6)
  # Remove all fragments with zero counts in a single condition
  count.cols <- c('cond1.1', 'cond1.2', 'cond2.1', 'cond2.2')
  filt.baits <- lapply(
    merged.baits,
    function(bait) {
      bait[rowSums(bait[,count.cols] > 0) == 4,]})
  # Store data
  names(filt.baits) <- all.baits
  common.counts[[cmp.name]] <- filt.baits
}

###############################################################################
## Collect correlation stats
###############################################################################
extract.correlation.metrics <- function(counts) {
  # Perform correlation
  log2.counts <- log2(counts)
  log2.cor <- cor(log2.counts)
  # Extract correlation data
  cor.df <- data.frame(
    'class' = c(
      'intra.replicate',
      'intra.replicate',
      'inter.replicate',
      'inter.replicate',
      'inter.replicate',
      'inter.replicate'),
    'correlation' = c(
      log2.cor['cond1.1', 'cond1.2'],
      log2.cor['cond2.1', 'cond2.2'],
      log2.cor['cond1.1', 'cond2.1'],
      log2.cor['cond1.1', 'cond2.2'],
      log2.cor['cond1.2', 'cond2.1'],
      log2.cor['cond1.2', 'cond2.2']
    ))
  return(cor.df)
}

###############################################################################
## Extract correlations
###############################################################################
calculate.normalisation.correlation <- function(bait.data) {
  # Extract normalisation factors
  replicates <- paste0(
    rep(c(bait.data$cond1[1], bait.data$cond2[1]), each=2),
    c('.Rep1', '.Rep2', '.Rep1', '.Rep2'))
  norm.factors <- calculate.normalisation(
    frequency.fits[replicates], abs(bait.data$fragDist))$factors
  # Extract correlations for raw and normalised data
  counts <- as.matrix(bait.data[,c('cond1.1', 'cond1.2', 'cond2.1', 'cond2.2')])
  raw.cor <- extract.correlation.metrics(counts)
  raw.cor$norm <- 'raw'
  norm.cor <- extract.correlation.metrics(counts / norm.factors)
  norm.cor$norm <- 'normalised'
  # Merge and return data
  comb.cor <- rbind(raw.cor, norm.cor)
  comb.cor$cond1 <- bait.data$cond1[1]
  comb.cor$cond2 <- bait.data$cond2[1]
  comb.cor$bait <- bait.data$baitID[1]
  return(comb.cor)
}

# Calculate correlation metrics
correlation.metrics <- lapply(
  common.counts,
  function(comp) {
    lapply(
      comp,
      calculate.normalisation.correlation)})
correlation.metrics <- rbindlist(do.call(c, correlation.metrics))
correlation.metrics$norm <- factor(
  correlation.metrics$norm,
  levels=c('raw', 'normalised'))
correlation.metrics$class <- factor(
  correlation.metrics$class,
  levels=c('intra.replicate', 'inter.replicate'))
pdf(file.path(outdir, 'count_correlation.pdf'), height=7, width=10)
ggplot(correlation.metrics, aes(x=class, y=correlation, fill=norm)) +
  geom_violin(draw_quantiles=c(0.25,0.5,0.75)) +
  labs(title='Effect Of Normalisation On Count Correlation')
dev.off()
















  


