rm(list=ls())
require('data.table')
require('ggplot2')
require('GenomicRanges')
source('~/github/CaptureC_Analysis/functions/capturec_peaks.R')
params <- list(
  indir='~/differential/results/',
  outdir='~/Desktop/LabMeeting/',
  norm.pattern='.deseq2_norm_results.txt',
  raw.pattern='.deseq2_raw_results.txt',
  dhs='~/differential/DHS-clusters_smooth80_refined50_300bp_min2peaks_annotated.txt',
  sig.padj=0.05,
  max.padj=0.1,
  window=1000)
# Extract input paths
input.paths <- list(
  'raw' = list.files(
    params$indir, pattern=params$raw.pattern, full.names=T),
  'norm' = list.files(
    params$indir, pattern=params$norm.pattern, full.names=T))
names(input.paths[['raw']]) <- gsub(params$raw.pattern, '', basename(input.paths[['raw']]))
names(input.paths[['norm']]) <- gsub(params$norm.pattern, '', basename(input.paths[['norm']]))
input.conditions <- unique(unlist(strsplit(names(input.paths[['raw']]), '\\.'))[c(F,T,T)])
# Read in chicgao data
chicago.dir <- '~/differential/chicago/'
chicago.files <- list.files(chicago.dir, pattern='_Rep1Rep2.Rds', full=T)
chicago.peaks <- lapply(
  chicago.files,
  function(path) {
    data <- readRDS(path)@x
    data[data$score > 5,]
  }
)
chicago.df <- data.frame(
  'baits' = gsub(
    'dm6_TS_Capture_(TSS|CRM)_.*?_Rep1Rep2.Rds', '\\1', basename(chicago.files)),
  'condition' = gsub(
    'dm6_TS_Capture_(TSS|CRM)_(.*?)_Rep1Rep2.Rds', '\\2', basename(chicago.files)),
  'counts' = sapply(chicago.peaks, nrow))
pdf(file.path(params$outdir, 'chicago_peak_no.pdf'), height=4, width=6)
ggplot(chicago.df, aes(x=condition, y=counts, fill=baits)) + 
  geom_col(position='dodge') +
  coord_flip() +
  labs(title='No. Of Peaks Identified By CHiCAGO')
dev.off()

# Extract fragment data
fragData <- read.table(
  '/g/furlong/project/37_Capture-C/analysis/TS_Capture/CRM_all/design/dm6_DpnII.rmap',
  sep='\t', col.names=c('chr', 'start', 'end', 'otherEndID'))
# Create chicago Gr
chicago.grange <- lapply(
  chicago.peaks,
  function(cp) {
    cp.merged <- merge(cp, fragData, by='otherEndID', all.x=T)
    cp.grange <- GRanges(
      seqnames = cp.merged$chr,
      ranges = IRanges(
        start=cp.merged$start,
        end=cp.merged$end),
      strand='*')
    return(cp.grange)})
chicago.condition <- gsub(
  'dm6_TS_Capture_(TSS|CRM)_(.*?)_Rep1Rep2.Rds', '\\2', basename(chicago.files))
# Merge all chicago peaks
chicago.grange.all <- reduce(do.call(c, chicago.grange))
# Create variables to store output
sig.fragments <- data.frame()
# Loop through comparisons and extract data
for (norm in names(input.paths)) {
  for (comp in names(input.paths[[norm]])) {
    # Read in data
    cat(paste0(norm, ' ', comp, '\n'))
    results.data <- fread(input.paths[[norm]][[comp]], showProgress=F)
    # Extract location of dhs for conditions
    cond1 <- strsplit(comp, '\\.')[[1]][2]
    cond2 <- strsplit(comp, '\\.')[[1]][3]
    # extract chicago datasets
    chicago.indices <- grepl(cond1, chicago.condition) | grepl(cond2, chicago.condition)
    chicago.grange.specific <- reduce(do.call(c, chicago.grange[chicago.indices]))
    # Extract significant peaks and convert to grange
    sig.results <- results.data[results.data$padj <= params$sig.padj,]
    sig.results$comp <- comp
    sig.results$norm <- norm
    sig.grange <- makeGRangesFromDataFrame(
      sig.results,
      seqnames.field='fragChr',
      start.field='fragStart',
      end.field='fragEnd')
    # Count overlap with all and specific dhs and store data
    sig.results$chicago <- 'None'
    sig.results$chicago[countOverlaps(sig.grange, chicago.grange.all) > 0] <- 'Any'
    sig.results$chicago[countOverlaps(sig.grange, chicago.grange.specific) > 0] <- 'Specific'
    # Add comparison
    sig.fragments <- rbind(sig.fragments, sig.results)
  }
}
sig.fragments$chicago <- factor(sig.fragments$chicago, levels=c('None', 'Any', 'Specific'))
sig.fragments$time <- as.character(
  grepl('6-8h', sig.fragments$cond1) + grepl('6-8h', sig.fragments$cond2))
sig.fragments$time <- gsub('^0$', '10-12h', sig.fragments$time)
sig.fragments$time <- gsub('^1$', '6-8h/10-12h', sig.fragments$time)
sig.fragments$time <- gsub('^2$', '6-8h', sig.fragments$time)
sig.fragments$time <- factor(sig.fragments$time, levels=c('6-8h', '10-12h', '6-8h/10-12h'))
sig.fragments$norm <- gsub('norm', 'normalised', sig.fragments$norm)
sig.fragments$distance <- cut(
    abs(sig.fragments$fragDist),
    breaks=c(2000, 4000, 8000, 16000, 32000, 64000, 128000, 256000, Inf),
    labels=c('2-4kb', '4-8kb', '8-16kb', '16-32kb', '32-64kb', '64-128kb', '128-256kb', '>256kb'))


pdf(file.path(params$outdir, 'number_significant_peaks_chicago.pdf'), height=4, width=5)
ggplot(sig.fragments, aes(x=norm, fill=chicago)) +
  geom_bar() +
labs(title='Significant Fragments Overlapping CHiCAGO')
dev.off()
pdf(file.path(params$outdir, 'number_significant_peaks_chicago_fill.pdf'), height=4, width=5)
ggplot(sig.fragments, aes(x=norm, fill=chicago)) +
  geom_bar(position='fill') +
  labs(title='Significant Fragments Overlapping CHiCAGO', y='ratio')
dev.off()
pdf(file.path(params$outdir, 'number_significant_peaks_chicago_time.pdf'), height=4, width=8)
ggplot(sig.fragments, aes(x=norm, fill=chicago)) +
  geom_bar() +
  facet_wrap(~time) +
  labs(title='Significant Fragments Overlapping CHiCAGO', y='count')
dev.off()
pdf(file.path(params$outdir, 'number_significant_peaks_chicago_time_fill.pdf'), height=4, width=8)
ggplot(sig.fragments, aes(x=norm, fill=chicago)) +
  geom_bar(position='fill') +
  facet_wrap(~time) +
  labs(title='Significant Fragments Overlapping CHiCAGO', y='ratio')
dev.off()
ggplot(sig.fragments, aes(x=abs(fragDist), col=norm)) +
  geom_density() +
  scale_x_log10(breaks=c(1e4, 1e5, 1e6, 1e7)) +
  























