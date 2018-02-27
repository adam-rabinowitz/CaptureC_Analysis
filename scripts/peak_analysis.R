rm(list=ls())
require('data.table')
require('GenomicRanges')
source('~/github/CaptureC_Analysis/functions/capturec_peaks.R')
params <- list(
  indir='~/differential/deseq2Results/',
  outdir='~/differential/deseq2Peaks/',
  pattern='.deseq2_results.txt$',
  dhs='~/differential/DHS-clusters_smooth80_refined50_300bp_min2peaks_annotated.txt',
  sig.padj=0.05,
  max.padj=0.1,
  window=1000)
# Extract input paths
input.paths <- list.files(
  params$indir, pattern=params$pattern, full.names=T)
names(input.paths) <- gsub(params$pattern, '', basename(input.paths))
input.conditions <- unique(unlist(strsplit(names(input.paths), '\\.'))[c(F,T,T)])
# Read dhs data
dhs.data <- read.table(params$dhs, header=T, sep='\t', check.names=F)
dhs.data <- setnames(
  dhs.data,
  c('6-8_mef_peak', '10-12_mef_peak', '6-8_elav_peak', '10-12_elav_peak', '6-8_DN_peak', '10-12_DN_peak'),
  c('6-8h_Mef2', '10-12h_Mef2', '6-8h_elav', '10-12h_elav', '6-8h_DN', '10-12h_DN'))
dhs.data <- dhs.data[,c('chr', 'start', 'end', input.conditions)]
dhs.data <- dhs.data[apply(dhs.data[,4:ncol(dhs.data)], 1, max) > 0,]
dhs.grange <- makeGRangesFromDataFrame(dhs.data, keep.extra.columns=T)
# Create variables to store output
sig.fragments <- data.frame()
sig.peaks <- data.frame()
# Loop through comparisons and extract data
for (comp in names(input.paths)) {
  # Read in data
  cat(paste0(comp, '\n'))
  results.data <- fread(input.paths[[comp]], showProgress=F)
  # Extract location of dhs for conditions
  cond.grange <- dhs.grange[
    apply(
      mcols(dhs.grange)[,strsplit(comp, '\\.')[[1]][2:3]],
      1, min) > 0]
  # Extract significant peaks and convert to grange
  sig.results <- results.data[results.data$padj <= params$sig.padj,]
  sig.grange <- makeGRangesFromDataFrame(
    sig.results,
    seqnames.field='fragChr',
    start.field='fragStart',
    end.field='fragEnd')
  # Count overlap with all and specific dhs and store data
  sig.results$dhs <- 'None'
  sig.results$dhs[countOverlaps(sig.grange, dhs.grange) > 0] <- 'Any'
  sig.results$dhs[countOverlaps(sig.grange, cond.grange) > 0] <- 'Specific'
  sig.fragments <- rbind(sig.fragments, sig.results)
  # Find results peaks
  results.list <- split(results.data, results.data$baitID)
  peaks <- find.results.peaks(
    results=results.list, sig.padj=params$sig.padj,
    max.padj=params$max.padj, window=1000)
  sig.peaks <- rbind(sig.peaks, peaks)
}
ggplot(sig.fragments, aes(x=cond1, fill=dhs)) +
  geom_bar(position='fill') +
  coord_flip()



