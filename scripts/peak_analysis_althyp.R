rm(list=ls())
require(ggplot2)
require('data.table')
require('GenomicRanges')
source('~/github/CaptureC_Analysis/functions/capturec_peaks.R')

###############################################################################
## Function to create dhs grange from input path
###############################################################################
create.dhs.grange <-  function(
  path
) {
  # Read in dhs data and change names
  dhs.data <- read.table(path, header=T, sep='\t', check.names=F)
  dhs.data <- setnames(
    dhs.data,
    c(
      '6-8_mef_peak', '10-12_mef_peak', '6-8_elav_peak', '10-12_elav_peak', 
      '6-8_DN_peak', '10-12_DN_peak'),
    c(
      '6-8h_Mef2', '10-12h_Mef2', '6-8h_elav', '10-12h_elav', '6-8h_DN',
      '10-12h_DN'))
  dhs.data <- dhs.data[,c(
    'chr', 'start', 'end', '6-8h_Mef2', '10-12h_Mef2', '6-8h_elav',
    '10-12h_elav', '6-8h_DN', '10-12h_DN')]
  dhs.data <- dhs.data[apply(dhs.data[,4:ncol(dhs.data)], 1, max) > 0,]
  makeGRangesFromDataFrame(dhs.data, keep.extra.columns=T)
}

###############################################################################
## Function to create chicago grange from input paths
###############################################################################
create.chicago.grange <- function(
  chicago.paths, fragfile, significance
) {
  # Read in frag file
  fragData <- read.table(
    fragfile, sep='\t', col.names=c('chr', 'start', 'end', 'otherEndID')
  )
  # Read in chicago data
  chicago.grange.indv <- lapply(
    chicago.paths,
    function(path) {
      data <- readRDS(path)@x
      sig.data <- data[data$score >= significance,]
      merge.data <- merge(sig.data, fragData, by='otherEndID', all.x=T)
      merge.grange <- makeGRangesFromDataFrame(
        merge.data,
        ignore.strand =T,
        seqnames.field='chr',
        start.field='start',
        end.field='end'
      )
      reduce(merge.grange, min.gapwidth=0)
    }
  )
  # Merge data by tissue and time
  chicago.grange.cond <- list()
  conditions <- gsub('^.*?(\\d+-\\d+h_[^_]+)_.*?$', '\\1', basename(chicago.paths))
  for (cond in unique(conditions)) {
    cond.grange <- chicago.grange.indv[grepl(cond, conditions)]
    cond.grange <- reduce(do.call(c, cond.grange), min.gapwidth=0)
    chicago.grange.cond[[cond]] <- cond.grange
  }
  # Find grange with all peaks
  chicago.grange.all <- chicago.grange.cond
  names(chicago.grange.all) <- NULL
  chicago.grange.all <- do.call(c, chicago.grange.all)
  chicago.grange.all <- reduce(chicago.grange.all, min.gapwidth=0)
  # Add overlaps to all peaks and return data
  for (cond in names(chicago.grange.cond)) {
    mcols(chicago.grange.all)[cond] <- countOverlaps(
      chicago.grange.all,
      chicago.grange.cond[[cond]])
  }
  return(chicago.grange.all)
}

###############################################################################
## Function to count specific and non-specifc overlaps between ranges
###############################################################################
count.overlap.type <- function(
  query, reference, cond1, cond2
) {
  # Check conditions
  cond.intersect <- intersect(
    c(cond1, cond2),
    colnames(mcols(reference))
  )
  if (length(cond.intersect) != 2) {
    stop('conditions not found')
  }
  # Split reference grange
  specific.index <- apply(
    mcols(reference)[,c(cond1, cond2)],
    1,
    sum
  ) > 0
  specific.reference <- reference[specific.index]
  # Extract overlaps
  output <- rep('none', length(query))
  output[countOverlaps(query, reference) > 0] <- 'general'
  output[countOverlaps(query, specific.reference) > 0] <- 'specific'
  return(output)
}

###############################################################################
## Function to read in bait data
###############################################################################
extract.bait.names <- function(path.list) {
  # Read in data
  data.list <- lapply(
    path.list,
    function(z) {
      read.table(
        z, header=F, sep='\t', col.names=c('chr', 'start', 'end', 'id', 'name')
      )
    }
  )
  # Merge data
  bait.data <- unique(do.call(rbind, data.list))
  if (length(bait.data$id) != length(unique(bait.data$id))) {
    stop('non-unique baitID are present')
  }
  # Extract data and return
  bait.names <- as.character(bait.data$name)
  names(bait.names) <- as.character(bait.data$id)
  return(bait.names)
}

###############################################################################
## Perform the analysis
###############################################################################
# Set parametrs
params <- list(
  workdir='/g/furlong/project/37_Capture-C/data/diffinter',
  indir='merged_results/fragments',
  outdir='merged_results/fragments',
  chicagoCRMdir='/g/furlong/project/37_Capture-C/analysis/TS_Capture/CRM_all',
  chicagoTSSdir='/g/furlong/project/37_Capture-C/analysis/TS_Capture/TSS_all',
  chicago.pattern='(Rep\\d+){2,4}\\.Rds',
  fragfile='/g/furlong/project/37_Capture-C/analysis/TS_Capture/CRM_all/design/dm6_DpnII.rmap',
  CRMbaits='/g/furlong/project/37_Capture-C/analysis/TS_Capture/CRM_all/design/dm6_DpnII.baitmap',
  TSSbaits='/g/furlong/project/37_Capture-C/analysis/TS_Capture/TSS_all/design/dm6_DpnII.baitmap',
  dhsfile='/g/furlong/project/53_DHS_variation/data/Charles_IDR/DHS-clusters_smooth80_refined50_300bp_min2peaks_expanded150_annotated.txt',
  chicago.sig=5,
  sig.padj=0.05,
  max.padj=0.1,
  window=1000
)
setwd(params$workdir)
# Set patterns
patterns <- list(
  'norm.above0.0' = 'above0.0.norm_deseq_results.txt$',
  'norm.above1.0' = 'above1.0.norm_deseq_results.txt$',
  'norm.below0.5' = 'below0.5.norm_deseq_results.txt$',
  'raw.above0.0' = 'above0.0.raw_deseq_results.txt$',
  'raw.above1.0' = 'above1.0.raw_deseq_results.txt$',
  'raw.below0.5' = 'below0.5.raw_deseq_results.txt$'
)
# Extract differential result paths
input.paths <- lapply(
  patterns,
  function(pattern) {
    list.files(path=params$indir, pattern=pattern, full=T)
  }
)
# Read in chicago and dhs data
chicago.paths <- c(
  list.files(params$chicagoTSSdir, pattern=params$chicago.pattern, full=T),
  list.files(params$chicagoCRMdir, pattern=params$chicago.pattern, full=T)
)
chicago.grange <- create.chicago.grange(
  chicago.paths, params$fragfile, params$chicago.sig)
dhs.grange <- create.dhs.grange(params$dhsfile)
# Loop through comparisons and extract data
sig.fragments <- data.frame()
sig.peaks <- data.frame()
for (analysis in names(input.paths)) {
  # Extract normalisation type and alt hypothesis
  normalisation <- gsub('(norm|raw).((above|below)\\d\\.\\d)', '\\1', analysis)
  althypothesis <- gsub('(norm|raw).((above|below)\\d\\.\\d)', '\\2', analysis)
  # Loop through all files
  for (path in input.paths[[analysis]]) {
    # Extract comparison
    comp <- gsub(
      '^.*?((CRM|TSS)\\.\\d+-\\d+h_(elav|Mef2|DN)\\.\\d+-\\d+h_(elav|Mef2|DN)).*?$',
      '\\1', basename(path)
    )
    conds <- unlist(strsplit(comp, '\\.'))[2:3]
    cat(paste(normalisation, althypothesis, comp, '\n'))
    # Read in data
    results.data <- fread(path, showProgress=F)
    sig.results <- results.data[results.data$padj <= params$sig.padj,]
    if (nrow(sig.results) == 0) {
      next
    }
    sig.results$comp <- comp
    sig.results$norm <- normalisation
    sig.results$altHyp <- althypothesis 
    sig.grange <- makeGRangesFromDataFrame(
      sig.results,
      seqnames.field='fragChr',
      start.field='fragStart',
      end.field='fragEnd'
    )
    # Find specific and non-specific overlaps
    conds <- unlist(strsplit(comp, '\\.'))[2:3]
    sig.results$dhs <- count.overlap.type(
      sig.grange,
      dhs.grange,
      conds[1],
      conds[2]
    )
    sig.results$chicago <- count.overlap.type(
      sig.grange,
      chicago.grange,
      conds[1],
      conds[2]
    )
    # Add comparison and time
    sig.results$comp <- comp
    time.count <- sum(grepl('6-8h', unlist(strsplit(comp, '\\.'))))
    if (time.count == 0) {
      sig.results$time <- '10-12h'
    } else if (time.count == 1) {
      sig.results$time <- '6-8h/10-12h'
    } else if (time.count == 2) {
      sig.results$time <- '6-8h'
    }
    sig.fragments <- rbind(sig.fragments, as.data.frame(sig.results))
    # # Find results peaks
    # results.list <- split(results.data, results.data$baitID)
    # peaks <- find.results.peaks(
    #   results=results.list, sig.padj=params$sig.padj,
    #   max.padj=params$max.padj, window=params$window)
    # peaks$comp <- comp
    # peaks$norm <- norm
    # sig.peaks <- rbind(sig.peaks, peaks)
  }
}
ggplot(
  sig.fragments[sig.fragments$norm == 'norm',],
  aes(x=abs(lfc), col=altHyp)
) +
  geom_density()

# Save data
write.table(
  sig.fragments, file.path(params$indir, 'significant_fragments.txt'),
  sep='\t', quote=F, row.names=F
)
# Format data for plotting
sig.fragments$conditions <- factor(
  gsub('(CRM|TSS).', '', sig.fragments$comp),
  levels=c(
    '6-8h_DN.6-8h_elav',
    '6-8h_DN.6-8h_Mef2',
    '6-8h_elav.6-8h_Mef2',
    '10-12h_DN.10-12h_elav',
    '10-12h_DN.10-12h_Mef2',
    '10-12h_elav.10-12h_Mef2',
    '6-8h_DN.10-12h_DN',
    '6-8h_elav.10-12h_elav',
    '6-8h_Mef2.10-12h_Mef2'
  )
)

chicago.overlap <- data.frame(
  'lfc' = c(1, 2, 1, 2),
  'padj' = c(0.05, 0.05, 0.01, 0.01),
  'count' = NA,
  'specific' = NA
)
for (i in 1:nrow(chicago.overlap)) {
  lfc <- chicago.overlap$lfc[i]
  padj <- chicago.overlap$padj[i]
  subset.sig.fragments <- sig.fragments[
    sig.fragments$norm == 'norm' &
    abs(sig.fragments$lfc) >= lfc &
    sig.fragments$padj <= padj,]
  chicago.overlap$count[i] <- nrow(subset.sig.fragments)
  chicago.overlap$specific[i] <- sum(subset.sig.fragments$chicago == 'specific')
}
chicago.overlap$ratio <- chicago.overlap$specific / chicago.overlap$count
ggplot(
  chicago.overlap,
  aes(x=as.character(padj), y=as.character(lfc), size=ratio)
) +
  geom_point() +
  scale_size_continuous(
    limits=c(0, max(chicago.overlap$ratio)),
    range=c(1,20)
  )



sig.fragments$dhs <- factor(sig.fragments$dhs, levels=c('none', 'general', 'specific'))
sig.fragments$chicago <- factor(sig.fragments$chicago, levels=c('none', 'general', 'specific'))

# Calculate enrichment of specfic chicago overlaps in specific dhs overlaps
phyper(
  q=sum(sig.fragments$dhs == 'specific' & sig.fragments$chicago == 'specific'),
  m=sum(sig.fragments$chicago == 'specific'),
  n=length(sig.fragments$chicago) - sum(sig.fragments$chicago == 'specific'),
  k=sum(sig.fragments$dhs == 'specific'),
  lower.tail=F
)
print('all')
print(nrow(sig.fragments))
print('dhs.specific')
print(sum(sig.fragments$dhs == 'specific'))
print('chicago.specific')
print(sum(sig.fragments$chicago == 'specific'))
print('both.specific')
print(sum(sig.fragments$chicago == 'specific' & sig.fragments$dhs == 'specific'))


###############################################################################
##
###############################################################################
pdf('dhs_chicago_overlap.pdf', height=4, width=8, onefile=T)
# Plot to display count of DHS overlaps for each comparison
ggplot(
  sig.fragments,
  aes(x=conditions, fill=dhs)
) +
  geom_bar() +
  coord_flip() +
  scale_y_continuous(expand=c(0,0)) +
  facet_wrap(~gsub('(TSS|CRM).*', '\\1', comp) + norm) +
  labs(title='Effect Of Normalisation On DHS Overlap')
# Plot to display ratio of DHS overlaps for each comparison
ggplot(
  sig.fragments,
  aes(x=conditions, fill=dhs)
) +
  geom_bar(position='fill') +
  coord_flip() +
  scale_y_continuous(expand=c(0,0)) +
  facet_wrap(~gsub('(TSS|CRM).*', '\\1', comp) + norm) +
  labs(title='Effect Of Normalisation On DHS Overlap')
# Plot to display count of chicago overlaps for each comparison
ggplot(
  sig.fragments,
  aes(x=conditions, fill=chicago)
) +
  geom_bar() +
  coord_flip() +
  facet_wrap(~gsub('(TSS|CRM).*', '\\1', comp) + norm) +
  labs(title='Effect Of Normalisation On Chicago Overlap')
# Plot to display ratio of chicago overlaps for each comparison
ggplot(
  sig.fragments,
  aes(x=conditions, fill=chicago)
) +
  geom_bar(position='fill') +
  coord_flip() +
  scale_y_continuous(expand=c(0,0)) +
  facet_wrap(~gsub('(TSS|CRM).*', '\\1', comp) + norm) +
  labs(title='Effect Of Normalisation On Chicago Overlap')
# Plot to display count of chicago overlaps for none and specific DHS overlaps
ggplot(
  sig.fragments[sig.fragments$dhs != 'general' & sig.fragments$norm == 'norm',],
  aes(x=conditions, fill=chicago)
) +
  geom_bar() +
  coord_flip() +
  facet_wrap(~dhs) +
  labs(title='Effect Of DHS Overlap On Chicago Overlap')
# Plot to display ratio of chicago overlaps for none and specific DHS overlaps
ggplot(
  sig.fragments[sig.fragments$dhs != 'general' & sig.fragments$norm == 'norm',],
  aes(x=conditions, fill=chicago)
) +
  geom_bar(position='fill') +
  coord_flip() +
  scale_y_continuous(expand=c(0,0)) +
  facet_wrap(~dhs) +
  labs(title='Effect Of DHS Overlap On Chicago Overlap')
dev.off()

sig.peaks <- sig.peaks[sig.peaks$norm == 'norm',]
sig.peaks$peakSigFragmentsBins <- cut(
  sig.peaks$peakSigFragments,
  breaks=c(1:10,Inf),
  right=F,
  labels=c(1:9, '10+'))
pdf(file.path(params$outdir, 'proportion_fragments_peak.pdf'), height=4, width=6)
ggplot(sig.peaks, aes(x=peakSigFragmentsBins)) +
  geom_col(aes(y=peakSigFragments /sum(sig.peaks$peakSigFragments)), fill='darkgreen') +
  labs(
    title='Proportion Of Fragments Within Peak Type',
    y='proportion of fragments',
    x='fragments per peak')
dev.off()
pdf(file.path(params$outdir, 'number_fragments_peak.pdf'), height=4, width=6)
ggplot(sig.peaks, aes(x=peakSigFragmentsBins)) +
  geom_bar(fill='darkgreen') +
  labs(
    title='Signficant Fragments Per Peak',
    y='number of peaks',
    x='fragments per peak')
dev.off()
pdf(file.path(params$outdir, 'size_significant_peaks.pdf'), height=4, width=6)
ggplot(sig.peaks, aes(x=peakEnd - peakStart)) +
  geom_density(fill='darkgreen') +
  scale_x_log10() +
  labs(title='Size Of Peak Regions', x='size of peak')
dev.off()


  scale_x_continuous(limits=c(0.5, 25))


peak.grange <- makeGRangesFromDataFrame(
  sig.peaks,
  seqnames.field='peakChr',
  start.field='peakStart',
  end.field='peakEnd',
  ignore.strand=T,
  keep.extra.columns=T)
peak.overlaps <- data.frame(findOverlaps(peak.grange, peak.grange))
peak.overlaps$queryName <- sig.peaks$peakID[peak.overlaps$queryHits]
peak.overlaps$subjectComp <- sig.peaks$comp[peak.overlaps$subjectHits]
queryOverlapComp <- tapply(peak.overlaps$subjectComp, peak.overlaps$queryName, paste, collapse=';')
queryOverlapCount <- tapply(
  peak.overlaps$subjectComp,
  peak.overlaps$queryName,
  function(z) {
    length(unique(z))})
ggplot(queryOverlapCount, )
sig.peaks$overlapCount <- queryOverlapCount[sig.peaks$peakID] - 1
sig.peaks$overlapCount[is.na(sig.peaks$overlapCount)]
sig.peaks$overlapComp <- queryOverlapComp[sig.peaks$peakID]
sig.peaks$overlapComp[is.na(sig.peaks$overlapComp)]

pdf(file.path(params$outdir, 'peak_overlaps.pdf'), height=4, width=8)
ggplot(sig.peaks, aes(x=overlapCount)) +
  geom_bar(fill='darkgreen') +
  scale_x_continuous(breaks=0:10) +
  labs(title='Overlap Between Peaks', x='overlap count', y='peak count')
dev.off()

min(sig.peaks$)


sig.fragments$dhs <- factor(sig.fragments$dhs, levels=c('None', 'Any', 'Specific'))
sig.fragments$time <- factor(sig.fragments$time, levels=c('6-8h', '10-12h', '6-8h/10-12h'))
sig.fragments$time <- as.character(
  grepl('6-8h', sig.fragments$cond1) + grepl('6-8h', sig.fragments$cond2))
sig.fragments$time <- gsub('^0$', '10-12h', sig.fragments$time)
sig.fragments$time <- gsub('^1$', '6-8h/10-12h', sig.fragments$time)
sig.fragments$time <- gsub('^2$', '6-8h', sig.fragments$time)
sig.fragments$time <- factor(sig.fragments$time, levels=c('6-8h', '10-12h', '6-8h/10-12h'))
sig.fragments$norm <- gsub('norm', 'normalised', sig.fragments$norm)

ggplot(sig.peaks, aes(x=peakSigFragments)) +
  geom_bar() +
  facet_wrap(~norm, scales='free_y')
ggplot(sig.peaks, aes(x=norm, y=peakEnd-peakStart)) +
  geom_violin() +
  scale_y_log10()
summary(sig.peaks$peakEnd - sig.peaks$peakStart)

pdf(file.path(params$outdir, 'number_significant_peaks.pdf'), height=4, width=4)
ggplot(sig.fragments, aes(x=norm, fill=norm)) +
  geom_bar() +
  guides(fill=F) +
  labs(title='Number Of Significant Fragments')
dev.off()
pdf(file.path(params$outdir, 'number_significant_peaks_time.pdf'), height=4, width=8)
ggplot(sig.fragments, aes(x=norm, fill=norm)) +
  geom_bar() +
  guides(fill=F) +
  facet_wrap(~time) +
  labs(title='Number Of Significant Fragments')
dev.off()

pdf(file.path(params$outdir, 'number_significant_peaks_dhs_fill.pdf'), height=4, width=5)
ggplot(sig.fragments, aes(x=norm, fill=dhs)) +
  geom_bar(position='fill') +
  labs(title='Significant Fragments Overlapping DHS', y='ratio')
dev.off()
pdf(file.path(params$outdir, 'number_significant_peaks_dhs.pdf'), height=4, width=5)
ggplot(sig.fragments, aes(x=norm, fill=dhs)) +
  geom_bar() +
  labs(title='Significant Fragments Overlapping DHS', y='count')
dev.off()
ggplot(sig.fragments, aes(x=comp, fill=dhs)) +
  geom_bar() +
  scale_y_continuous() +
  coord_flip() +
  facet_wrap(~norm) +
ggplot(sig.fragments, aes(x=norm, fill=dhs)) +
  geom_bar(position='fill') +
  scale_y_continuous(expand=c(0,0)) +
  facet_wrap(~time)
pdf(file.path(params$outdir, 'number_significant_peaks_dhs_time.pdf'), height=4, width=8)
ggplot(sig.fragments, aes(x=norm, fill=dhs)) +
  geom_bar() +
  scale_y_continuous() +
  facet_wrap(~time) +
  labs(title='Significant Fragments Overlapping DHS', y='count')
dev.off()
pdf(file.path(params$outdir, 'number_significant_peaks_dhs_time_fill.pdf'), height=4, width=8)
ggplot(sig.fragments, aes(x=norm, fill=dhs)) +
  geom_bar(position='fill') +
  scale_y_continuous() +
  facet_wrap(~time) +
  labs(title='Significant Fragments Overlapping DHS', y='ratio')
dev.off()
pdf(file.path(params$outdir, 'distance_significant_peaks_time.pdf'), height=4, width=8)
ggplot(sig.fragments, aes(x=abs(fragDist), col=norm)) +
  geom_density() +
  scale_x_log10(breaks=c(1e4, 1e5, 1e6, 1e7)) +
  facet_wrap(~time) +
  labs(title='Distance Of Significant Fragments', x='distance from bait')
dev.off()











comp.lost <- (tapply(sig.fragments$norm, sig.fragments$comp, function(z) {sum(z=='norm')}) - 
tapply(sig.fragments$norm, sig.fragments$comp, function(z) {sum(z=='raw')}))
comp.lost[order(comp.lost)]









