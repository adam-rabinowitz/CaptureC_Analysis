# Script to generate distance-decay plots for replicates of all samples
rm(list=ls())
require(ggplot2)
require(reshape2)
source('~/github/CaptureC_Analysis/functions/capturec_normalisation.R')
# Find input files
indir <- '~/differential/newCounts/'
outdir <- '~/Desktop/LabMeeting/'
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
# Extract data
fitted.data <- lapply(
  frequency.fits,
  function(fit) {
  data.frame(
    'distance' = fit$knots,
    'frequency' = fit$penalties)})
fitted.mean <- Reduce('+', fitted.data) / length(fitted.data)
# Create plot data
fitted.plot <- data.frame()
for (fit.name in names(fitted.data)) {
  sample <- unlist(strsplit(fit.name, '\\.'))[1]
  replicate <- unlist(strsplit(fit.name, '\\.'))[2]
  if (replicate == 'Rep1') {
    fitted.plot <- rbind(
      fitted.plot,
      data.frame(
        'sample'=sample,
        'replicate'='Mean',
        'distance'=fitted.mean$distance,
        'frequency'=fitted.mean$frequency))}
  fitted.plot <- rbind(
    fitted.plot,
    data.frame(
      'sample'=sample,
      'replicate'=replicate,
      'distance'=fitted.data[[fit.name]]$distance,
      'frequency'=fitted.data[[fit.name]]$frequency))
}
plot.list <- list()
y.lim <- c(min(fitted.plot$frequency), max(fitted.plot$frequency))
pdf(file.path(outdir, 'combined_fit.pdf'), height=6, width=12)
ggplot(fitted.plot, aes(x=distance, y=frequency, col=replicate)) +
  geom_line() +
  facet_wrap(~sample) +
  scale_x_log10() +
  scale_y_log10(limits=y.lim) +
  scale_color_manual(values=c('black', 'red', 'green')) +
  labs(title='All Samples')
dev.off()
for (sample in unique(fitted.plot$sample)) {
  subset.plot <- fitted.plot[fitted.plot$sample == sample,]
  plot <- ggplot(subset.plot, aes(x=distance, y=frequency, colour=replicate)) +
    geom_line() +
    scale_x_log10() +
    scale_y_log10(limits=y.lim) +
    labs(title=sample) +
    scale_color_manual(values=c('black', 'red', 'green'))
  pdf(file.path(outdir, paste0(sample, '_distance_decay.pdf')), height=4, width=6)
  print(plot)
  dev.off()
}
# Save plots to file
pdf <- pdf(
  file.path(outdir, 'individual_fit_plot.pdf'), height=7, width=10, onefile=T, )
x <- lapply(plot.list, print)
dev.off()










