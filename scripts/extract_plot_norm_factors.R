# Script to generate distance-decay plots for replicates of all samples
rm(list=ls())
require(ggplot2)
require(reshape2)
source('~/github/CaptureC_Analysis/functions/capturec_normalisation.R')
source('~/github/CaptureC_Analysis/functions/capturec_differential.R')
# Find input files
indir <- '~/differential/newCounts/'
outdir <- '~/Desktop/LabMeeting/'
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
# Create breaks for plotting
breaks <- 10^(seq(
  log10(min(frequency.fits[[1]]$knots)),
  log10(max(frequency.fits[[1]]$knots)),
  length.out=100
))
# Calculate factors for each comparison
comparisons <- read.table(comp.file, sep='\t', stringsAsFactors=F)
plot.data <- list()
for (cmp.no in 1:nrow(comparisons)) {
  # Extract comparison name
  cond1 <- comparisons[cmp.no, 1]
  cond2 <- comparisons[cmp.no, 2]
  if (substring(cond1, 1, 3) != substring(cond2, 1, 3)) {
    stop('unexpected comparison')}
  cmp.name <- paste(
    substring(cond1, 1, 3),
    substring(cond1, 5),
    substring(cond2, 5),
    sep='.')
  # Extract replicates
  replicates <- paste0(
    rep(comparisons[cmp.no,], each=2),
    c('.Rep1', '.Rep2', '.Rep1', '.Rep2'))
  # Extract fitted values and store
  factors <- calculate.normalisation(frequency.fits[replicates], x=breaks)
  # Generate data for plots
  output <- lapply(
    factors,
    function(d) {
      d <- as.data.frame(d)
      d$distance <- breaks
      d <- reshape2::melt(d, id.vars='distance', variable.name ='sample')
      d$comparison <- cmp.name
      condition <- as.character(d$sample)
      condition <- replace(condition, condition == replicates[1], 'Cond1.1')
      condition <- replace(condition, condition == replicates[2], 'Cond1.2')
      condition <- replace(condition, condition == replicates[3], 'Cond2.1')
      condition <- replace(condition, condition == replicates[4], 'Cond2.2')
      d$condition <- condition
      return(d)})
  plot.data[[cmp.name]] <- output
}
# Seperate out fit and pplot data
plot.data <- do.call(c, plot.data)
fitted.data <- rbindlist(plot.data[c(T,F)])
fitted.data$time <- NA
fitted.data$time[grepl('6-8h', fitted.data$comparison)] <- '6-8h'
fitted.data$time[grepl('10-12h', fitted.data$comparison)] <- '10-12h'
fitted.data$time[
  grepl('6-8h', fitted.data$comparison) & 
  grepl('10-12h', fitted.data$comparison)] <- '6-8h/10-12h'
factor.data <- rbindlist(plot.data[c(F,T)])
factor.data$time <- NA
factor.data$time[grepl('6-8h', factor.data$comparison)] <- '6-8h'
factor.data$time[grepl('10-12h', factor.data$comparison)] <- '10-12h'
factor.data$time[
  grepl('6-8h', factor.data$comparison) & 
  grepl('10-12h', factor.data$comparison)] <- '6-8h/10-12h'
# Create plots for comparison fits
fit.y.lim <- c(min(fitted.data$value), max(fitted.data$value))
fit.plots <- list()
fit.plots[['combined']] <- ggplot(fitted.data, aes(x=distance, y=value, col=condition)) +
  geom_line() +
  facet_wrap(~comparison) +
  scale_x_log10() +
  scale_y_log10(limits=fit.y.lim) +
  scale_color_manual(values=c('brown1', 'brown3', 'dodgerblue1', 'dodgerblue3')) +
  labs(title='All Comparisons', y='probability')
for (comp in unique(fitted.data$comparison)) {
  comp.data <- fitted.data[fitted.data$comparison == comp,]
  fit.plots[[comp]] <- ggplot(comp.data, aes(x=distance, y=value, col=sample)) +
    geom_line() +
    scale_x_log10() +
    scale_y_log10(limits=fit.y.lim) +
    scale_color_manual(values=c('brown1', 'brown3', 'dodgerblue1', 'dodgerblue3')) +
    labs(title=comp, y='probability')
}
# Create plots for comparison fits
factor.y.lim <- c(
  1/(2^max(abs(log2(factor.data$value)))),
  2^max(abs(log2(factor.data$value))))
factor.y.breaks <- 2^seq(log2(factor.y.lim[1]), log2(factor.y.lim[2]), length.out=7)
factor.y.breaks <- round(factor.y.breaks, 3)
factor.plots <- list()
factor.plots[['combined']] <- ggplot(factor.data, aes(x=distance, y=value, col=condition)) +
  geom_line() +
  facet_wrap(~comparison) +
  scale_x_log10() +
  scale_y_log10(limits=factor.y.lim, breaks=factor.y.breaks) +
  scale_color_manual(values=c('brown1', 'brown3', 'dodgerblue1', 'dodgerblue3')) +
  labs(title='All Comparisons', y='normalisation factor')
for (comp in unique(factor.data$comparison)) {
  comp.data <- factor.data[factor.data$comparison == comp,]
  factor.plots[[comp]] <- ggplot(comp.data, aes(x=distance, y=value, col=sample)) +
    geom_line() +
    scale_x_log10() +
    scale_y_log10(limits=factor.y.lim, breaks=factor.y.breaks) +
    scale_color_manual(values=c('brown1', 'brown3', 'dodgerblue1', 'dodgerblue3')) +
    labs(title=comp, y='normalisation factor')
}
pdf(file.path(outdir, 'example_distance_decay_comparison.pdf'), height=4, width=8)
print(fit.plots[[2]])
dev.off()
pdf(file.path(outdir, 'all_distance_decay_comparison.pdf'), height=6, width=10)
print(fit.plots[[1]])
dev.off()
pdf(file.path(outdir, 'example_distance_decay_factors.pdf'), height=4, width=8)
print(factor.plots[[2]])
dev.off()
pdf(file.path(outdir, 'all_distance_decay_factors.pdf'), height=6, width=10)
print(factor.plots[[1]])
dev.off()
pdf <- pdf(
  file.path(outdir, 'comparison_fit_plot.pdf'), height=7, width=10, onefile=T, )
  x <- lapply(fit.plots, print)
dev.off()
pdf <- pdf(
  file.path(outdir, 'comparison_factor_plot.pdf'), height=7, width=10, onefile=T, )
  x <- lapply(factor.plots, print)
dev.off()
f

head(factor.data)


ggplot(factor.data[factor.data$time == '6-8h',], aes(x=distance, y=value, col=condition)) +
  geom_line() +
  facet_wrap(~comparison) +
  scale_x_log10() +
  scale_y_log10(limits=factor.y.lim, breaks=factor.y.breaks) +
  scale_color_manual(values=c('brown1', 'brown3', 'dodgerblue1', 'dodgerblue3')) +
  labs(title='6-8h', y='normalisation factor')
ggplot(factor.data[factor.data$time == '10-12h',], aes(x=distance, y=value, col=condition)) +
  geom_line() +
  facet_wrap(~comparison) +
  scale_x_log10() +
  scale_y_log10(limits=factor.y.lim, breaks=factor.y.breaks) +
  scale_color_manual(values=c('brown1', 'brown3', 'dodgerblue1', 'dodgerblue3')) +
  labs(title='10-12h', y='normalisation factor')
pdf('')
ggplot(factor.data[factor.data$time == '6-8h/10-12h',], aes(x=distance, y=value, col=condition)) +
  geom_line() +
  facet_wrap(~comparison) +
  scale_x_log10() +
  scale_y_log10(limits=factor.y.lim, breaks=factor.y.breaks) +
  scale_color_manual(values=c('brown1', 'brown3', 'dodgerblue1', 'dodgerblue3')) +
  labs(title='6-8h/10-12h', y='normalisation factor')
pdf(file.path(outdir, '6-8h.distance_decay_fits.pdf'), height=4, width=8)
ggplot(fitted.data[fitted.data$time == '6-8h',], aes(x=distance, y=value, col=condition)) +
  geom_line() +
  facet_wrap(~comparison) +
  scale_x_log10() +
  scale_y_log10(limits=fit.y.lim) +
  scale_color_manual(values=c('brown1', 'brown3', 'dodgerblue1', 'dodgerblue3')) +
  labs(title='6-8h', y='probability')
dev.off()
pdf(file.path(outdir, '10-12h.distance_decay_fits.pdf'), height=4, width=8)
ggplot(fitted.data[fitted.data$time == '10-12h',], aes(x=distance, y=value, col=condition)) +
  geom_line() +
  facet_wrap(~comparison) +
  scale_x_log10() +
  scale_y_log10(limits=fit.y.lim) +
  scale_color_manual(values=c('brown1', 'brown3', 'dodgerblue1', 'dodgerblue3')) +
  labs(title='10-12h', y='probability')
dev.off()
pdf(file.path(outdir, '6-8h.10-12h.distance_decay_fits.pdf'), height=4, width=8)
ggplot(fitted.data[fitted.data$time == '6-8h/10-12h',], aes(x=distance, y=value, col=condition)) +
  geom_line() +
  facet_wrap(~comparison) +
  scale_x_log10() +
  scale_y_log10(limits=fit.y.lim) +
  scale_color_manual(values=c('brown1', 'brown3', 'dodgerblue1', 'dodgerblue3')) +
  labs(title='6-8h/10-12h', y='probability')
dev.off()











