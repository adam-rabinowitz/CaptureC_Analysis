# Script to generate distance-decay plots for replicates of all samples
rm(list=ls())
require(ggplot2)
require(reshape2)
source('/g/furlong/adamr/github/CaptureC_Analysis/functions/capturec_normalisation.R')
source('/g/furlong/adamr/github/CaptureC_Analysis/functions/capturec_input.R')
source('~/github/CaptureC_Analysis/functions/capturec_differential.R')
# Find input files
params <- list(
  'workdir' = '/g/furlong/project/37_Capture-C/data/diffinter',
  'countdir' = 'input',
  'chrsizes' = 'chr_sizes.txt',
  'comparisons' = 'comparisons.txt',
  'outdir' = 'correlation'
)
# Extract input file paths
setwd('/g/furlong/project/37_Capture-C/data/diffinter/')
input.paths <- list.files(
  params$countdir,
  pattern='_Rep\\dRep\\d.counts.txt',
  full.names=T
)
# Extract chromosome sizes
chr.sizes <- read.table(
  params$chrsizes, sep='\t', col.names=c('name', 'length')
)
chr.sizes <- split(chr.sizes$length, chr.sizes$name)
# Read comparisons
comparisons <- read.table(
  params$comparisons, header=F, sep='\t', stringsAsFactors=F
)
# Create fits
frequency.fits <- lapply(
  input.paths,
  distance.decay.fit,
  min.dist=2000,
  bin.size=1000,
  k=20,
  chr.sizes=chr.sizes,
  cores=6
)
frequency.fits <- do.call(c, frequency.fits)
# Read input data
input.data <- lapply(
  input.paths,
  function(z) {
    data <- fread(z)
    data.filt <- data[data$baitChr == data$fragChr,]
    data.list <- split(data.filt, data.filt$baitID)
    return(data.list)
  }
)
names(input.data) <- gsub('_Rep1Rep2.counts.txt', '', basename(input.paths))
# Loop through comparisons and merge counts
common.counts <- list()
for (cmp.no in 1:nrow(comparisons)) {
  # Extract data
  cond1 <- as.character(comparisons[cmp.no, 1])
  cond2 <- as.character(comparisons[cmp.no, 2])
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
  # Check names
  if (
    colnames(counts)[1] != colnames(counts)[2] |
    colnames(counts)[3] != colnames(counts)[4]
  ) {
    stop('mismatched replicates')
  }
  # Perform correlation
  log2.counts <- log2(counts)
  log2.cor <- cor(log2.counts)
  if (any(is.na(log2.cor))) {
    print(log2.cor)
    stop('Could not calculate all corrlations')
  }
  # Extract correlation data
  diag(log2.cor) <- NA
  cor.df <- melt(
    log2.cor,
    na.rm=T,
    varnames=c('cond1', 'cond2'),
    value.name='correlation'
  )
  row.names(cor.df) <- 1:nrow(cor.df)
  # Add information on correlation type and return
  cor.df$class <- 'condition'
  cor.df$class[cor.df$cond1 == cor.df$cond2] <- 'replicate'
  cor.df$comp <- paste(colnames(counts)[c(1,3)], collapse='.')
  return(cor.df)
}

###############################################################################
## Extract correlations
###############################################################################
calculate.normalisation.correlation <- function(bait.data) {
  # Extract normalisation factors
  samples <- rep(c(bait.data$cond1[1], bait.data$cond2[1]), each=2)
  #samples <- gsub('elav', 'Elav', samples)
  replicates <- paste0(
    samples,
    c('.Rep1', '.Rep2', '.Rep1', '.Rep2'))
  norm.factors <- calculate.normalisation(
    frequency.fits[replicates], abs(bait.data$fragDist))$factors
  # Extract correlations for raw and normalised data
  counts <- as.matrix(bait.data[,c('cond1.1', 'cond1.2', 'cond2.1', 'cond2.2')])
  colnames(counts) <- samples
  raw.cor <- extract.correlation.metrics(counts)
  raw.cor$norm <- 'raw'
  norm.cor <- extract.correlation.metrics(counts / norm.factors)
  norm.cor$norm <- 'normalised'
  # Merge and return data
  comb.cor <- rbind(raw.cor, norm.cor)
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
# Format and add additional columns to metrics
correlation.metrics$cond1 <- gsub('elav', 'Elav', correlation.metrics$cond1)
correlation.metrics$cond2 <- gsub('elav', 'Elav', correlation.metrics$cond2)
correlation.metrics$time <- sapply(
  strsplit(correlation.metrics$comp, '_'),
  function(z) {
    sum(z=='6-8h')
  }
)
correlation.metrics$time <- factor(
  correlation.metrics$time,
  levels=c(2,0,1),
  labels=c('6-8h', '10-12h', '6-8h/10-12h')
)
correlation.metrics$time <- gsub('^0$', '10-12h',  correlation.metrics$time)
correlation.metrics$time <- gsub('^1$', '6-8h/10-12h',  correlation.metrics$time)
correlation.metrics$time <- gsub('^2$', '6-8h',  correlation.metrics$time)
correlation.metrics$time <- factor(correlation.metrics$time, levels=c('6-8h', '10-12h', '6-8h/10-12h'))
correlation.metrics$class <- factor(correlation.metrics$class, levels=c('replicate', 'condition'))
correlation.metrics$class <- factor(correlation.metrics$class, levels=c('replicate', 'condition'))
correlation.metrics$norm <- factor(
  correlation.metrics$norm,
  levels=c('raw', 'normalised'))
# 
pdf(file.path(outdir, 'raw_count_correlation.pdf'), height=4, width=4)
ggplot(correlation.metrics[correlation.metrics$norm == 'raw',], aes(x=class, y=correlation, fill=norm)) +
  geom_violin(draw_quantiles=c(0.25,0.5,0.75)) +
  scale_y_continuous(limits=c(0.25, 1))
dev.off()
pdf(file.path(outdir, 'raw_count_correlation_time.pdf'), height=4, width=8)
ggplot(correlation.metrics[correlation.metrics$norm == 'raw',], aes(x=class, y=correlation, fill=norm)) +
  geom_violin(draw_quantiles=c(0.25,0.5,0.75)) +
  scale_y_continuous(limits=c(0.25, 1)) +
  facet_wrap(~time)
dev.off()
pdf(file.path(outdir, 'count_correlation_time.pdf'), height=4, width=9)
ggplot(correlation.metrics, aes(x=class, y=correlation, fill=norm)) +
  geom_violin(draw_quantiles=c(0.25,0.5,0.75)) +
  facet_wrap(~time) +
  scale_y_continuous(limits=c(0.25, 1))
dev.off()
pdf(file.path(outdir, 'count_correlation.pdf'), height=4, width=6)
ggplot(correlation.metrics, aes(x=class, y=correlation, fill=norm)) +
  geom_violin(draw_quantiles=c(0.25,0.5,0.75)) +
  scale_y_continuous(limits=c(0.25, 1))
dev.off()

###############################################################################
## Create plot showing raw inter-replicate correlation for all tissues
###############################################################################
rawReplicateCorrelation <- ggplot(
  correlation.metrics[
    (
      correlation.metrics$class == 'replicate' &
      correlation.metrics$norm == 'raw' &
      (
        correlation.metrics$time == '6-8h' |
        correlation.metrics$time == '10-12h'
      )
    )
  ,],
  aes(
    x=gsub('(CRM|TSS)_\\d+-\\d+h_', '', cond1),
    y=correlation,
    fill=substring(cond1, 1, 3)
  )
) +
geom_violin(draw_quantiles=c(0.25,0.5,0.75)) +
facet_wrap(~time) +
scale_y_continuous(limits=c(0.25,1), expand=c(0,0)) +
labs(
  title='raw inter-replicate correlation',
  fill='datastet',
  x='tissue'
)

###############################################################################
## Create plot showing raw inter-sample correlation for all tissues
###############################################################################
rawConditionCorrelationAll <- ggplot(
  correlation.metrics[
    (
      correlation.metrics$class == 'condition' &
      correlation.metrics$norm == 'raw'
    )
  ,],
  aes(
    x=gsub('(CRM|TSS)_\\d+-\\d+h_', '', cond1),
    y=correlation,
    fill=substring(cond1, 1, 3)
  )
) +
geom_violin(draw_quantiles=c(0.25,0.5,0.75)) +
facet_wrap(~time) +
scale_y_continuous(limits=c(0.25,1), expand=c(0,0)) +
labs(
  title='raw inter-sample correlation',
  fill='datastet',
  x='tissue'
)

###############################################################################
## Create plot showing raw inter-sample correlation for Elav and Mef2 tissues
###############################################################################
rawConditionCorrelationElavMef2 <- ggplot(
  correlation.metrics[
    (
      correlation.metrics$class == 'condition' &
      correlation.metrics$norm == 'raw' &
      !grepl('DN$', correlation.metrics$cond1) &
      !grepl('DN$', correlation.metrics$cond2)
    ),
  ],
  aes(
    x=gsub('(CRM|TSS)_\\d+-\\d+h_', '', cond1),
    y=correlation,
    fill=substring(cond1, 1, 3)
  )
) +
  geom_violin(draw_quantiles=c(0.25,0.5,0.75)) +
  facet_wrap(~time) +
  scale_y_continuous(limits=c(0.25,1), expand=c(0,0)) +
  labs(
    title='raw inter-sample correlation (-DN)',
    fill='datastet',
    x='tissue'
)

###############################################################################
## Create plot showing normalised inter-replicate correlation for all tissues
###############################################################################
normReplicateCorrelation <- ggplot(
  correlation.metrics[
    (
      correlation.metrics$class == 'replicate' &
      correlation.metrics$norm == 'normalised' &
      (
        correlation.metrics$time == '6-8h' |
          correlation.metrics$time == '10-12h'
      )
    )
    ,],
  aes(
    x=gsub('(CRM|TSS)_\\d+-\\d+h_', '', cond1),
    y=correlation,
    fill=substring(cond1, 1, 3)
  )
) +
  geom_violin(draw_quantiles=c(0.25,0.5,0.75)) +
  facet_wrap(~time) +
  scale_y_continuous(limits=c(0.25,1), expand=c(0,0)) +
  labs(
    title='normalised inter-replicate correlation',
    fill='datastet',
    x='tissue'
  )

###############################################################################
## Create plot showing normalised inter-sample correlation for all tissues
###############################################################################
normConditionCorrelationAll <- ggplot(
  correlation.metrics[
    (
      correlation.metrics$class == 'condition' &
        correlation.metrics$norm == 'normalised'
    )
    ,],
  aes(
    x=gsub('(CRM|TSS)_\\d+-\\d+h_', '', cond1),
    y=correlation,
    fill=substring(cond1, 1, 3)
  )
) +
  geom_violin(draw_quantiles=c(0.25,0.5,0.75)) +
  facet_wrap(~time) +
  scale_y_continuous(limits=c(0.25,1), expand=c(0,0)) +
  labs(
    title='normalised inter-sample correlation',
    fill='datastet',
    x='tissue'
  )

###############################################################################
## Create plot showing normalised inter-sample correlation for Elav and Mef2 tissues
###############################################################################
normConditionCorrelationElavMef2 <- ggplot(
  correlation.metrics[
    (
      correlation.metrics$class == 'condition' &
        correlation.metrics$norm == 'normalised' &
        !grepl('DN$', correlation.metrics$cond1) &
        !grepl('DN$', correlation.metrics$cond2)
    ),
    ],
  aes(
    x=gsub('(CRM|TSS)_\\d+-\\d+h_', '', cond1),
    y=correlation,
    fill=substring(cond1, 1, 3)
  )
) +
  geom_violin(draw_quantiles=c(0.25,0.5,0.75)) +
  facet_wrap(~time) +
  scale_y_continuous(limits=c(0.25,1), expand=c(0,0)) +
  labs(
    title='normalised inter-sample correlation (-DN)',
    fill='datastet',
    x='tissue'
  )

###############################################################################
## Create plot showing affect of normalisation at 6-8h 
###############################################################################
effectNormReplicateCorrelation <- ggplot(
  correlation.metrics[
    (
      correlation.metrics$class == 'replicate' &
      (
        correlation.metrics$time == '6-8h' |
        correlation.metrics$time == '10-12h'
      )
    ),
    ],
  aes(
    x=gsub('(CRM|TSS)_\\d+-\\d+h_', '', cond1),
    y=correlation,
    fill=norm
  )
) +
  geom_violin(draw_quantiles=c(0.25,0.5,0.75)) +
  facet_wrap(~time) +
  scale_y_continuous(limits=c(0.25,1), expand=c(0,0)) +
  labs(
    title='affect of normalisation on inter-replicate correlation',
    fill='normalisation',
    x='tissue'
  )

###############################################################################
## Create plot showing affect of normalisation at 10-12h
###############################################################################
effectNormConditionCorrelation <- ggplot(
  correlation.metrics[
    (
      correlation.metrics$class == 'condition'
    ),
    ],
  aes(
    x=gsub('(CRM|TSS)_\\d+-\\d+h_', '', cond1),
    y=correlation,
    fill=norm
  )
) +
  geom_violin(draw_quantiles=c(0.25,0.5,0.75)) +
  facet_wrap(~time) +
  scale_y_continuous(limits=c(0.25,1), expand=c(0,0)) +
  labs(
    title='affect of normalisation on inter-condition correlation',
    fill='normalisation',
    x='tissue'
  )

###############################################################################
## Create plot showing affect of normalisation at 6-8h/10-12h
###############################################################################
effectNormConditionCorrelationElavMef2 <- ggplot(
  correlation.metrics[
    (
      correlation.metrics$class == 'condition' &
      !grepl('DN', correlation.metrics$cond1) &
      !grepl('DN', correlation.metrics$cond2)
    ),
    ],
  aes(
    x=gsub('(CRM|TSS)_\\d+-\\d+h_', '', cond1),
    y=correlation,
    fill=norm
  )
) +
  geom_violin(draw_quantiles=c(0.25,0.5,0.75)) +
  facet_wrap(~time) +
  scale_y_continuous(limits=c(0.25,1), expand=c(0,0)) +
  labs(
    title='affect of normalisation on inter-condition correlation (-DN)',
    fill='normalisation',
    x='tissue'
  )


pdf(
  file.path(params$outdir, 'correlation_plots.pdf'),
  height=4, width=8, onefile=T, pointsize=6
)
print(rawReplicateCorrelation)
print(rawConditionCorrelationAll)
print(rawConditionCorrelationElavMef2)
print(normReplicateCorrelation)
print(normConditionCorrelationAll)
print(normConditionCorrelationElavMef2)
print(effectNormReplicateCorrelation)
print(effectNormConditionCorrelation)
print(effectNormConditionCorrelationElavMef2)
dev.off()

cor.df
ggplot(
  cor.df[cor.df$class=='intra',],
  aes(x=cond2, y=correlation)
) +
geom_point()

correlation.metrics
replicate.metrics <- correlation.metrics[correlation.metrics$class=='replicate',]
table(replicate.metrics$cond1)
table(replicate.metrics$cond2)

summary(
  correlation.metrics$correlation[
    !grepl('DN', correlation.metrics$cond1) &
    #!grepl('DN', correlation.metrics$cond2) &
    grepl('Elav', correlation.metrics$cond1) &
    correlation.metrics$norm == 'raw'
  ]
)







  


