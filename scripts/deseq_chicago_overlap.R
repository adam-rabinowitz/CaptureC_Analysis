require(ggplot2)
require(data.table)

###############################################################################
## Function to read in chicago data
###############################################################################
read.chicago.paths <- function(
  chicago.paths
) {
  # Loop through chicago paths
  chicago.data <- lapply(
    chicago.paths,
    function(paths) {
      output <- lapply(
        paths,
        function(path) {
          # Read in chicago data and select and rename columns
          cdata <- as.data.frame(readRDS(path)@x)
          cdata <- cdata[,c('baitID', 'otherEndID', 'log.p', 'log.q', 'score')]
          colnames(cdata) <- c('baitID', 'fragID', 'log.p', 'log.q', 'score')
          # Convert bait and frag ID to characters
          cdata$baitID <- as.character(cdata$baitID)
          cdata$fragID <- as.character(cdata$fragID)
          # Split data by bait and return
          cdata.list <- split(cdata, cdata$baitID)
          return(cdata.list)
        }
      )
      # Rename data with tissue and time and return
      names(output) <- gsub('^.*?(\\d+-\\d+h_[^_]+).*?$', '\\1', basename(paths))
      return(output)
    }
  )
  return(chicago.data)
}

###############################################################################
## Function to extract data from chicago datasets
###############################################################################
extract.chicago.metrics <- function(
  chicago.data, dataset, cond1, cond2, baitID, fragID, metric, calculation
) {
  # Extract bait data
  bait.data <- list(
    'cond1' = chicago.data[[dataset]][[cond1]][[baitID]],
    'cond2' = chicago.data[[dataset]][[cond2]][[baitID]]
  )
  if (!all(sapply(bait.data, is.data.frame))) {
    stop('bait data has not been succesfully extracted')
  }
  # Extract metric data
  if (!metric %in% colnames(bait.data[[1]])) {
    stop('metric not found in data')
  }
  metric.data <- lapply(
    bait.data,
    function(z) {
      z[base::match(fragID, z$fragID), metric]
    }
  )
  metric.data <- unlist(metric.data)
  if (all(is.na(metric.data))) {
    print(c(dataset, cond1, cond2, baitID, fragID))
    stop('all metrics are NA')
  }
  # Extract fragment data
  result <- calculation(metric.data, na.rm=T)
  return(result)
}

###############################################################################
## Extract chicago metrics for results file
###############################################################################
extract.chicago.metrics.results <- function(
  results, chicago.data, metric, calculation
) {
  # Create output variable
  output <- numeric(nrow(results))
  # Loop through results, populate output and return
  for (i in 1:nrow(results)) {
    output[i] <- extract.chicago.metrics(
      chicago.data=chicago.data,
      dataset=results$dataset[i],
      cond1=results$cond1[i],
      cond2=results$cond2[i],
      baitID=as.character(results$baitID[i]),
      fragID=as.character(results$fragID[i]),
      metric=metric,
      calculation=calculation
    )
  }
  return(output)
}

###############################################################################
## Functions to extract distance between bait and fragment by ID
###############################################################################
find.distance <- function(
  chr1, start1, end1, chr2, start2, end2
) {
  # Set distance as infinite for different chromosomes
  if (chr1 != chr2) {
    distance = Inf
  # Calculate distance for identical chromsome
  } else {
    # Calculate distance
    distance = max(
      start1 - end2,
      start2 - end1,
      1
    ) - 1
    # Adjust for sign
    if (start1 < start2) {
      distance <- distance * -1
    }
  }
  # Return data
  return(distance)
}

###############################################################################
## Function to extract distances between multiple baits and fragments
###############################################################################
find.distances <- function(
  id1, id2, frag.data
) {
  # Check id1 and id2 are the same length
  if (length(id1) != length(id2)) {
    stop('id1 and id2 must be the same length')
  }
  # Extract data
  id1.indices <- match(id1, frag.data$id)
  id2.indices <- match(id2, frag.data$id)
  id1.chr <- frag.data$chr[id1.indices]
  id1.start <- frag.data$start[id1.indices]
  id1.end <- frag.data$end[id1.indices]
  id2.chr <- frag.data$chr[id2.indices]
  id2.start <- frag.data$start[id2.indices]
  id2.end <- frag.data$end[id2.indices]
  # Create output variable, populate and return
  distances <- numeric(length(id1))
  for (i in seq_along(distances)) {
    distances[i] <- find.distance(
      id1.chr[i], id1.start[i], id1.end[i],
      id2.chr[i], id2.start[i], id2.end[i]
    )
  }
  return(distances)
}

###############################################################################
## Extract chicago distance metrics
###############################################################################
extract.chicago.distance.metrics <- function(
  chicago.data, min.score, frag.data
) {
  # Order scores
  scores <- scores[order(scores)]
  # Loop thrpugh data
  lapply(
    chicago.data,
    function(dataset) {
      lapply(
        dataset,
        function(tissue) {
          lapply(
            tissue,
            function(bait) {
              bait <- bait[bait$score >= min.score,]
              bait$dist <- find.distances(
                id1=bait$baitID,
                id2=bait$fragID,
                frag.data=frag.data
              )
              return(bait)
            }
          )
        }
      )
    }
  )
}

###############################################################################
## Perform analysis
###############################################################################
# Extract chicago file paths and read in data
chicago.paths <- list(
  'CRM' = list.files(
    '/g/furlong/project/37_Capture-C/analysis/TS_Capture/CRM_all/',
    pattern='_(Rep\\d+){2,}.Rds$', full.names=T
  ),
  'TSS' = list.files(
    '/g/furlong/project/37_Capture-C/analysis/TS_Capture/TSS_all/',
    pattern='_(Rep\\d+){2,}.Rds$', full.names=T
  )
)
chicago.data <- read.chicago.paths(chicago.paths)
# Import frgament data
frag.data <- read.table(
  '/g/furlong/project/37_Capture-C/analysis/TS_Capture/CRM_all/design/dm6_DpnII.rmap',
  header=F, sep='\t', col.names=c('chr', 'start', 'end', 'id')
)
# Extract distance data for chicago
chicago.dist.data <- extract.chicago.distance.metrics(
  chicago.data=chicago.data, min.score=3, frag.data=frag.data
)
chicag.dist.df <- data.table::rbindlist(
  do.call(c, do.call(c, chicago.dist.data))
)
chicag.dist.df <- chicag.dist.df[is.finite(chicag.dist.df$dist),]
# Read in significant fragments and format
sig.results <- read.table(
  '/g/furlong/project/37_Capture-C/data/diffinter/five_hypotheses_merged/fragments_old/significant_fragments_fdr_0.1.txt',
  header=T, sep='\t', stringsAsFactors=F
)
sig.results <- sig.results[sig.results$padj < 0.05,]
sig.results$altHyp = paste(sig.results$altHyp, sig.results$lfcThr, sep='_')
sig.results$lfcThr <- NULL
# Add chicago score data
sig.results$chicagoScore <- extract.chicago.metrics.results(
  sig.results, chicago.data, metric='score', calculation=max
)
sig.results$chicagoLogp <- extract.chicago.metrics.results(
  sig.results, chicago.data, metric='log.p', calculation=min
)
sig.results$chicagoLogq <- extract.chicago.metrics.results(
  sig.results, chicago.data, metric='log.q', calculation=min
)

###############################################################################
## Create plot comparing distance profiles for chicago and deseq
###############################################################################
deseq.chicago.dist.df <- data.frame(
  rbind(
    data.frame(
      'dataset' = 'chicago_>3',
      'dist' = chicag.dist.df$dist[chicag.dist.df$score > 3]
    ),
    data.frame(
      'dataset' = 'chicago_>5',
      'dist' = chicag.dist.df$dist[chicag.dist.df$score > 5]
    ),
    data.frame(
      'dataset' = 'chicago_>10',
      'dist' = chicag.dist.df$dist[chicag.dist.df$score > 10]
    ),
    data.frame(
      'dataset' = 'deseq_>0',
      'dist' = sig.results$fragDist[
        sig.results$altHyp == 'greaterAbs_0' & sig.results$normalisation == 'distance_size']
    ),
    data.frame(
      'dataset' = 'deseq_<1',
      'dist' = sig.results$fragDist[
        sig.results$altHyp == 'lessAbs_1' & sig.results$normalisation == 'distance_size']
    )
  )
)
pdf(
  '/g/furlong/project/37_Capture-C/data/diffinter/five_hypotheses_merged/plots/distribution_of_significant_fragments.pdf',
  'height'=7, width=10, onefile=T
)
ggplot(
  deseq.chicago.dist.df,
  aes(x=dataset)
) +
  geom_bar(fill='darkgreen', stat='count') +
  geom_text(stat='count', aes(label=..count..), vjust=-0.2) +
  labs(
    'title'='Number Of Significant Fragments: DESeq2 and Chicago'
  )
ggplot(
  deseq.chicago.dist.df,
  aes(x=abs(dist), col=dataset)
) +
  geom_density() +
  scale_x_log10(limits=c(1000, 4e7)) +
  labs(
    'x' = 'absolute distance from bait',
    'title' = 'Distribution Of Significant Fragments: DESeq2 and Chicago'
  )
ggplot(
  sig.results[
    sig.results$altHyp %in% c('greaterAbs_0', 'lessAbs_1') &
      sig.results$normalisation == 'distance_size'
    ,],
  aes(x=abs(fragDist), y=chicagoScore)
) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~altHyp) +
  geom_smooth() +
  labs(
    x='absolute distance',
    y='chicago score',
    title='Effect Of Distance On Chicago Score'
  )
ggplot(
  sig.results[
    sig.results$altHyp == 'greaterAbs_0' &
      sig.results$normalisation == 'distance_size' &
      abs(sig.results$fragDist) > 1e6
    ,],
  aes(x=baitName)
) +
  geom_bar(fill='darkgreen') +
  labs(
    x='bait name',
    title='Number Of Distal (>1e6) greaterAbs_0 Significant Fragments Per Bait'
  )
ggplot(
  sig.results[
    sig.results$altHyp %in% c('greaterAbs_0', 'lessAbs_1') &
    sig.results$normalisation == 'distance_size'
  ,],
  aes(x=altHyp, fill=chicagoScore >= 3)
) +
  geom_bar(position='fill') +
  labs(
    x = 'alternative hypothesis',
    y = 'ratio',
    fill='overlap',
    title='Chicago Overlap: Score >= 3'
  )
ggplot(
  sig.results[
    sig.results$altHyp %in% c('greaterAbs_0', 'lessAbs_1') &
      sig.results$normalisation == 'distance_size'
    ,],
  aes(x=altHyp, fill=chicagoScore >= 5)
) +
  geom_bar(position='fill') +
  labs(
    x = 'alternative hypothesis',
    y = 'ratio',
    fill='overlap',
    title='Chicago Overlap: Score >= 5'
  )

###############################################################################
## Show how chicago score varies for significant deseq interactions
###############################################################################
pdf(
  '/g/furlong/project/37_Capture-C/data/diffinter/five_hypotheses_merged/plots/deseq_chicago_overlap_score.pdf',
  height=7, width=10, onefile=T
)
ggplot(
  sig.results[
    sig.results$altHyp %in% c('greaterAbs_0', 'lessAbs_1') &
      sig.results$normalisation == 'distance_size'
    ,],
  aes(x=altHyp, fill=cut(chicagoScore, breaks=c(0,3,5,10,Inf), include.lowest=T))
) +
  geom_bar(stat='count') +
  scale_y_continuous(breaks=c(seq(0, 10000, 1000), 10000, 20000, 30000)) +
  theme(panel.grid.minor = element_blank()) +
  labs(
    x = 'alternative hypothesis',
    y = 'count',
    fill='chicago score',
    title='Distribution Of Chicago Score For Significant DESeq Fragments'
  )
ggplot(
  sig.results[
    sig.results$altHyp %in% c('greaterAbs_0', 'lessAbs_1') &
      sig.results$normalisation == 'distance_size'
    ,],
  aes(x=altHyp, fill=cut(chicagoScore, breaks=c(0,3,5,10,Inf), include.lowest=T))
) +
  geom_bar(position='fill') +
  labs(
    x = 'alternative hypothesis',
    y = 'ratio',
    fill='chicago score',
    title='Distribution Of Chicago Score For Significant DESeq Fragments'
  )
dev.off()



###############################################################################
## Plot distances of chicago data
###############################################################################
chicago.dist.count <- function(chicago.data, scores) {
  
}

################################
ggplot(
  sig.results,
  aes(x=altHyp),
) +
  geom_bar() +
  facet_wrap(~normalisation) +
  labs('')

###############################################################################
## Plot Effect of normalisation on signifcant interactions
###############################################################################
pdf(
  '/g/furlong/project/37_Capture-C/data/diffinter/five_hypotheses_merged/plots/alternate_hypothesis_plots.pdf',
  height=7, width=10, onefile=T
)
ggplot(
  sig.results,
  aes(x=altHyp, fill=normalisation)
) +
  geom_bar(position='dodge') +
  labs(
    x='alternative hypothesis',
    title='Significant Fragments For Normalisation Procedures'
  )
ggplot(
  sig.results[sig.results$normalisation == 'distance_size',],
  aes(x=abs(fragDist)),
) +
  geom_density() +
  facet_wrap(~altHyp, scales='free_y') +
  scale_x_log10(breaks=c(1e4, 1e5, 1e6, 1e7)) +
  labs(
    x='absolute distance',
    title='Distance Of Significant Fragments For Alternative Hypotheses'
  )
ggplot(
  sig.results[
    sig.results$altHyp %in% c('greaterAbs_0', 'lessAbs_1') &
    sig.results$normalisation == 'distance_size'
  ,],
  aes(x=repMean, col=altHyp),
) +
  geom_density() +
  scale_x_log10() +
  labs(
    x='mean replicate counts',
    title='Mean Counts For Significant Fragments Of Alternate Hypotheses'
  )
ggplot(
  sig.results[
    sig.results$altHyp %in% c('greaterAbs_0', 'lessAbs_1') &
    sig.results$normalisation == 'distance_size'
  ,],
  aes(x=altHyp, y=chicagoScore)
) +
  geom_violin() +
  scale_y_log10(breaks=c(1,3,10,30,100)) +
  labs(
    x='alternative hypothesis',
    y='chicago score',
    title='Distibution Of Chicago Scores For Alternate Hypotheses'
  )
ggplot(
  sig.results[
    sig.results$altHyp %in% c('greaterAbs_0', 'lessAbs_1') &
      sig.results$normalisation == 'distance_size'
    ,],
  aes(x=altHyp, y=abs(chicagoLogp))
) +
  geom_violin() +
  scale_y_log10(breaks=c(1,3,10,30,100)) +
  labs(
    x='alternative hypothesis',
    y='chicago -log10(pvalue)',
    title='Distibution Of Chicago P-Values For Alternate Hypotheses'
  )
dev.off()

distal.greaterAbs_0 <- with(
  sig.results[
    sig.results$altHyp == 'greaterAbs_0' &
    sig.results$normalisation == 'distance_size' &
    abs(sig.results$fragDist) > 1e6
    ,],
  table(baitName)
)


baitNames.greaterAbs.1 <- with(
  sig.results[
    sig.results$altHyp == 'greaterAbs.1' &
    sig.results$normalisation == 'distance_size'
  ,],
  table(baitName)
)
baitNames.greaterAbs.1
baitNames.greaterAbs.1['GluRIA_155003'] / sum(baitNames.greaterAbs.1)

sig.results[
  sig.results$altHyp == 'greaterAbs.1' &
  sig.results$baitName == 'VT42009.1_257616' 
,]






