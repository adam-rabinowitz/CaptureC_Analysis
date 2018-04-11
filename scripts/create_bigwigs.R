rm(list=ls())
require(data.table)
source('/g/furlong/adamr/github/CaptureC_Analysis/functions/capturec_input.R')
source('/g/furlong/adamr/github/CaptureC_Analysis/functions/capturec_normalisation.R')
source('/g/furlong/adamr/github/CaptureC_Analysis/functions/capturec_differential.R')

###############################################################################
## Function to create bigwig from bedgraph data
###############################################################################
createBigWig <- function(
  bg, bg.path, chr.path, bw.script
) {
  # Check bedgraph path
  if (!endsWith(bg.path, '.bedGraph')) {
    stop('Unrecognised suffix for bedgraph file path')
  }
  # Escape brackets
  bg.path <- gsub('(\\(|\\))', '\\\\\\1', bg.path)
  # Save bedgraph to file
  write.table(
    bg, file=bg.path, row.names=F, col.names=F, sep='\t', quote=F
  )
  # Create path to bigwig file and create command to convert
  bw.path <- gsub('.bedGraph$', '.bw', bg.path)
  command <- paste(bw.script, bg.path, chr.path, bw.path, '&&', 'rm', bg.path)
  system(command, wait=T)
}

###############################################################################
## Function to calculate size factors
###############################################################################
calculate.size.factors <- function(
  bait.data, distance.fits
) {
  # Extract count data for bait and find common fragments
  counts <- gen.frag.counts(bait.data)
  common.indices <- which(apply(counts, 1, min) > 0)
  common.counts <- counts[common.indices,]
  common.distances <- abs(bait.data$fragDist)[common.indices]
  # Calculate size factors and return
  common.norm.factors <- calculate.normalisation(
    distance.fits[colnames(common.counts)],
    common.distances
  )$factors
  size.factors <- DESeq2::estimateSizeFactorsForMatrix(
    common.counts / common.norm.factors
  )
  return(size.factors)
}

###############################################################################
## Main script to generate bigwigs
###############################################################################
# Set parameters
params <- list(
  mindist=2000,
  binsize=1000,
  cores=4,
  k=20,
  indir='/g/furlong/project/37_Capture-C/data/diffinter/input',
  outdir='/g/furlong/project/37_Capture-C/data/diffinter/tracks',
  pattern='_(Rep\\d+){2,}.counts.txt',
  chr.file='/g/furlong/project/37_Capture-C/data/diffinter/chr_sizes.txt',
  wig.script='/g/furlong/adamr/bedGraphToBigWig'
)
# Read in chromosome files
chr.sizes <- read.csv(
params$chr.file, sep='\t', header=F, col.names=c('chr', 'size')
)
chr.sizes <- split(chr.sizes$size, chr.sizes$chr)
# Extract inut paths
input.paths <- list(
  'TSS' = list.files(
    params$indir, pattern='TSS_.*?_(Rep\\d+){2,}\\.counts\\.txt$', full.names=T
  ),
  'CRM' = list.files(
    params$indir, pattern='CRM_.*?_(Rep\\d+){2,}\\.counts\\.txt$', full.names=T
  )
)
# Process baits sequentially for each dataset
for (dataset in names(input.paths)) {
  # Create fits
  distance.fits <- distance.decay.fit.all(
    input.paths[[dataset]],
    min.dist=params$mindist,
    bin.size=params$binsize,
    k=params$k,
    chr.sizes=chr.sizes,
    cores=params$cores
  )
  # Read in counts for each condition
  merged.counts <- extract.merged.counts(
    input.paths[[dataset]],
    suffix=params$pattern,
    min.dist=params$mindist,
    intra.only=T,
    cores=params$cores
  )
  # Process each bait sequentially
  for (bait in names(merged.counts)) {
    # Extract count data for bait
    bait.data <- merged.counts[[bait]]
    bait.name <- unique(bait.data$baitName)
    if (length(bait.name) != 1) {
      stop('non-unique bait name')
    }
    reps <- colnames(bait.data)[grepl('\\.Rep\\d+$', colnames(bait.data))]
    counts <- gen.frag.counts(bait.data)
    # Calculate size factors
    size.factors <- calculate.size.factors(
      bait.data=bait.data,
      distance.fits=distance.fits
    )
    # Calculate distance normalisation and normalise counts
    norm.factors <- calculate.normalisation(
      distance.fits[colnames(counts)],
      abs(bait.data$fragDist)
    )$factors
    norm.counts <- counts / norm.factors
    size.counts <- t(t(norm.counts) / size.factors[colnames(norm.counts)])
    # Extract interval data
    interval.data <- data.frame(
      'chr' = bait.data$fragChr,
      'start' = bait.data$fragStart - 1,
      'end' = bait.data$fragEnd
    )
    # Extract all conditions and spli replicate by conditions
    condition.list <- split(
      colnames(counts),
      gsub('\\.Rep\\d+$', '', colnames(counts))
    )
    # Loop through conditions and save data to bedgraph
    for (condition in names(condition.list)) {
      # Create output file prefix
      prefix <- file.path(
        params$outdir,
        paste(
          bait.name, condition, sep='_'
        )
      )
      # Extract counts for condition
      condition.reps <- condition.list[[condition]]
      condition.counts <- size.counts[,condition.reps]
      # Create mean bedgraph
      mean.bg <- interval.data
      mean.bg$score <- apply(condition.counts, 1, mean)
      mean.bg <- mean.bg[mean.bg$score > 0,]
      mean.path <- paste0(prefix, '_Mean.bedGraph')
      createBigWig(
        bg=mean.bg, bg.path=mean.path, chr.path=params$chr.file,
        bw.script=params$wig.script
      )
      # Loop through replicates and create individual bedgraphs
      for (cr in condition.reps) {
        rep.no <- gsub('^.*?\\.(Rep\\d+$)', '\\1', cr)
        rep.bg <- interval.data
        rep.bg$score <- condition.counts[,cr]
        rep.bg <- rep.bg[rep.bg$score > 0,]
        rep.path <- paste0(prefix, '_', rep.no, '.bedGraph')
        createBigWig(
          bg=rep.bg, bg.path=rep.path, chr.path=params$chr.file,
          bw.script=params$wig.script
        )
      }
    }
  }
}



