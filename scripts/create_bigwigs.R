rm(list=ls())
require(data.table)
source('~/github/CaptureC_Analysis/functions/capturec_input.R')
source('~/github/CaptureC_Analysis/functions/capturec_normalisation.R')
source('~/github/CaptureC_Analysis/functions/capturec_differential.R')
# Set parameters
params <- list(
  min.dist=2000,
  bin.size=1000,
  cores=8,
  k=20,
  indir='/g/furlong/project/37_Capture-C/data/diffinter/input',
  outdir='/g/furlong/project/37_Capture-C/data/diffinter/tracks',
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
    params$indir, pattern='TSS_.*?_Rep1Rep2.counts.txt$', full.names=T
  ),
  'CRM' = list.files(
    params$indir, pattern='CRM_.*?_Rep1Rep2.counts.txt$', full.names=T
  )
)
# Create fits
dataset.fits <- lapply(
  input.paths,
  function(dataset) {
    fits <- lapply(
      dataset,
      distance.decay.fit,
      min.dist=params$min.dist,
      bin.size=params$bin.size,
      k=params$k,
      chr.sizes=chr.sizes,
      cores=params$cores
    )
    do.call(c, fits)
  }
)
# Extract counts
dataset.counts <- lapply(
  input.paths,
  function(dataset) {
    lapply(
      dataset,
      read.split.replicates,
      min.dist=params$min.dist,
      intra.only=T
    )
  }
)
# Process baits sequentially for each dataset
for (dataset in names(dataset.counts)) {
  # Process each bait sequentially
  baits <- unique(c(sapply(dataset.counts[[dataset]], names)))
  for (bait in baits) {
    bait.list <- lapply(dataset.counts[[dataset]], function(z) z[[bait]])
    # Find common fragments for each bait
    bait.frags <- lapply(
      bait.list,
      function(z) {
        z$fragID[pmin(z$N1, z$N2) > 0]
      }
    )
    common.frags <- Reduce(intersect, bait.frags)
    # Extract matrix of common counts for each fragment
    common.counts <- lapply(
      bait.list,
      function(z) {
        counts <- data.frame(z[,c('N1', 'N2')])
        row.names(counts) <- as.character(z$fragID)
        colnames(counts) <- paste0(
          rep(z$dataset[1], 2),
          c('.Rep1', '.Rep2')
        )
        counts <- counts[as.character(common.frags),]
      }
    )
    common.counts <- do.call(cbind, common.counts)
    # Calculate norm matrix
    common.distances <- bait.list[[1]]$fragDist[
      match(common.frags, bait.list[[1]]$fragID)]
    norm.factors <- calculate.normalisation(
      dataset.fits[[dataset]][colnames(common.counts)],
      abs(common.distances)
    )
    # Calculate size factors
    size.factors <- DESeq2::estimateSizeFactorsForMatrix(
      common.counts / norm.factors$factors
    )
    # Process each bait sequentially
    for (bait.data in bait.list) {
      # Extract replicate names
      replicates <- paste0(
        rep(bait.data$dataset[1], 2),
        c('.Rep1', '.Rep2')
      )
      # Extract replicate counts and normalise to distance
      replicate.counts <- bait.data[,c('N1', 'N2')]
      # Calculate and apply distance normalisation 
      norm.factors <- calculate.normalisation(
        dataset.fits[[dataset]],
        abs(bait.data$fragDist)
      )$factors
      norm.counts <- replicate.counts / norm.factors[,replicates]
      # Apply size factors
      size.counts <- t(t(norm.counts) / size.factors[replicates])
      colnames(size.counts) <- replicates
      # Extract interval data
      interval.data <- data.frame(
        'chr' = bait.data$fragChr,
        'start' = bait.data$fragStart - 1,
        'end' = bait.data$fragEnd
      )
      # Create output file prefix and header
      prefix <- file.path(
        params$outdir,
        paste(
          gsub('_', '.', bait.data$dataset[1]),
          bait.data$baitID[1],
          sep='.'
        )
      )
      # Save replicate 1 bedgraph
      rep1.bg <- interval.data
      rep1.bg$score <- size.counts[,1]
      rep1.bg <- rep1.bg[rep1.bg$score > 0,]
      write.table(
        rep1.bg, file=paste0(prefix, '.Rep1.bedGraph'), row.names=F,
        col.names=F, sep='\t', quote=F)
      system(
        paste(
          params$wig.script,
          paste0(prefix, '.Rep1.bedGraph'),
          params$chr.file,
          paste0(prefix, '.Rep1.bw')
        )
      )
      # Save replicate 2 bedgraph
      rep2.bg <- interval.data
      rep2.bg$score <- size.counts[,2]
      rep2.bg <- rep2.bg[rep2.bg$score > 0,]
      write.table(
        rep2.bg, file=paste0(prefix, '.Rep2.bedGraph'), row.names=F,
        col.names=F, sep='\t', quote=F)
      system(
        paste(
          params$wig.script,
          paste0(prefix, '.Rep2.bedGraph'),
          params$chr.file,
          paste0(prefix, '.Rep2.bw')
        )
      )
      # Save combined bedgraph
      comb.bg <- interval.data
      comb.bg$score <- apply(size.counts, 1, mean)
      write.table(
        comb.bg, file=paste0(prefix, '.Comb.bedGraph'), row.names=F,
        col.names=F, sep='\t', quote=F)
      system(
        paste(
          params$wig.script,
          paste0(prefix, '.Comb.bedGraph'),
          params$chr.file,
          paste0(prefix, '.Comb.bw')
        )
      )
    }
  }
}





