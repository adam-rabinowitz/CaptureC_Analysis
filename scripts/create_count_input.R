# Import required modules
rm(list=ls())
require(data.table)
require(parallel)

# Set input file paths
rds.files <- list(
  'TSS' = list.files(
    '/g/furlong/project/37_Capture-C/analysis/TS_Capture/TSS_all',
    pattern='Rep\\dRep\\d.Rds$', full.names=T),
  'CRM' = list.files(
    '/g/furlong/project/37_Capture-C/analysis/TS_Capture/CRM_all',
    pattern='Rep\\dRep\\d.Rds$', full.names=T))
bait.files <- list(
  'TSS' = '/g/furlong/project/37_Capture-C/analysis/TS_Capture/TSS_all/design/dm6_DpnII.baitmap',
  'CRM' = '/g/furlong/project/37_Capture-C/analysis/TS_Capture/CRM_all/design/dm6_DpnII.baitmap')
frag.files <- list(
  'TSS' = '/g/furlong/project/37_Capture-C/analysis/TS_Capture/TSS_all/design/dm6_DpnII.rmap',
  'CRM' = '/g/furlong/project/37_Capture-C/analysis/TS_Capture/CRM_all/design/dm6_DpnII.rmap')
# Read in files
outdir <- '/g/furlong/project/37_Capture-C/data/diffinter/input/'

###############################################################################
## Function to create files
###############################################################################
read.rds <- function(rds, bait, frag, outdir) {
  # Extract dataset and create output file
  dataset <- gsub(
    '^.*?((TSS|CRM)_\\d+-\\d+h_[^_]+)_Rep\\dRep\\d.Rds$', '\\1', rds)
  prefix <- gsub(
    '^.*?((TSS|CRM)_\\d+-\\d+h_[^_]+_Rep\\dRep\\d).Rds$', '\\1', rds)
  replicates <- gsub(
    '^.*?(TSS|CRM)_\\d+-\\d+h_[^_]+_(Rep\\dRep\\d).Rds$', '\\2', rds)
  replicates <- c(
    substring(replicates, 1, 4),
    substring(replicates, 5, 8)
  )
  # Read in counts and filter desired columns
  cd <- readRDS(rds)@x
  cols <- c('baitID', 'otherEndID', 'distSign', 'N.1', 'N.2')
  cd <- cd[,cols,with=F]
  names(cd) <- c('baitID', 'fragID', 'fragDist', replicates)
  # Read in bait data and merge
  baitmap <- fread(bait)
  names(baitmap) <- c('baitChr', 'baitStart', 'baitEnd', 'baitID', 'baitName')
  cd.bait <- merge(cd, baitmap, by='baitID', all.x=T)
  # Read in fragment data and merge
  fragmap <- fread(frag)
  names(fragmap) <- c('fragChr', 'fragStart', 'fragEnd', 'fragID')
  cd.frag <- merge(cd.bait, fragmap, by='fragID', all.x=T)
  # Add additional columns
  cd.frag$dataset <- dataset
  # Reorder data and save
  output <- cd.frag[,c(
    'dataset', 'baitName', 'baitID', 'baitChr', 'baitStart', 'baitEnd',
    'fragID', 'fragChr', 'fragStart', 'fragEnd', 'fragDist', replicates),
    with=F]
  output <- output[order(output$baitID, output$fragID),]
  outfile <- file.path(outdir, paste0(prefix, '.counts.txt'))
  write.table(output, outfile, row.names=F, col.names=T, quote=F, sep='\t')
}

###############################################################################
## Calculate background
###############################################################################
# Create files
for (DS in c('CRM', 'TSS')) {
  rds <- rds.files[[DS]]
  bait <- bait.files[[DS]]
  frag <- frag.files[[DS]]
  cores <- length(rds)
  mclapply(
    rds,
    read.rds,
    bait=bait,
    frag=frag,
    outdir=outdir,
    mc.cores=length(rds)
  )
}
