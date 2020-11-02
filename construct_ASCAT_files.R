packages <- c('data.table', 'tidyr', 'getopt', 'parallel', 'Morpho', 'matrixStats')

invisible(suppressMessages(suppressWarnings(lapply(packages, library, character.only = TRUE))))

testCriteriaFile <- function(opt, file, type) {
  
  if(is.null(opt[[file]])) stop(sprintf('Please supply a %s file', type))
  
  if(!file.exists(opt[[file]])) stop(sprintf('The supplied %s file does not exist', type))
  
  if(file.access(opt[[file]], mode = 4) == -1) stop(sprintf('Missing read permission for the %s file.', type))
  
}

testCriteriaTokenized <- function(opt, variable) {
  
  if(is.null(opt[[variable]])) stop(sprintf('Please supply %s', variable))
  
  if(!grepl(',', opt[[variable]], fixed = TRUE)) stop(sprintf('Please supply a CSV string for %s', variable))
  
  strsplit(gsub(' ', '', opt[[variable]]), split = ',')[[1]]
  
}

prcRegression <- function(covVec, sysVar, chrMeta) {
  
  unlist(Map(function(chrCov, chrSysVar) {mu <- mean(chrCov); mu + lsfit(chrSysVar, chrCov)$residuals}, split(covVec, chrMeta), split(sysVar, chrMeta)))
  
}

concatenateMeta <- function(DT, meta) {
  
  invisible(lapply(colnames(meta), function(x, DT, meta) {set(DT, j = x, value = meta[[x]])}, DT, meta))
  
  setcolorder(DT, colnames(meta))

}

writeResults <- function(idx, cov, baf, controlBaf, indIdx, outdir) {
  
  loc <- c(1:3, idx + 3)
  
  fwrite(cov[, ..loc], file = file.path(outdir, sprintf('%s_tumor_logR.txt', colnames(cov)[idx + 3])), sep = '\t', quote = F)
  
  fwrite(baf[, ..loc], file = file.path(outdir, sprintf('%s_tumor_baf.txt', colnames(baf)[idx + 3])), sep = '\t', quote = F)
  
  fwrite(cov[, ..loc], file = file.path(outdir, sprintf('%s_control_logR.txt', colnames(cov)[idx + 3])), sep = '\t', quote = F)
  
  loc <- c(colnames(controlBaf)[1:3], indIdx[idx])
  
  fwrite(setnames(controlBaf[, ..loc],c(colnames(controlBaf)[1:3], colnames(baf)[idx + 3])), file = file.path(outdir, sprintf('%s_control_baf.txt', colnames(baf)[idx + 3])), sep = '\t', quote = F)
  
}

mcPredefined <- function(threads) {
  function(...) {
    mclapply(..., mc.cores = threads)
  }
}

mapPredefined <- function(threads) {
  function(...) {
    mcmapply(..., SIMPLIFY = F, mc.cores= threads)
  }
}

spec = matrix(c(
  'covFile','c',1,'character',
  'bafFile','b',1,'character',
  'outputDir','o',1,'character',
  'is_control','i',1,'character',
  'is_male','m',1,'character',
  'individual_index','d',1,'character',
  'is_baf_control','l',1,'character',
  'prcomp','p',1,'integer',
  'min_coverage','mc',1,'integer',
  'min_samples','ms',1,'integer',
  'threads','t',1,'integer',
  'help','h',0,'logical'
), byrow=TRUE, ncol=4)

opt <- getopt(spec)

if(!is.null(opt$help) || length(opt) < 7) {
  
  cat(getopt(spec, usage = TRUE))
  
  q(status = 1)
  
}

print('QC...')

testCriteriaFile(opt, 'covFile', 'coverage')

testCriteriaFile(opt, 'bafFile', 'BAF')

if(!file.exists(opt$outputDir)) dir.create(opt$outputDir, recursive = T, mode = '750')

is_control <- ifelse(is.na(as.numeric(is_control <- testCriteriaTokenized(opt, 'is_control'))), as.logical(is_control), as.logical(as.numeric(is_control)))

is_male <- ifelse(is.na(as.numeric(is_male <- testCriteriaTokenized(opt, 'is_male'))), as.logical(is_male), as.logical(as.numeric(is_male)))

individual_index <- testCriteriaTokenized(opt, 'individual_index')

is_baf_control <- ifelse(is.na(as.numeric(is_baf_control <- testCriteriaTokenized(opt, 'is_baf_control'))), as.logical(is_baf_control), as.logical(as.numeric(is_baf_control)))

if(!any(sapply(lapply(list(is_male, individual_index, is_baf_control), length), identical, length(is_control)))) stop(sprintf('Information vectors for is_control, is_male, individual_index, is_baf_control are of unequal length. Same length as is_controls: %s', paste0(outcome, collapse = ',')))

prc <- if(!is.null(opt$prcomp)) max(opt$prcomp, 1)

mc <- ifelse(is.null(opt$min_coverage), 20, max(opt$min_coverage, 1))

ms <- ifelse(is.null(opt$min_samples), 1, max(opt$min_samples, 1))

threads <- ifelse(is.null(opt$threads), 1, max(opt$threads, 1))

print('Reading data...')

cov <- fread(opt$covFile, sep = '\t', header = T)

baf <- fread(opt$bafFile, sep = '\t', header = T)

print('Reviewing data...')

if(!identical(dim(cov), dim(baf))) stop('Dimensions of the Coverage and BAF matrices are different')

if(!identical(colnames(cov), colnames(baf))) if(identical(sort(colnames(cov)), sort(colnames(baf)))) setorderv(baf, colnames(cov)) else stop(sprintf('Difference in the column names: %s', paste0(setdiff(colnames(cov), colnames(baf)), collapse = ',')))

if(!identical((snpIdents <- ifelse(threads > 1, get('mcPredefined')(2), get('lapply'))(meta <- list(baf = baf[, 1:3], cov = cov[, 1:3]), unite, col = 'comb', 1:3))$baf, snpIdents$cov)) stop('Coverage and BAF matrix are not in the same order')

print('Prefiltering data...')

baf[, (1:3) := NULL]

cov[, (1:3) := NULL]

yIdx <- grepl('y|Y|chry|chrY', meta$baf[, chr], perl = T)

if(sum(yIdx) > 0) {
  
  baf <- baf[!yIdx]
  
  cov <- cov[!yIdx]
  
  meta <- lapply(meta, `[`, !yIdx)
  
}

xIdx <- meta$baf[, grep('x|X|chrx|chrX', chr, perl = T)]

if(sum(is_male) > 0 && length(xIdx) > 0) set(cov, i = xIdx, j = colnames(cov)[is_male], value = cov[xIdx, lapply(.SD, `*`, 2), .SDcols = is_male])

covIdx <- cov[, rowSums(.SD >= mc, na.rm = T) >= ms, .SDcols = is_control]

if(sum(!covIdx) > 0) {
  
  baf <- baf[covIdx]
  
  cov <- cov[covIdx]
  
  meta <- lapply(meta, `[`, covIdx)
  
}

print('Normalizing data...')

covNormalized <- ifelse(threads > 1, get('mapPredefined')(threads), get('Map'))(function(bafVec, covVec){ dens <- density(covVec[(bafVec >= 0.4) & (bafVec <= 0.6)], bw = 0.01, from = 0, to = 1000, n = 4096); covVec / dens$x[which.max(dens$y)]}, baf, cov)

invisible(lapply(names(covNormalized), function(x, covNorm, cov) set(cov, j = x, value = covNorm[[x]]), covNormalized, cov))

rm(covNormalized)

if(!is.null(opt$prcomp)) {

  print('PCA analysis...')

  sysVar <- as.data.table((prcompfast(cov[, ..is_control], scale. = T))$x[,1:prc])

  covNormalized <- ifelse(threads > 1, get('mcPredefined')(threads), get('lapply'))(cov, prcRegression, sysVar, meta$baf[['chr']])
  
  invisible(lapply(names(covNormalized), function(x, covNorm, cov) set(cov, j = x, value = covNorm[[x]]), covNormalized, cov))
  
  rm(covNormalized)
  
}

print('Calculating LogR...')

medControl <- cov[, ..is_control][,rowMedians(as.matrix(.SD))]

cov <- cov[, lapply(.SD, function(x) log2(pmax(x / get('medControl', parent.frame()), 1e-300)))]

if(sum(is_male) > 0 && length(xIdx) > 0) {

  print('Correcting X-chromosome...')
	
  xIdx <- meta$baf[, grep('x|X|chrx|chrX', chr, perl = T)]
  
  set(cov, i = xIdx, j = colnames(cov)[is_male], value = cov[xIdx, lapply(.SD, `-`, 1), .SDcols = is_male])

}

print('Calculating control BAF profile...')

controlBaf <- data.table()

invisible(lapply(unique(individual_index), function(x, ii, baf, controlBaf){indIdx <- ii == x; controlBaf[, (x) := baf[, rowMedians(as.matrix(.SD)), .SDcols = indIdx]]}, individual_index, baf, controlBaf))

concatenateMeta(baf, meta$baf)

concatenateMeta(controlBaf, meta$baf)

concatenateMeta(cov, meta$cov)

print('Writing results...')

invisible(ifelse(threads > 1, get('mcPredefined')(threads), get('lapply'))(seq_len(length(individual_index)), writeResults, cov, baf, controlBaf, individual_index, opt$outputDir))
