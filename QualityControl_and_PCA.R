packages <- c('data.table', 'tidyr', 'getopt', 'parallel', 'rtracklayer', 'Morpho', 'ggplot2', 'naturalsort', 'GenomeInfoDb')

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

plotContribution <- function(outDir, pca) {
  
  pdf(file = file.path(outDir, 'PC_BarPlot.pdf'), width = 2 * length(pca$sdev), height = 10)
  
  print(ggplot(data.frame(id = factor(colnames(pca$rotation), levels = naturalsort(colnames(pca$rotation))), contr = pca$sdev / sum(pca$sdev)), aes(x = id, y = contr)) + geom_bar(stat = 'identity') + labs(x = 'PC number', y = 'Proportion contribution to variance'))
  
  dev.off()
  
}

mcPredefined <- function(threads) {
  function(...) {
    mclapply(..., mc.cores = threads)
  }
}

mapPredefined <- function(threads) {
  function(...) {
    mcmapply(..., mc.cores= threads)
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
  'genome','g',1,'character',
  'threads','t',1,'integer',
  'help','h',0,'logical'
), byrow=TRUE, ncol=4)

opt <- getopt(spec)

if(!is.null(opt$help) || length(opt) < 7) {
  
  cat(getopt(spec, usage = TRUE))
  
  q(status = 1)
  
}

testCriteriaFile(opt, 'covFile', 'coverage')

testCriteriaFile(opt, 'bafFile', 'BAF')

if(!file.exists(opt$outputDir)) dir.create(opt$outputDir, recursive = T, mode = '0750')

is_control <- ifelse(is.na(as.numeric(is_control <- testCriteriaTokenized(opt, 'is_control'))), as.logical(is_control), as.logical(as.numeric(is_control)))

is_male <- ifelse(is.na(as.numeric(is_male <- testCriteriaTokenized(opt, 'is_male'))), as.logical(is_male), as.logical(as.numeric(is_male)))

individual_index <- testCriteriaTokenized(opt, 'individual_index')

is_baf_control <- ifelse(is.na(as.numeric(is_baf_control <- testCriteriaTokenized(opt, 'is_baf_control'))), as.logical(is_baf_control), as.logical(as.numeric(is_baf_control)))

if(!any(sapply(lapply(list(is_male, individual_index, is_baf_control), length), identical, length(is_control)))) stop(sprintf('Information vectors for is_control, is_male, individual_index, is_baf_control are of unequal length. Same length as is_controls: %s', paste0(outcome, collapse = ',')))

prc <- ifelse(is.null(opt$prcomp), 2, max(opt$prcomp, 2))

mc <- ifelse(is.null(opt$min_coverage), 20, opt$min_coverage)

ms <- ifelse(is.null(opt$min_samples), 1, max(opt$min_samples, 1))

threads <- ifelse(is.null(opt$threads), 1, max(opt$threads, 1))

cov <- fread(opt$covFile, sep = '\t', header = T)

baf <- fread(opt$bafFile, sep = '\t', header = T)

if(!identical(dim(cov), dim(baf))) stop('Dimensions of the Coverage and BAF matrices are different')

if(!identical(colnames(cov), colnames(baf))) if(identical(sort(colnames(cov)), sort(colnames(baf)))) setorderv(baf, colnames(cov)) else stop(sprintf('Difference in the column names: %s', paste0(setdiff(colnames(cov), colnames(baf)), collapse = ',')))

if(!identical((snpIdents <- ifelse(threads > 1, get('mcPredefined')(2), get('lapply'))(meta <- list(baf = baf[, 1:3], cov = cov[, 1:3]), unite, col = 'comb', 1:3))$baf, snpIdents$cov)) stop('Coverage and BAF matrix are not in the same order')

baf[, (1:3) := NULL]

cov[, (1:3) := NULL]

yIdx <- grepl('y|Y|chry|chrY', meta$baf[, chr], perl = T)

if(sum(yIdx) > 0) {
  
  baf <- baf[!yIdx]
  
  cov <- cov[!yIdx]
  
  meta <- lapply(meta, `[`, !yIdx)
  
}

xIdx <- grepl('x|X|chrx|chrX', meta$baf[, chr], perl = T)

if(sum(is_male) > 0 && sum(xIdx) > 0) set(cov, j = colnames(cov)[is_male], value = cov[, lapply(.SD, `*`, as.numeric(xIdx) + 1), .SDcols = is_male])

covIdx <- cov[, rowSums(.SD >= mc, na.rm = T) >= ms, .SDcols = is_control]

if(sum(!covIdx) > 0) {
  
  baf <- baf[covIdx]
  
  cov <- cov[covIdx]
  
  meta <- lapply(meta, `[`, covIdx)
  
}

covNormalized <- ifelse(threads > 1, get('mapPredefined')(threads), get('Map'))(function(bafVec, covVec){dens <- density(covVec[(bafVec >= 0.4) & (bafVec <= 0.6)], bw = 0.01, from = 0, to = 1000, n = 4096); covVec / dens$x[which.max(dens$y)]}, baf, cov)

pca <- prcompfast(covNormalized[, is_control], scale. = T)

invisible(plotContribution(opt$outputDir, pca))

invisible(lapply(as.data.frame(combn(1:prc, 2)), function(x, sdev, pca) {pdf(file = file.path(opt$outputDir, sprintf('PC%s_vs_PC%s.pdf', x[1], x[2])), width = 10, height = 10); print(ggplot(as.data.frame(pca$rotation), aes(x = pca$rotation[, sprintf('PC%s', x[1])], y = pca$rotation[, sprintf('PC%s', x[2])])) + geom_point() + labs(x = sprintf('PC%s (%f%%)', x[1], sdev[x[1]] * 100), y = sprintf('PC%s (%f%%)', x[2], sdev[x[2]] * 100))); dev.off()}, pca$sdev / sum(pca$sdev), pca))

gList <- GRangesList(lapply(1:prc, function(idx, pca) {gr <- GRanges(seqnames = meta$baf[, chr], ranges = IRanges(start = meta$baf[, pos], width = 1), score = pca$x[, idx]); seqlevels(gr) <- naturalsort(seqlevels(gr)); return(gr)}, pca))

names(gList) <- as.character(1:prc)

invisible(lapply(names(gList), function(x, gList) export.wig(gList[[x]], sprintf('PC%s.wig', x)), gList))

if(!is.null(opt$genome)) {
  
  library(sprintf('BSgenome.Hsapiens.UCSC.%s', opt$genome), character.only = T)
  
  gList <- GRangesList(lapply(gList, function(x) {seqlevels(x) <- mapSeqlevels(seqlevels(x), 'UCSC'); return(x)}))
  
  gList <- GRangesList(lapply(gList, function(x, Hsapiens) {seqlengths(x) <- seqlengths(Hsapiens)[seqlevels(x)]; genome(x) <- genome(Hsapiens)[seqlevels(x)]; return(x)}, Hsapiens))
  
  invisible(lapply(names(gList), function(x, gList){rtracklayer::export.bw(gList[[x]], sprintf('PC%s.bw', x))}, gList))

}
