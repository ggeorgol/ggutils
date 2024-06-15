#' Imports `computeMatrix` output into R
#'
#' @description
#' This script imports a txt.gz matrix generated
#' from deepTools's `computeMatrix` function
#' and parses the infomation and density matrix
#' into a list object
#'
#' @details
#' depends on the data.table::fread function and GenomicRanges
#' It converts the regions into a GRanges object. For conversion
#' to sparse it requires Matrix
#'
#' @param file Path to `computeMatrix` output txt.gz matrix
#' @param sparse Logical. Convert matrix to sparse matrix.
#'
#' @return A list object
#'

import_computeMatrix <- function(file, sparse = FALSE) {
    suppressMessages(require(data.table))
    suppressMessages(require(GenomicRanges))
    suppressMessages(require(Matrix))
    dens.matrix <- fread(file, data.table = FALSE)
    parsed.list <- colnames(dens.matrix)[1]
    parsed.list <- gsub("@|\\{|\\}","",parsed.list)
    parsed.list <- gsub('\\:','\\=',parsed.list)
    parsed.list <- gsub('\\[','c(',parsed.list)
    parsed.list <- gsub('\\]','\\)', parsed.list)
    parsed.list <- gsub('false','FALSE',parsed.list)
    parsed.list <- gsub('true','TRUE',parsed.list)
    parsed.list <- gsub('null','NULL',parsed.list)
    parsed.list <- eval(parse(text = paste0('list(',parsed.list,')')))
    granges <- dens.matrix[,1:4]
    colnames(granges) <- c('chr','start','end','ID')
    parsed.list$granges <- makeGRangesFromDataFrame(df = granges,
            seqnames.field = 'chr',
            start.field = 'start',
            end.field = 'end',
            keep.extra.columns = TRUE)
    data.matrix <- as.matrix(dens.matrix[,7:ncol(dens.matrix)])
    rownames(data.matrix) <- parsed.list$granges$ID
    if (sparse) {
        data.matrix <- Matrix::Matrix(data.matrix, sparse = TRUE)
    }

    # Define sample boundaries
    sample_boundaries <- parsed.list$sample_boundaries
    sample_boundaries <- as.data.frame(cbind((sample_boundaries + 1)[1:(length(sample_boundaries)-1)], sample_boundaries[-1]))
    colnames(sample_boundaries) <- c('start','end')
    sample_boundaries <- apply(sample_boundaries, 1, function(x){as.list(x)})
    names(sample_boundaries) <- parsed.list$sample_labels

    # Define group boundaries
    group_boundaries <- parsed.list$group_boundaries
    group_boundaries <- as.data.frame(cbind((group_boundaries + 1)[1:(length(group_boundaries)-1)], group_boundaries[-1]))
    colnames(group_boundaries) <- c('start','end')
    group_boundaries <- apply(group_boundaries, 1, function(x){as.list(x)})
    names(group_boundaries) <- parsed.list$group_labels

    # Split data.matrix by sample and group

    data.matrix.list <- list()

    for (i in seq_along(sample_boundaries)) {
        data.matrix.list[[i]] <- list()
        names(data.matrix.list)[i] <- names(sample_boundaries)[i]
        for (j in seq_along(group_boundaries)) {
            data.matrix.list[[i]][[j]] <- data.matrix[group_boundaries[[j]]$start:group_boundaries[[j]]$end, sample_boundaries[[i]]$start:sample_boundaries[[i]]$end]
            names(data.matrix.list[[i]])[j] <- names(group_boundaries)[j]
        }
    }


    parsed.list$data <- data.matrix.list
    return(parsed.list)
}
