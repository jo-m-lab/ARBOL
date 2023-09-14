#' @rdname data
#' @name data
#' @aliases mockSC
#' @title Synthetic single-cell data
#'
#' @description
#' \code{mockSC()} is designed to generate synthetic single-cell data.
#'   These data are not meant to represent biologically meaningful use-cases,
#'   but are solely intended for use in examples, for unit-testing, and to
#'   demonstrate \code{SPOTlight}'s general functionality.
#'
#' @param ng,nc,nt integer scalar specifying the number
#'   of genes, cells, types (groups) and spots to simulate.
#'
#' @return
#' \itemize{
#' \item{\code{mockSC} returns a \code{Seurat} object
#'   with rows = genes, columns = single cells, and cell metadata
#'   (\code{colData}) column \code{type} containing group identifiers.}
#' }
#'
#' @examples
#' sce <- mockSC()
NULL

#' @rdname data
#' @importFrom SeuratObject CreateSeuratObject
#' @importFrom stats rnbinom runif
#' @export
mockSC <- function(ng = 200, nc = 50, nt = 3) {
    z <- lapply(seq_len(nt), \(t) {
        ms <- 2^runif(ng, 2, 10)
        ds <- 0.5 + 100 / ms
        y <- rnbinom(ng * nc, mu = ms, size = 1 / ds)
        y <- matrix(y, nrow = ng, ncol = nc)
        dimnames(y) <- list(
            paste0("gene", seq_len(ng)),
            paste0("cell", seq_len(nc))
        )
        x <- CreateSeuratObject(counts = y)
        x$type <- factor(
            paste0("type", t),
            paste0("type", seq_len(nt))
        )
        return(x)
    })
    # Return merged seurat object
    suppressWarnings(merge(z[[1]], z[2:length(z)]))
}
