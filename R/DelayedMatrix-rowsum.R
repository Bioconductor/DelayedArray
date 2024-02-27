### =========================================================================
### rowsum() and colsum() methods for DelayedMatrix objects
### -------------------------------------------------------------------------
###
### These methods are block processed.
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Simple helpers to support verbose operations
###

.summarize_matrix <- function(x)
{
    ans <- paste0("<", paste0(dim(x), collapse=" x "), "> ", class(x)[[1L]])
    if (is.object(x))
        ans <- paste0(ans, " object")
    ans
}

final_hstrip_noop <- function(init, i, grid) {
    what <- .summarize_matrix(init)
    #message("  ", wmsg("| result for horizontal strip ",
    #                   i, "/", nrow(grid), ": ", what))
    message("  ", wmsg("| result is a ", what, " --> returning as-is"))
    init
}

final_vstrip_noop <- function(init, j, grid) {
    what <- .summarize_matrix(init)
    #message("  ", wmsg("| result for vertical strip ",
    #                   j, "/", ncol(grid), ": ", what))
    message("  ", wmsg("| result is a ", what, " --> returning as-is"))
    init
}

realize_matrix <- function(x, BACKEND, verbose)
{
    if (verbose) {
        what <- .summarize_matrix(x)
        message("  ", wmsg("| realizing ", what, " as ",
                           BACKEND, " object ..."),
                " ", appendLF=FALSE)
    }
    ans <- realize(x, BACKEND=BACKEND)
    if (verbose)
        message("ok")
    ans
}

write_full_sink_rows <- function(sink, sink_grid, i, block, verbose)
{
    if (verbose) {
        what <- .summarize_matrix(block)
        message("  ", wmsg("| writing ", what, " to ",
                           class(sink), " object ..."),
                " ", appendLF=FALSE)
    }
    sink <- write_block(sink, sink_grid[[i, 1L]], block)
    if (verbose)
        message("ok")
    sink
}

write_full_sink_cols <- function(sink, sink_grid, j, block, verbose)
{
    if (verbose) {
        what <- .summarize_matrix(block)
        message("  ", wmsg("| writing ", what, " to ",
                           class(sink), " object ..."),
                " ", appendLF=FALSE)
    }
    sink <- write_block(sink, sink_grid[[1L, j]], block)
    if (verbose)
        message("ok")
    sink
}

### 'strip_results' is guaranteed to be a list of length >= 1.
combine_strip_results <- function(fname, strip_results, verbose)
{
    res1 <- strip_results[[1L]]
    if (length(strip_results) == 1L)
        return(res1)
    if (verbose) {
        message(wmsg("=== FINAL STEP ==="))
        if (is.matrix(res1)) {
            what <- "matrices"
        } else {
            what <- paste0(class(res1)[[1L]], " objects")
        }
        message("  ", wmsg("| ", fname, "()'ing strip results ",
                           "(", length(strip_results), " ", what, ") ",
                           "together ..."),
                " ", appendLF=FALSE)
    }
    FUN <- match.fun(fname)
    ans <- do.call(FUN, strip_results)
    if (verbose) {
        message("ok")
        message("=== DONE ===\n")
    }
    ans
}

shared_sink_as_DelayedArray <- function(sink, verbose)
{
    if (verbose) {
        message(wmsg("=== FINAL STEP ==="))
        message("  ", wmsg("| turning ", class(sink), " object ",
                           "into DelayedArray object ..."),
                " ", appendLF=FALSE)
    }
    close(sink)
    ans <- as(sink, "DelayedArray")
    if (verbose) {
        message("ok")
        message("=== DONE ===\n")
    }
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Helpers for BLOCK_rowsum() and BLOCK_colsum()
###

### Whether 'BACKEND' is compatible with the "shared sink" route (see below
### in this file for what the "shared sink" route is).
.compatible_BACKEND <- function(BACKEND)
{
    if (is.null(BACKEND))
        return(FALSE)
    ## Same check as in load_BACKEND_package().
    if (!isSingleString(BACKEND))
        stop(wmsg("'BACKEND' must be a single string or NULL"))
    ## write_block() method for RleRealizationSink objects is broken (it
    ## ignores the 'viewport' argument!) so, until this is fixed, the
    ## RleArray realization backend is not compatible.
    BACKEND != "RleArray"
}

### 'input_grid' must be a 2D grid.
### Returns a 2D grid 'sink_grid' that verifies:
###   (a) refdim(sink_grid)[[1]] == refdim(input_grid)[[1]];
###   (b) refdim(sink_grid)[[2]] == sink_ncol;
###   (c) the blocks on 'sink_grid' are made of full rows and they align
###       with the horizontal strips on 'input_grid'.
### The consequences of (c) are that ncol(sink_grid) == 1 and
### nrow(sink_grid) == length(sink_grid) == nrow(input_grid).
.make_sink_grid_of_hstrips <- function(input_grid, sink_ncol)
{
    stopifnot(is(input_grid, "ArrayGrid"),
              length(dim(input_grid)) == 2L,
              isSingleInteger(sink_ncol))
    if (is(input_grid, "ArbitraryArrayGrid")) {
        tickmarks <- list(input_grid@tickmarks[[1L]], sink_ncol)
        ArbitraryArrayGrid(tickmarks)
    } else {
        refdim <- c(refdim(input_grid)[[1L]], sink_ncol)
        spacings <- c(nrow(input_grid[[1L]]), sink_ncol)
        RegularArrayGrid(refdim, spacings=spacings)
    }
}

### 'input_grid' must be a 2D grid.
### Returns a 2D grid 'sink_grid' that verifies:
###   (a) refdim(sink_grid)[[1]] == sink_nrow;
###   (b) refdim(sink_grid)[[2]] == refdim(input_grid)[[2]];
###   (c) the blocks on 'sink_grid' are made of full columns and they align
###       with the vertical strips on 'input_grid'.
### The consequences of (c) are that nrow(sink_grid) == 1 and
### ncol(sink_grid) == length(sink_grid) == ncol(input_grid).
.make_sink_grid_of_vstrips <- function(input_grid, sink_nrow)
{
    stopifnot(is(input_grid, "ArrayGrid"),
              length(dim(input_grid)) == 2L,
              isSingleInteger(sink_nrow))
    if (is(input_grid, "ArbitraryArrayGrid")) {
        tickmarks <- list(sink_nrow, input_grid@tickmarks[[2L]])
        ArbitraryArrayGrid(tickmarks)
    } else {
        refdim <- c(sink_nrow, refdim(input_grid)[[2L]])
        spacings <- c(sink_nrow, ncol(input_grid[[1L]]))
        RegularArrayGrid(refdim, spacings=spacings)
    }
}

### Note that each block on 'sink_grid' is a horizontal strip made of one or
### more rows.
.sink_chunking_is_compatible_with_hstrips <- function(sink_chunkdim, sink_grid)
{
    stopifnot(is(sink_grid, "ArrayGrid"),
              length(dim(sink_grid)) == 2L,
              ncol(sink_grid) == 1L)
    if (is.null(sink_chunkdim))  # no-chunking
        return(TRUE)
    stopifnot(is.integer(sink_chunkdim), length(sink_chunkdim) == 2L)
    if (all(sink_chunkdim == refdim(sink_grid)))
        return(TRUE)
    ## Dumb heuristic: We consider incompatible chunks that are taller than
    ## the first block in 'sink_grid'.
    ## FIXME: This could certainly be improved/refined.
    ## Anyway, the most important thing for now is that it covers the
    ## worst-case scenario, which is when the sink uses a storage layout
    ## that is column-oriented (e.g. TENxRealizationSink object),
    ## and 'sink_grid' has more than one horizontal strip.
    ## So whatever heuristic we use, we want to make sure that it returns
    ## FALSE in this case.
    if (sink_chunkdim[[1L]] <= nrow(sink_grid[[1L]]))
        return(TRUE)
    FALSE
}

### Note that each block on 'sink_grid' is a vertical strip made of one or
### more columns.
.sink_chunking_is_compatible_with_vstrips <- function(sink_chunkdim, sink_grid)
{
    stopifnot(is(sink_grid, "ArrayGrid"),
              length(dim(sink_grid)) == 2L,
              nrow(sink_grid) == 1L)
    if (is.null(sink_chunkdim))  # no-chunking
        return(TRUE)
    stopifnot(is.integer(sink_chunkdim), length(sink_chunkdim) == 2L)
    ## We treat the "single big chunk" case like the no-chunking case.
    ## Note that the "single big chunk" situation only happens for very
    ## small sinks in which case the chunking does not significantly impact
    ## the writing performance. However, treating this situation as compatible
    ## with the sink vertical strips is convenient when testing things like
    ## BLOCK_rowsum(..., BACKEND="HDF5Array") on a small toy dataset.
    if (all(sink_chunkdim == refdim(sink_grid)))
        return(TRUE)
    ## Dumb heuristic: We consider incompatible chunks that are wider than
    ## the first block in 'sink_grid'.
    ## FIXME: This could certainly be improved/refined.
    ## Anyway, the most important thing for now is that it covers the
    ## worst-case scenario, which is when the sink uses a storage layout
    ## that is row-oriented (i.e. is the transposed of what is used by a
    ## TENxRealizationSink object), and 'sink_grid' has more than one
    ## vertical strip. Whatever heuristic we use, we want to make sure that
    ## it returns FALSE in this case.
    if (sink_chunkdim[[2L]] <= ncol(sink_grid[[1L]]))
        return(TRUE)
    FALSE
}

### Returns a "shared sink" if we can take the "shared sink" route, or NULL
### if we can't.
.make_shared_sink_along_hstrips <- function(BACKEND, sink_grid,
                                            sink_rownames, sink_colnames,
                                            input_grid, BPPARAM)
{
    stopifnot(nrow(sink_grid) == nrow(input_grid), ncol(sink_grid) == 1L)

    ## Note that, at the moment, we don't try the "shared sink" route if
    ## parallel processing is enabled because there's no guarantee that the
    ## realization sink will support concurrent writes (e.g. HDF5 does not).
    ## TODO (maybe):
    ## - For registered realization backends, we could register
    ##   their ability to do concurrent writes, and decide based on that.
    ## - Alternatively, we could introduce a new generic (e.g.
    ##   supports_concurrent_writing() or concurrent_writes(), to define
    ##   in RealizationSink-class.R) with a method defined for RealizationSink
    ##   objects that returns FALSE. Then concrete subclasses that support
    ##   concurrent writes (e.g. TileDBRealizationSink?) would overwrite it
    ##   with a method that returns TRUE.
    if (!is.null(BPPARAM))
        return(NULL)
    if (!.compatible_BACKEND(BACKEND) || nrow(input_grid) == 1L)
        return(NULL)
    sink_dimnames <- list(sink_rownames, sink_colnames)
    sink <- RealizationSink(BACKEND, refdim(sink_grid), dimnames=sink_dimnames,
                                     type="double")
    ## We take the "shared sink" route only if the chunks are "compatible"
    ## with the writing of full sink rows by callback function FINAL()
    ## below (this callback function will get called at the end of processing
    ## each horizontal strip).
    ok <- .sink_chunking_is_compatible_with_hstrips(chunkdim(sink), sink_grid)
    if (ok) sink else NULL
}

### Returns a "shared sink" if we can take the "shared sink" route, or NULL
### if we can't.
.make_shared_sink_along_vstrips <- function(BACKEND, sink_grid,
                                            sink_rownames, sink_colnames,
                                            input_grid, BPPARAM)
{
    stopifnot(nrow(sink_grid) == 1L, ncol(sink_grid) == ncol(input_grid))

    if (!is.null(BPPARAM))
        return(NULL)
    if (!.compatible_BACKEND(BACKEND) || ncol(input_grid) == 1L)
        return(NULL)
    sink_dimnames <- list(sink_rownames, sink_colnames)
    sink <- RealizationSink(BACKEND, refdim(sink_grid), dimnames=sink_dimnames,
                                     type="double")
    ## We take the "shared sink" route only if the chunks are "compatible"
    ## with the writing of full sink columns by callback function FINAL()
    ## below (this callback function will get called at the end of processing
    ## each vertical strip).
    ok <- .sink_chunking_is_compatible_with_vstrips(chunkdim(sink), sink_grid)
    if (ok) sink else NULL
}

### Returns a RealizationSink + its associated grid in a named list if we can
### take the "shared sink" route, or NULL if we can't.
make_shared_sink_and_grid_along_hstrips <-
    function(BACKEND, input_grid, sink_ncol,
             sink_rownames, sink_colnames, BPPARAM)
{
    sink_grid <- .make_sink_grid_of_hstrips(input_grid, sink_ncol)
    sink <- .make_shared_sink_along_hstrips(BACKEND, sink_grid,
                                            sink_rownames, sink_colnames,
                                            input_grid, BPPARAM)
    if (is.null(sink))
        return(NULL)
    list(sink=sink, sink_grid=sink_grid)
}

### Returns a RealizationSink + its associated grid in a named list if we can
### take the "shared sink" route, or NULL if we can't.
make_shared_sink_and_grid_along_vstrips <-
    function(BACKEND, input_grid, sink_nrow,
             sink_rownames, sink_colnames, BPPARAM)
{
    sink_grid <- .make_sink_grid_of_vstrips(input_grid, sink_nrow)
    sink <- .make_shared_sink_along_vstrips(BACKEND, sink_grid,
                                            sink_rownames, sink_colnames,
                                            input_grid, BPPARAM)
    if (is.null(sink))
        return(NULL)
    list(sink=sink, sink_grid=sink_grid)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### BLOCK_rowsum() and BLOCK_colsum()
###

### x: a matrix-like object (typically a DelayedMatrix).
### Walks on the grid defined on matrix-like object 'x'.
### If 'BACKEND' is NULL, returns an ordinary matrix. Otherwise, returns
### a DelayedMatrix object that is either pristine or the result of cbind'ing
### several pristine DelayedMatrix objects together (delayed cbind()).
BLOCK_rowsum <- function(x, group, reorder=TRUE, na.rm=FALSE,
                         grid=NULL, as.sparse=NA,
                         BACKEND=getAutoRealizationBackend(),
                         BPPARAM=getAutoBPPARAM(), verbose=NA)
{
    stopifnot(length(dim(x)) == 2L)  # matrix-like object
    verbose <- normarg_verbose(verbose)

    ugroup <- as.character(S4Arrays:::compute_ugroup(group, nrow(x), reorder))
    if (!isTRUEorFALSE(na.rm))
        stop(wmsg("'na.rm' must be TRUE or FALSE"))

    ## --- define INIT() ---

    ## INIT() must return a matrix of type "double" rather than "integer".
    ## This is to avoid integer overflows during the within-strip walks.
    INIT <- function(j, grid, ugroup, x_colnames) {
        vp <- grid[[1L, j]]
        dn <- list(ugroup, extractROWS(x_colnames, ranges(vp)[2L]))
        matrix(0.0, nrow=length(ugroup), ncol=ncol(vp), dimnames=dn)
    }
    INIT_MoreArgs <- list(ugroup=ugroup, x_colnames=colnames(x))

    ## --- define FUN() ---

    FUN <- function(init, block, group, ugroup, na.rm=FALSE) {
        if (is(block, "SparseArraySeed"))
            block <- as(block, "CsparseMatrix")  # to dgCMatrix or lgCMatrix
        vp <- currentViewport()
        group2 <- extractROWS(group, ranges(vp)[1L])
        block_ans <- rowsum(block, group2, reorder=FALSE, na.rm=na.rm)
        if (!is.matrix(block_ans))
            block_ans <- as.matrix(block_ans)
        m <- match(rownames(block_ans), ugroup)
        init[m, ] <- init[m, ] + block_ans
        init
    }
    FUN_MoreArgs <- list(group=group, ugroup=ugroup, na.rm=na.rm)

    ## --- define FINAL() ---

    if (is.null(BACKEND)) {
        FINAL <- if (verbose) final_vstrip_noop else NULL
        FINAL_MoreArgs <- list()
    } else {
        ## The "shared sink" route consists in using a single realization sink
        ## shared across all strips. Can we take this route?
        ## make_shared_sink_and_grid_along_vstrips() will figure it out and
        ## return a RealizationSink + its associated grid in a named list if
        ## it turns out that we can take the "shared sink" route, or NULL if
        ## we can't.
        grid <- best_grid_for_vstrip_apply(x, grid)
        sink_and_grid <- make_shared_sink_and_grid_along_vstrips(BACKEND,
                                                   grid, length(ugroup),
                                                   ugroup, colnames(x),
                                                   BPPARAM)
        if (is.null(sink_and_grid)) {
            FINAL <- function(init, j, grid, BACKEND, verbose) {
                realize_matrix(init, BACKEND, verbose)
            }
            FINAL_MoreArgs <- list(BACKEND=BACKEND, verbose=verbose)
        } else {
            ## "shared sink" route.
            FINAL <- function(init, j, grid, sink, sink_grid, verbose) {
                write_full_sink_cols(sink, sink_grid, j, init, verbose)
            }
            FINAL_MoreArgs <- c(sink_and_grid, list(verbose=verbose))
        }
    }

    ## --- block processing ---

    strip_results <- vstrip_apply(x, INIT, INIT_MoreArgs,
                                     FUN, FUN_MoreArgs,
                                     FINAL, FINAL_MoreArgs,
                                     grid=grid, as.sparse=as.sparse,
                                     BPPARAM=BPPARAM, verbose=verbose)

    ## --- turn output of block processing into object and return it ---

    if (is.null(BACKEND) || is.null(sink_and_grid)) {
        combine_strip_results("cbind", strip_results, verbose)
    } else {
        ## "shared sink" route.
        shared_sink_as_DelayedArray(sink_and_grid$sink, verbose)
    }
}

### x: a matrix-like object (typically a DelayedMatrix).
### Walks on the grid defined on matrix-like object 'x'.
### If 'BACKEND' is NULL, returns an ordinary matrix. Otherwise, returns
### a DelayedMatrix object that is either pristine or the result of rbind'ing
### several pristine DelayedMatrix objects together (delayed rbind()).
BLOCK_colsum <- function(x, group, reorder=TRUE, na.rm=FALSE,
                         grid=NULL, as.sparse=NA,
                         BACKEND=getAutoRealizationBackend(),
                         BPPARAM=getAutoBPPARAM(), verbose=NA)
{
    stopifnot(length(dim(x)) == 2L)  # matrix-like object
    verbose <- normarg_verbose(verbose)

    ugroup <- as.character(S4Arrays:::compute_ugroup(group, ncol(x), reorder))
    if (!isTRUEorFALSE(na.rm))
        stop(wmsg("'na.rm' must be TRUE or FALSE"))

    ## --- define INIT() ---

    ## INIT() must return a matrix of type "double" rather than "integer".
    ## This is to avoid integer overflows during the within-strip walks.
    INIT <- function(i, grid, ugroup, x_rownames) {
        vp <- grid[[i, 1L]]
        dn <- list(extractROWS(x_rownames, ranges(vp)[1L]), ugroup)
        matrix(0.0, nrow=nrow(vp), ncol=length(ugroup), dimnames=dn)
    }
    INIT_MoreArgs <- list(ugroup=ugroup, x_rownames=rownames(x))

    ## --- define FUN() ---

    FUN <- function(init, block, group, ugroup, na.rm=FALSE) {
        if (is(block, "SparseArraySeed"))
            block <- as(block, "CsparseMatrix")  # to dgCMatrix or lgCMatrix
        vp <- currentViewport()
        group2 <- extractROWS(group, ranges(vp)[2L])
        block_ans <- colsum(block, group2, reorder=FALSE, na.rm=na.rm)
        if (!is.matrix(block_ans))
            block_ans <- as.matrix(block_ans)
        m <- match(colnames(block_ans), ugroup)
        init[ , m] <- init[ , m] + block_ans
        init
    }
    FUN_MoreArgs <- list(group=group, ugroup=ugroup, na.rm=na.rm)

    ## --- define FINAL() ---

    if (is.null(BACKEND)) {
        FINAL <- if (verbose) final_hstrip_noop else NULL
        FINAL_MoreArgs <- list()
    } else {
        ## The "shared sink" route consists in using a single realization sink
        ## shared across all strips. Can we take this route?
        ## make_shared_sink_and_grid_along_hstrips() will figure it out and
        ## return a RealizationSink + its associated grid in a named list if
        ## it turns out that we can take the "shared sink" route, or NULL if
        ## we can't.
        grid <- best_grid_for_hstrip_apply(x, grid)
        sink_and_grid <- make_shared_sink_and_grid_along_hstrips(BACKEND,
                                                   grid, length(ugroup),
                                                   rownames(x), ugroup,
                                                   BPPARAM)
        if (is.null(sink_and_grid)) {
            FINAL <- function(init, i, grid, BACKEND, verbose) {
                realize_matrix(init, BACKEND, verbose)
            }
            FINAL_MoreArgs <- list(BACKEND=BACKEND, verbose=verbose)
        } else {
            ## "shared sink" route.
            FINAL <- function(init, i, grid, sink, sink_grid, verbose) {
                write_full_sink_rows(sink, sink_grid, i, init, verbose)
            }
            FINAL_MoreArgs <- c(sink_and_grid, list(verbose=verbose))
        }
    }

    ## --- block processing ---

    strip_results <- hstrip_apply(x, INIT, INIT_MoreArgs,
                                     FUN, FUN_MoreArgs,
                                     FINAL, FINAL_MoreArgs,
                                     grid=grid, as.sparse=as.sparse,
                                     BPPARAM=BPPARAM, verbose=verbose)

    ## --- turn output of block processing into object and return it ---

    if (is.null(BACKEND) || is.null(sink_and_grid)) {
        combine_strip_results("rbind", strip_results, verbose)
    } else {
        ## "shared sink" route.
        shared_sink_as_DelayedArray(sink_and_grid$sink, verbose)
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### rowsum() and colsum() methods
###

### S3/S4 combo for rowsum.DelayedMatrix
rowsum.DelayedMatrix <- function(x, group, reorder=TRUE, ...)
    BLOCK_rowsum(x, group, reorder=reorder, ...)
setMethod("rowsum", "DelayedMatrix", BLOCK_rowsum)

setMethod("colsum", "DelayedMatrix", BLOCK_colsum)

