### =========================================================================
### DelayedMatrix %*%, crossprod(), and tcrossprod()
### -------------------------------------------------------------------------
###
### The %*%, crossprod(), and tcrossprod() methods for DelayedMatrix objects
### are block processed.
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Helpers for BLOCK_mult_Lgrid() and BLOCK_mult_Rgrid()
###

.make_shared_sink_and_grid_for_Lgrid_apply <-
    function(x, y, Lgrid, BACKEND, BPPARAM,
                   transpose.x=FALSE, transpose.y=FALSE)
{
    if (transpose.x) {
        input_grid <- t(best_grid_for_vstrip_apply(x, Lgrid))
        sink_rownames <- colnames(x)
    } else {
        input_grid <- best_grid_for_hstrip_apply(x, Lgrid)
        sink_rownames <- rownames(x)
    }
    if (transpose.y) {
        sink_ncol <- nrow(y)
        sink_colnames <- rownames(y)
    } else {
        sink_ncol <- ncol(y)
        sink_colnames <- colnames(y)
    }
    make_shared_sink_and_grid_along_hstrips(BACKEND,
                                            input_grid, sink_ncol,
                                            sink_rownames, sink_colnames,
                                            BPPARAM)
}

.make_shared_sink_and_grid_for_Rgrid_apply <-
    function(x, y, Rgrid, BACKEND, BPPARAM,
                   transpose.x=FALSE, transpose.y=FALSE)
{
    if (transpose.y) {
        input_grid <- t(best_grid_for_hstrip_apply(y, Rgrid))
        sink_colnames <- rownames(y)
    } else {
        input_grid <- best_grid_for_vstrip_apply(y, Rgrid)
        sink_colnames <- colnames(y)
    }
    if (transpose.x) {
        sink_nrow <- ncol(x)
        sink_rownames <- colnames(x)
    } else {
        sink_nrow <- nrow(x)
        sink_rownames <- rownames(x)
    }
    make_shared_sink_and_grid_along_vstrips(BACKEND,
                                            input_grid, sink_nrow,
                                            sink_rownames, sink_colnames,
                                            BPPARAM)
}

### INIT, BLOCK_OP: callback functions.
### INIT() must take 3 arguments: i (or j), grid, y.
### BLOCK_OP() must take 3 arguments: x_block, y, vp_ranges.
### See BLOCK_mult_Lgrid() below for the other arguments.
### Walks on the "left grid" which is defined on matrix-like object 'x'.
.Lgrid_apply <- function(x, y, Lgrid, as.sparse, BACKEND, BPPARAM, verbose,
                               INIT, BLOCK_OP,
                               transpose.x=FALSE, transpose.y=FALSE)
{
    verbose <- normarg_verbose(verbose)

    ## --- define FUN() ---

    FUN <- function(init, block, y, BLOCK_OP) {
        if (is(block, "SparseArraySeed"))
            block <- as(block, "CsparseMatrix")  # to dgCMatrix or lgCMatrix
        vp <- currentViewport()
        block_ans <- BLOCK_OP(block, y, ranges(vp))
        if (!is.matrix(block_ans))
            block_ans <- as.matrix(block_ans)
        init + block_ans
    }
    FUN_MoreArgs <- list(y=y, BLOCK_OP=BLOCK_OP)

    ## --- define FINAL() ---

    if (is.null(BACKEND)) {
        if (verbose) {
            FINAL <- if (transpose.x) final_vstrip_noop else final_hstrip_noop
        } else {
            FINAL <- NULL
        }
        FINAL_MoreArgs <- list()
    } else {
        ## The "shared sink" route consists in using a single realization sink
        ## shared across all strips. Can we take this route?
        ## .make_shared_sink_and_grid_for_Lgrid_apply() will figure it out and
        ## return a RealizationSink + its associated grid in a named list if
        ## it turns out that we can take the "shared sink" route, or NULL if
        ## we can't.
        sink_and_grid <- .make_shared_sink_and_grid_for_Lgrid_apply(x, y,
                                                    Lgrid, BACKEND, BPPARAM,
                                                    transpose.x, transpose.y)
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

    STRIP_APPLY <- if (transpose.x) vstrip_apply else hstrip_apply
    INIT_MoreArgs <- list(y=y)
    strip_results <- STRIP_APPLY(x, INIT, INIT_MoreArgs,
                                    FUN, FUN_MoreArgs,
                                    FINAL, FINAL_MoreArgs,
                                    grid=Lgrid, as.sparse=as.sparse,
                                    BPPARAM=BPPARAM, verbose=verbose)

    ## --- turn output of block processing into object and return it ---

    if (is.null(BACKEND) || is.null(sink_and_grid)) {
        combine_strip_results("rbind", strip_results, verbose)
    } else {
        ## "shared sink" route.
        shared_sink_as_DelayedArray(sink_and_grid$sink, verbose)
    }
}

### INIT and BLOCK_OP: callback functions.
### INIT() must take 3 arguments: j (or i), grid, x.
### BLOCK_OP() must take 3 arguments: x, y_block, vp_ranges.
### See BLOCK_mult_Rgrid() below for the other arguments.
### Walks on the "right grid" which is defined on matrix-like object 'y'.
.Rgrid_apply <- function(x, y, Rgrid, as.sparse, BACKEND, BPPARAM, verbose,
                               INIT, BLOCK_OP,
                               transpose.x=FALSE, transpose.y=FALSE)
{
    verbose <- normarg_verbose(verbose)

    ## --- define FUN() ---

    FUN <- function(init, block, x, BLOCK_OP) {
        if (is(block, "SparseArraySeed"))
            block <- as(block, "CsparseMatrix")  # to dgCMatrix or lgCMatrix
        vp <- currentViewport()
        block_ans <- BLOCK_OP(x, block, ranges(vp))
        if (!is.matrix(block_ans))
            block_ans <- as.matrix(block_ans)
        init + block_ans
    }
    FUN_MoreArgs <- list(x=x, BLOCK_OP=BLOCK_OP)

    ## --- define FINAL() ---

    if (is.null(BACKEND)) {
        if (verbose) {
            FINAL <- if (transpose.y) final_hstrip_noop else final_vstrip_noop
        } else {
            FINAL <- NULL
        }
        FINAL_MoreArgs <- list()
    } else {
        ## The "shared sink" route consists in using a single realization sink
        ## shared across all strips. Can we take this route?
        ## .make_shared_sink_and_grid_for_Rgrid_apply() will figure it out and
        ## return a RealizationSink + its associated grid in a named list if
        ## it turns out that we can take the "shared sink" route, or NULL if
        ## we can't.
        sink_and_grid <- .make_shared_sink_and_grid_for_Rgrid_apply(x, y,
                                                    Rgrid, BACKEND, BPPARAM,
                                                    transpose.x, transpose.y)
        if (is.null(sink_and_grid)) {
            FINAL <- function(init, j, grid, BACKEND, verbose) {
                realize_matrix(init, BACKEND=BACKEND, verbose)
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

    STRIP_APPLY <- if (transpose.y) hstrip_apply else vstrip_apply
    INIT_MoreArgs <- list(x=x)
    strip_results <- STRIP_APPLY(y, INIT, INIT_MoreArgs,
                                    FUN, FUN_MoreArgs,
                                    FINAL, FINAL_MoreArgs,
                                    grid=Rgrid, as.sparse=as.sparse,
                                    BPPARAM=BPPARAM, verbose=verbose)

    ## --- turn output of block processing into object and return it ---

    if (is.null(BACKEND) || is.null(sink_and_grid)) {
        combine_strip_results("cbind", strip_results, verbose)
    } else {
        ## "shared sink" route.
        shared_sink_as_DelayedArray(sink_and_grid$sink, verbose)
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### BLOCK_mult_Lgrid() and BLOCK_mult_Rgrid()
###
### These are the 2 workhorses behind block matrix multiplication between:
###   - **any** matrix-like object (typically a DelayedMatrix),
###   - an ordinary matrix (or other supported matrix-like object,
###     see .is_supported() above).
###
### Should be able to handle any type() supported by base::`%*%`, that is,
### integer, double, and complex. However, the realization backend specified
### via `BACKEND` might introduce some restrictions e.g. will it support
### realization of a matrix of type complex?
###

### Supported matrix-like objects must support [ as well as native %*%,
### crossprod(), and tcrossprod() with the blocks returned by read_block()
### (i.e. which are either ordinary matrices or SparseMatrix objects).
.is_supported <- function(x)
    is.matrix(x) || is(x, "sparseMatrix") || is(x, "SparseMatrix")

### x: a matrix-like object (typically a DelayedMatrix).
### y: an ordinary matrix or other supported object (see .is_supported above).
### If 'BACKEND' is NULL, returns an ordinary matrix. Otherwise, returns
### a DelayedMatrix object that is either pristine or the result of rbind'ing
### several pristine DelayedMatrix objects together (delayed rbind()).
BLOCK_mult_Lgrid <- function(x, y, Lgrid=NULL, as.sparse=NA,
                                   BACKEND=getAutoRealizationBackend(),
                                   BPPARAM=getAutoBPPARAM(), verbose=NA,
                                   op=c("mult", "crossprod", "tcrossprod"))
{
    ## See the BLOCK_OP() callback functions below for what operations will
    ## effectively be performed on 'y'.
    if (!.is_supported(y))
        stop(wmsg("this operation does not support 'y' ",
                  "of class ", class(y)[[1L]]))
    stopifnot(length(dim(x)) == 2L)
    op <- match.arg(op)
    transpose.x <- transpose.y <- FALSE

    ## All INIT() callback functions must return a matrix of type "double"
    ## rather than "integer". This is to avoid integer overflows during the
    ## within-strip walks.
    switch(op,
        mult={
            stopifnot(ncol(x) == nrow(y))
            INIT <- function(i, grid, y) {
                matrix(0.0, nrow=nrow(grid[[i, 1L]]), ncol=ncol(y))
            }
            BLOCK_OP <- function(x_block, y, vp_ranges) {
                idx <- (start(vp_ranges)[[2L]]):(end(vp_ranges)[[2L]])
                base::`%*%`(x_block, y[idx, , drop=FALSE])
            }
        },
        crossprod={
            transpose.x <- TRUE
            stopifnot(nrow(x) == nrow(y))
            INIT <- function(j, grid, y) {
                matrix(0.0, nrow=ncol(grid[[1L, j]]), ncol=ncol(y))
            }
            BLOCK_OP <- function(x_block, y, vp_ranges) {
                idx <- (start(vp_ranges)[[1L]]):(end(vp_ranges)[[1L]])
                base::crossprod(x_block, y[idx, , drop=FALSE])
            }
        },
        tcrossprod={
            transpose.y <- TRUE
            stopifnot(ncol(x) == ncol(y))
            INIT <- function(i, grid, y) {
                matrix(0.0, nrow=nrow(grid[[i, 1L]]), ncol=nrow(y))
            }
            BLOCK_OP <- function(x_block, y, vp_ranges) {
                idx <- (start(vp_ranges)[[2L]]):(end(vp_ranges)[[2L]])
                base::tcrossprod(x_block, y[ , idx, drop=FALSE])
            }
        },
        stop(wmsg("invalid 'op'"))  # should never happen
    )

    .Lgrid_apply(x, y, Lgrid, as.sparse, BACKEND, BPPARAM, verbose,
                       INIT, BLOCK_OP, transpose.x, transpose.y)
}

### x: an ordinary matrix or other supported object (see .is_supported above).
### y: a matrix-like object (typically a DelayedMatrix).
### If 'BACKEND' is NULL, returns an ordinary matrix. Otherwise, returns
### a DelayedMatrix object that is either pristine or the result of cbind'ing
### several pristine DelayedMatrix objects together (delayed cbind()).
BLOCK_mult_Rgrid <- function(x, y, Rgrid=NULL, as.sparse=NA,
                                   BACKEND=getAutoRealizationBackend(),
                                   BPPARAM=getAutoBPPARAM(), verbose=NA,
                                   op=c("mult", "crossprod", "tcrossprod"))
{
    ## See the BLOCK_OP() callback functions below for what operations will
    ## effectively be performed on 'x'.
    if (!.is_supported(x))
        stop(wmsg("this operation does not support 'x' ",
                  "of class ", class(x)[[1L]]))
    stopifnot(length(dim(y)) == 2L)
    op <- match.arg(op)
    transpose.x <- transpose.y <- FALSE

    ## All INIT() callback functions must return a matrix of type "double"
    ## rather than "integer". This is to avoid integer overflows during the
    ## within-strip walks.
    switch(op,
        mult={
            stopifnot(ncol(x) == nrow(y))
            INIT <- function(j, grid, x) {
                matrix(0.0, nrow=nrow(x), ncol=ncol(grid[[1L, j]]))
            }
            BLOCK_OP <- function(x, y_block, vp_ranges) {
                idx <- (start(vp_ranges)[[1L]]):(end(vp_ranges)[[1L]])
                base::`%*%`(x[ , idx, drop=FALSE], y_block)
            }
        },
        crossprod={
            transpose.x <- TRUE
            stopifnot(nrow(x) == nrow(y))
            INIT <- function(j, grid, x) {
                matrix(0.0, nrow=ncol(x), ncol=ncol(grid[[1L, j]]))
            }
            BLOCK_OP <- function(x, y_block, vp_ranges) {
                idx <- (start(vp_ranges)[[1L]]):(end(vp_ranges)[[1L]])
                base::crossprod(x[idx, , drop=FALSE], y_block)
            }
        },
        tcrossprod={
            transpose.y <- TRUE
            stopifnot(ncol(x) == ncol(y))
            INIT <- function(i, grid, x) {
                matrix(0.0, nrow=nrow(x), ncol=nrow(grid[[i, 1L]]))
            }
            BLOCK_OP <- function(x, y_block, vp_ranges) {
                idx <- (start(vp_ranges)[[2L]]):(end(vp_ranges)[[2L]])
                base::tcrossprod(x[ , idx, drop=FALSE], y_block)
            }
        },
        stop(wmsg("invalid 'op'"))  # should never happen
    )

    .Rgrid_apply(x, y, Rgrid, as.sparse, BACKEND, BPPARAM, verbose,
                       INIT, BLOCK_OP, transpose.x, transpose.y)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### %*%, crossprod(), and tcrossprod() methods between:
###   - a DelayedMatrix object,
###   - an ordinary matrix or vector (or other supported matrix-like
###     object, see .is_supported() above)
###

setMethod("%*%", c("DelayedMatrix", "ANY"),
    function(x, y)
    {
        if (is.atomic(y) && is.vector(y)) {
            ## Returns a 1-col ordinary matrix (like base::`%*%` does).
            y <- cbind(y, deparse.level=0)
            BLOCK_mult_Lgrid(x, y, BACKEND=NULL)
        } else {
            BLOCK_mult_Lgrid(x, y)
        }
    }
)
setMethod("%*%", c("ANY", "DelayedMatrix"),
    function(x, y)
    {
        if (is.atomic(x) && is.vector(x)) {
            ## Returns a 1-row ordinary matrix (like base::`%*%` does).
            x <- rbind(x, deparse.level=0)
            BLOCK_mult_Rgrid(x, y, BACKEND=NULL)
        } else {
            BLOCK_mult_Rgrid(x, y)
        }
    }
)

setMethod("crossprod", c("DelayedMatrix", "ANY"),
    function(x, y)
    {
        if (is.atomic(y) && is.vector(y)) {
            ## Returns a 1-col ordinary matrix (like base::crossprod() does).
            y <- cbind(y, deparse.level=0)
            BLOCK_mult_Lgrid(x, y, BACKEND=NULL, op="crossprod")
        } else {
            BLOCK_mult_Lgrid(x, y, op="crossprod")
        }
    }
)
setMethod("crossprod", c("ANY", "DelayedMatrix"),
    function(x, y)
    {
        if (is.atomic(x) && is.vector(x)) {
            ## Returns a 1-row ordinary matrix (like base::crossprod() does).
            x <- cbind(x, deparse.level=0)
            BLOCK_mult_Rgrid(x, y, BACKEND=NULL, op="crossprod")
        } else {
            BLOCK_mult_Rgrid(x, y, op="crossprod")
        }
    }
)

setMethod("tcrossprod", c("DelayedMatrix", "ANY"),
    function(x, y)
    {
        if (is.atomic(y) && is.vector(y)) {
            ## Note that base::tcrossprod() does not work with a vector on
            ## the right!
            ## Returns a 1-col ordinary matrix (like base::tcrossprod() would
            ## probably do if it were supporting a vector on the right).
            y <- rbind(y, deparse.level=0)
            BLOCK_mult_Lgrid(x, y, BACKEND=NULL, op="tcrossprod")
        } else {
            BLOCK_mult_Lgrid(x, y, op="tcrossprod")
        }
    }
)
setMethod("tcrossprod", c("ANY", "DelayedMatrix"),
    function(x, y)
    {
        if (is.atomic(x) && is.vector(x)) {
            ## Returns a 1-row ordinary matrix (like base::tcrossprod() does).
            x <- rbind(x, deparse.level=0)
            BLOCK_mult_Rgrid(x, y, BACKEND=NULL, op="tcrossprod")
        } else {
            BLOCK_mult_Rgrid(x, y, op="tcrossprod")
        }
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Parallelized schemes for matrix multiplication.
###
### by Aaron Lun
###
### This splits one or both matrices into blocks according to the
### desired parallelization scheme, and distributes them to workers.
### This also requires care to respect the maximum block size.
###

### A thin wrapper around read_block().
### TODO: Get rid of this when we switch from SparseArraySeed to
### SVT_SparseArray (then callers should just call read_block() directly).
.read_matrix_block <- function(...) {
    block <- read_block(..., as.sparse=NA)
    if (is(block, "SparseArraySeed"))
        block <- as(block, "CsparseMatrix")  # to dgCMatrix or lgCMatrix
    block
}

.grid_by_dimension <- function(x, nworkers)
# Splits a dimension of the matrix into at least 'nworkers' blocks.
# If the block size is too large, it is reduced to obtain the desired
# number of blocks in order for parallelization to be effective.
{
    old <- getAutoBlockLength(type(x))

    ideal_size_by_row <- max(1, ceiling(nrow(x)/nworkers) * ncol(x))
    if (old > ideal_size_by_row) {
        row_grid <- rowAutoGrid(x, block.length=ideal_size_by_row)
    } else {
        row_grid <- rowAutoGrid(x)
    }

    ideal_size_by_col <- max(1, ceiling(ncol(x)/nworkers) * nrow(x))
    if (old > ideal_size_by_col) {
        col_grid <- colAutoGrid(x, block.length=ideal_size_by_col)
    } else {
        col_grid <- colAutoGrid(x)
    }

    list(row=row_grid, col=col_grid)
}

.left_mult <- function(bid, grid, x, y, MULT) {
    # this, and all other calls, had better yield a non-DA, otherwise MULT will recurse endlessly.
    block <- .read_matrix_block(x, grid[[bid]])
    MULT(block, y)
}

.right_mult <- function(bid, grid, x, y, MULT) {
    block <- .read_matrix_block(y, grid[[bid]])
    MULT(x, block)
}

.super_BLOCK_mult <- function(x, y, MULT, transposed.x=FALSE, transposed.y=FALSE, BPPARAM=getAutoBPPARAM())
# Controller function that split jobs for a multiplication function "MULT".
# This accommodates %*%, crossprod and tcrossprod for two arguments.
{
    if (is.null(BPPARAM)) {
        nworkers <- 1L
    } else {
        nworkers <- BiocParallel::bpnworkers(BPPARAM)
    }

    # Choosing the right dimension to iterate over, depending on MULT.
    x_grid <- .grid_by_dimension(x, nworkers)
    if (transposed.x) {
        x_grid <- x_grid$col
    } else {
        x_grid <- x_grid$row
    }

    y_grid <- .grid_by_dimension(y, nworkers)
    if (transposed.y) {
        y_grid <- y_grid$row
    } else {
        y_grid <- y_grid$col
    }

    # Always iterating over the 'larger' matrix, to better split up the work.
    # In the context of file-backed matrices, this operates under the heuristic
    # that the larger matrix is the file-backed one.
    if (length(x) > length(y)) {
        chosen_scheme <- "x"
    } else {
        chosen_scheme <- "y"
    }

    # Switch to iteration over the other argument if the chosen one is
    # single-block and non-DA (at which point you might as well iterate
    # over the other argument anyway). This avoids infinite recursion
    # when 'x' or 'y' fail to get realized via read_block().
    if (chosen_scheme=="x" && length(x_grid)==1L && !is(x, "DelayedMatrix")) {
        chosen_scheme <- "y"
    } else if (chosen_scheme=="y" && length(y_grid)==1L && !is(y, "DelayedMatrix")) {
        chosen_scheme <- "x"
    }

    old <- getAutoBPPARAM()
    on.exit(setAutoBPPARAM(old))
    setAutoBPPARAM(NULL) # Avoid re-parallelizing in further calls to 'MULT'.

    if (chosen_scheme=="x") {
        out <- S4Arrays:::bplapply2(seq_len(length(x_grid)),
                         FUN=.left_mult,
                         x=x, y=y, grid=x_grid,
                         MULT=MULT,
                         BPPARAM=BPPARAM)
        ans <- do.call(rbind, out)
    } else if (chosen_scheme=="y") {
        out <- S4Arrays:::bplapply2(seq_len(length(y_grid)),
                         FUN=.right_mult,
                         x=x, y=y, grid=y_grid,
                         MULT=MULT,
                         BPPARAM=BPPARAM)
        ans <- do.call(cbind, out)
    }

    realize(ans)
}

setMethod("%*%", c("DelayedMatrix", "DelayedMatrix"), function(x, y) .super_BLOCK_mult(x, y, MULT=`%*%`))

setMethod("crossprod", c("DelayedMatrix", "DelayedMatrix"), function(x, y)
    .super_BLOCK_mult(x, y, MULT=crossprod, transposed.x=TRUE)
)

setMethod("tcrossprod", c("DelayedMatrix", "DelayedMatrix"), function(x, y)
    .super_BLOCK_mult(x, y, MULT=tcrossprod, transposed.y=TRUE)
)

.solo_mult <- function(bid, grid, x, MULT) {
    block <- read_block(x, grid[[bid]])
    MULT(block)
}

.super_BLOCK_self <- function(x, MULT, transposed=FALSE, BPPARAM=getAutoBPPARAM())
# Controller function that split jobs for a multiplication function "MULT".
# This accommodates crossprod and tcrossprod for single arguments.
{
    if (is.null(BPPARAM)) {
        nworkers <- 1L
    } else {
        nworkers <- BiocParallel::bpnworkers(BPPARAM)
    }

    # Choosing the right dimension to iterate over, depending on MULT.
    grid <- .grid_by_dimension(x, nworkers)
    if (transposed) {
        fast <- grid$col
        slow <- grid$row
    } else {
        fast <- grid$row
        slow <- grid$col
    }

    old <- getAutoBPPARAM()
    on.exit(setAutoBPPARAM(old))
    setAutoBPPARAM(NULL) # Avoid re-parallelizing in further calls to 'MULT'.

    if (getAutoMultParallelAgnostic()) {
        out <- S4Arrays:::bplapply2(seq_len(length(slow)),
                         FUN=.left_mult,
                         x=x, y=x, grid=slow,
                         MULT=MULT,
                         BPPARAM=BPPARAM)
        ans <- do.call(rbind, out)

    } else {
        ans <- S4Arrays:::bplapply2(seq_len(length(fast)),
                         FUN=.solo_mult,
                         x=x, grid=fast,
                         MULT=MULT,
                         BPPARAM=BPPARAM)
        ans <- Reduce("+", ans)
    }

    DelayedArray(realize(ans))
}

setMethod("crossprod", c("DelayedMatrix", "missing"), function(x, y)
    .super_BLOCK_self(x, MULT=crossprod)
)

setMethod("tcrossprod", c("DelayedMatrix", "missing"), function(x, y)
    .super_BLOCK_self(x, MULT=tcrossprod, transposed=TRUE)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### User-visible global settings for parallelized matrix multiplication.
###
### by Aaron Lun
###
### This allows the user to specify whether or not they want to guarantee
### the identical matrix products regardless of the number of workers.
### This is because splitting by the common dimension does not preserve the
### order of addition operations, which changes the output due to numerical
### imprecision in the inner products of each vector.
###

setAutoMultParallelAgnostic <- function(agnostic=TRUE) {
    S4Arrays:::set_user_option("auto.mult.parallel.agnostic", agnostic)
}

getAutoMultParallelAgnostic <- function() {
    S4Arrays:::get_user_option("auto.mult.parallel.agnostic")
}

