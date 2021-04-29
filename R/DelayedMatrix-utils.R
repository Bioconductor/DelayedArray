### =========================================================================
### Common operations on DelayedMatrix objects
### -------------------------------------------------------------------------
###


.read_matrix_block <- function(...) {
    block <- read_block(..., as.sparse=NA)
    if (is(block, "SparseArraySeed"))
        block <- as(block, "CsparseMatrix")  # to dgCMatrix or lgCMatrix
    block 
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### rowsum() / colsum()
###

.compute_rowsum_for_block <- function(x, grid, i, j, group, na.rm=FALSE)
{
    viewport <- grid[[i, j]]
    block <- .read_matrix_block(x, viewport)
    group2 <- extractROWS(group, ranges(viewport)[1L])
    rowsum(block, group2, reorder=FALSE, na.rm=na.rm)
}
.compute_colsum_for_block <- function(x, grid, i, j, group, na.rm=FALSE)
{
    viewport <- grid[[i, j]]
    block <- .read_matrix_block(x, viewport)
    group2 <- extractROWS(group, ranges(viewport)[2L])
    colsum(block, group2, reorder=FALSE, na.rm=na.rm)
}

.compute_rowsum_for_grid_col <- function(x, grid, j, group, ugroup,
                                         na.rm=FALSE, verbose=FALSE)
{
    grid_nrow <- nrow(grid)
    grid_ncol <- ncol(grid)
    ans <- matrix(0L, nrow=length(ugroup), ncol=ncol(grid[[1L, j]]))
    ## Inner loop on the grid rows. Sequential.
    for (i in seq_len(grid_nrow)) {
        if (verbose)
            message("Processing block [[", i, "/", grid_nrow, ", ",
                                           j, "/", grid_ncol, "]] ... ",
                    appendLF=FALSE)
        block_ans <- .compute_rowsum_for_block(x, grid, i, j,
                                               group, na.rm=na.rm)
        m <- match(rownames(block_ans), ugroup)
        ans[m, ] <- ans[m, ] + block_ans
        if (verbose)
            message("OK")
    }
    ans
}
.compute_colsum_for_grid_row <- function(x, grid, i, group, ugroup,
                                         na.rm=FALSE, verbose=FALSE)
{
    grid_nrow <- nrow(grid)
    grid_ncol <- ncol(grid)
    ans <- matrix(0L, nrow=nrow(grid[[i, 1L]]), ncol=length(ugroup))
    ## Inner loop on the grid cols. Sequential.
    for (j in seq_len(grid_ncol)) {
        if (verbose)
            message("Processing block [[", i, "/", grid_nrow, ", ",
                                           j, "/", grid_ncol, "]] ... ",
                    appendLF=FALSE)
        block_ans <- .compute_colsum_for_block(x, grid, i, j,
                                               group, na.rm=na.rm)
        m <- match(colnames(block_ans), ugroup)
        ans[ , m] <- ans[ , m] + block_ans
        if (verbose)
            message("OK")
    }
    ans
}

.BLOCK_rowsum <- function(x, group, reorder=TRUE, na.rm=FALSE, grid=NULL)
{
    ugroup <- as.character(compute_ugroup(group, nrow(x), reorder))
    if (!isTRUEorFALSE(na.rm))
        stop(wmsg("'na.rm' must be TRUE or FALSE"))
    grid <- normarg_grid(grid, x)

    ## Outer loop on the grid columns. Parallelized.
    block_results <- bplapply2(seq_len(ncol(grid)),
        function(j) {
            .compute_rowsum_for_grid_col(x, grid, j, group, ugroup,
                                         na.rm=na.rm,
                                         verbose=get_verbose_block_processing())
        },
        BPPARAM=getAutoBPPARAM()
    )

    ans <- do.call(cbind, block_results)
    dimnames(ans) <- list(ugroup, colnames(x))
    ans
}
.BLOCK_colsum <- function(x, group, reorder=TRUE, na.rm=FALSE, grid=NULL)
{
    ugroup <- as.character(compute_ugroup(group, ncol(x), reorder))
    if (!isTRUEorFALSE(na.rm))
        stop(wmsg("'na.rm' must be TRUE or FALSE"))
    grid <- normarg_grid(grid, x)

    ## Outer loop on the grid rows. Parallelized.
    block_results <- bplapply2(seq_len(nrow(grid)),
        function(i) {
            .compute_colsum_for_grid_row(x, grid, i, group, ugroup,
                                         na.rm=na.rm,
                                         verbose=get_verbose_block_processing())
        },
        BPPARAM=getAutoBPPARAM()
    )

    ans <- do.call(rbind, block_results)
    dimnames(ans) <- list(rownames(x), ugroup)
    ans
}

### S3/S4 combo for rowsum.DelayedMatrix
rowsum.DelayedMatrix <- function(x, group, reorder=TRUE, ...)
    .BLOCK_rowsum(x, group, reorder=reorder, ...)
setMethod("rowsum", "DelayedMatrix", .BLOCK_rowsum)

setMethod("colsum", "DelayedMatrix", .BLOCK_colsum)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Matrix multiplication
###
### We only support multiplication of an ordinary matrix (typically
### small) by a DelayedMatrix object (typically big). Multiplication of 2
### DelayedMatrix objects is not supported.
###

.BLOCK_mult_by_left_matrix <- function(x, y)
{
    stopifnot(is.matrix(x),
              is(y, "DelayedMatrix") || is.matrix(y),
              ncol(x) == nrow(y))

    ans_dim <- c(nrow(x), ncol(y))
    ans_dimnames <- simplify_NULL_dimnames(list(rownames(x), colnames(y)))
    ans_type <- typeof(vector(type(x), 1L) * vector(type(y), 1L))
    sink <- AutoRealizationSink(ans_dim, ans_dimnames, ans_type)
    on.exit(close(sink))

    y_grid <- colAutoGrid(y)
    ans_spacings <- c(ans_dim[[1L]], ncol(y_grid[[1L]]))
    ans_grid <- RegularArrayGrid(ans_dim, ans_spacings)  # parallel to 'y_grid'
    nblock <- length(y_grid)  # same as 'length(ans_grid)'
    for (bid in seq_len(nblock)) {
        if (get_verbose_block_processing())
            message("Processing block ", bid, "/", nblock, " ... ",
                    appendLF=FALSE)
        y_viewport <- y_grid[[bid]]
        block <- .read_matrix_block(y, y_viewport)
        block_ans <- x %*% block
        write_block(sink, ans_grid[[bid]], block_ans)
        if (get_verbose_block_processing())
            message("OK")
    }
    as(sink, "DelayedArray")
}

setMethod("%*%", c("ANY", "DelayedMatrix"),
    function(x, y)
    {
        if (!is.matrix(x)) {
            if (!is.vector(x))
                stop(wmsg("matrix multiplication of a ", class(x), " object ",
                          "by a ", class(y), " object is not supported"))
            x_len <- length(x)
            y_nrow <- nrow(y)
            if (x_len != 0L && x_len != y_nrow)
                stop(wmsg("non-conformable arguments"))
            x <- matrix(x, ncol=y_nrow)
        }
        .BLOCK_mult_by_left_matrix(x, y)
    }
)

setMethod("%*%", c("DelayedMatrix", "ANY"),
    function(x, y)
    {
        if (!is.matrix(y)) {
            if (!is.vector(y))
                stop(wmsg("matrix multiplication of a ", class(x), " object ",
                          "by a ", class(y), " object is not supported"))
            y_len <- length(y)
            x_ncol <- ncol(x)
            if (y_len != 0L && y_len != x_ncol)
                stop(wmsg("non-conformable arguments"))
            y <- matrix(y, nrow=x_ncol)
        }
        t(t(y) %*% t(x))
    }
)

.BLOCK_matrix_mult <- function(x, y)
{
    stop(wmsg("Matrix multiplication of 2 DelayedMatrix derivatives is not ",
              "supported at the moment. Only matrix multiplication between ",
              "a DelayedMatrix derivative and an ordinary matrix or vector ",
              "is supported for now."))

    x_dim <- dim(x)
    y_dim <- dim(y)
    stopifnot(length(x_dim) == 2L, length(y_dim) == 2L, ncol(x) == nrow(y))

    ans_dim <- c(nrow(x), ncol(y))
    ans_dimnames <- simplify_NULL_dimnames(list(rownames(x), colnames(y)))
    ans_type <- typeof(vector(type(x), 1L) * vector(type(y), 1L))
    sink <- AutoRealizationSink(ans_dim, ans_dimnames, ans_type)
    on.exit(close(sink))
}

setMethod("%*%", c("DelayedMatrix", "DelayedMatrix"), .BLOCK_matrix_mult)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Parallelized schemes for matrix multiplication.
###
### by Aaron Lun
###
### This splits one or both matrices into blocks according to the
### desired parallelization scheme, and distributes them to workers.
### This also requires care to respect the maximum block size.
###

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
        out <- bplapply2(seq_len(length(x_grid)),
                         FUN=.left_mult, 
                         x=x, y=y, grid=x_grid, 
                         MULT=MULT, 
                         BPPARAM=BPPARAM)
        ans <- do.call(rbind, out)
    } else if (chosen_scheme=="y") {
        out <- bplapply2(seq_len(length(y_grid)),
                         FUN=.right_mult, 
                         x=x, y=y, grid=y_grid, 
                         MULT=MULT, 
                         BPPARAM=BPPARAM)
        ans <- do.call(cbind, out)
    }

    realize(ans)
}

setMethod("%*%", c("DelayedMatrix", "ANY"), function(x, y) {
    if (is.null(dim(y))) y <- cbind(y)
    .super_BLOCK_mult(x, y, MULT=`%*%`)
})

setMethod("%*%", c("ANY", "DelayedMatrix"), function(x, y) {
    if (is.null(dim(x))) x <- rbind(x)
    .super_BLOCK_mult(x, y, MULT=`%*%`)
})

setMethod("%*%", c("DelayedMatrix", "DelayedMatrix"), function(x, y) .super_BLOCK_mult(x, y, MULT=`%*%`))

setMethod("crossprod", c("DelayedMatrix", "ANY"), function(x, y) {
    if (is.null(dim(y))) y <- cbind(y)
    .super_BLOCK_mult(x, y, MULT=crossprod, transposed.x=TRUE)
})

setMethod("crossprod", c("ANY", "DelayedMatrix"), function(x, y) {
    if (is.null(dim(x))) x <- rbind(x)
    .super_BLOCK_mult(x, y, MULT=crossprod, transposed.x=TRUE)
})

setMethod("crossprod", c("DelayedMatrix", "DelayedMatrix"), function(x, y) 
    .super_BLOCK_mult(x, y, MULT=crossprod, transposed.x=TRUE)
)

# tcrossprod with vector 'y' doesn't work in base, and so it won't work here either.
setMethod("tcrossprod", c("DelayedMatrix", "ANY"), function(x, y) 
    .super_BLOCK_mult(x, y, MULT=tcrossprod, transposed.y=TRUE)
)

setMethod("tcrossprod", c("ANY", "DelayedMatrix"), function(x, y) {
    if (is.null(dim(x))) x <- rbind(x)
    .super_BLOCK_mult(x, y, MULT=tcrossprod, transposed.y=TRUE)
})

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
        out <- bplapply2(seq_len(length(slow)),
                         FUN=.left_mult, 
                         x=x, y=x, grid=slow,
                         MULT=MULT,
                         BPPARAM=BPPARAM)
        ans <- do.call(rbind, out)

    } else {
        ans <- bplapply2(seq_len(length(fast)),
                         FUN=.solo_mult,
                         x=x, grid=fast,
                         MULT=MULT,
                         BPPARAM=BPPARAM)
        ans <- Reduce("+", ans)
    }

    realize(ans)
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
    set_user_option("auto.mult.parallel.agnostic", agnostic)
}

getAutoMultParallelAgnostic <- function() {
    get_user_option("auto.mult.parallel.agnostic")
}
