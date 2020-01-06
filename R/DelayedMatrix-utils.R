### =========================================================================
### Common operations on DelayedMatrix objects
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### rowsum() / colsum()
###

colsum <- function(x, group, reorder=TRUE, na.rm=FALSE, ...)
{
    t(rowsum(t(x), group, reorder=reorder, na.rm=na.rm, ...))
}

.compute_rowsum_for_block <- function(x, grid, i, j, group, na.rm=FALSE)
{
    viewport <- grid[[i, j]]
    block <- read_block(x, viewport)
    group2 <- extractROWS(group, ranges(viewport)[1L])
    rowsum(block, group2, reorder=FALSE, na.rm=na.rm)
}
.compute_colsum_for_block <- function(x, grid, i, j, group, na.rm=FALSE)
{
    viewport <- grid[[i, j]]
    block <- read_block(x, viewport)
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

.compute_ugroup <- function(group, expected_group_len, reorder=TRUE)
{
    if (!(is.vector(group) || is.factor(group)))
        stop(wmsg("'group' must be a vector or factor"))
    if (length(group) != expected_group_len)
        stop(wmsg("incorrect length for 'group'"))
    if (!isTRUEorFALSE(reorder))
        stop(wmsg("'reorder' must be TRUE or FALSE"))
    ## Taken from base::rowsum.default().
    ugroup <- unique(group)
    if (anyNA(ugroup))
        warning(wmsg("missing values for 'group'"))
    if (reorder)
        ugroup <- sort(ugroup, na.last=TRUE, method="quick")
    as.character(ugroup)
}

.BLOCK_rowsum <- function(x, group, reorder=TRUE, na.rm=FALSE, grid=NULL)
{
    ugroup <- .compute_ugroup(group, nrow(x), reorder)
    if (!isTRUEorFALSE(na.rm))
        stop(wmsg("'na.rm' must be TRUE or FALSE"))
    grid <- normarg_grid(grid, x)

    ## Outer loop on the grid columns. Parallelized.
    block_results <- bplapply(seq_len(ncol(grid)),
        function(j) {
            .compute_rowsum_for_grid_col(x, grid, j, group, ugroup,
                                         na.rm=na.rm,
                                         verbose=get_verbose_block_processing())
        },
        BPPARAM=getAutoBPPARAM()
    )

    ans <- do.call(cbind, block_results)
    dimnames(ans) <- list(as.character(ugroup), colnames(x))
    ans
}
.BLOCK_colsum <- function(x, group, reorder=TRUE, na.rm=FALSE, grid=NULL)
{
    ugroup <- .compute_ugroup(group, ncol(x), reorder)
    if (!isTRUEorFALSE(na.rm))
        stop(wmsg("'na.rm' must be TRUE or FALSE"))
    grid <- normarg_grid(grid, x)

    ## Outer loop on the grid rows. Parallelized.
    block_results <- bplapply(seq_len(nrow(grid)),
        function(i) {
            .compute_colsum_for_grid_row(x, grid, i, group, ugroup,
                                         na.rm=na.rm,
                                         verbose=get_verbose_block_processing())
        },
        BPPARAM=getAutoBPPARAM()
    )

    ans <- do.call(rbind, block_results)
    dimnames(ans) <- list(rownames(x), as.character(ugroup))
    ans
}

setGeneric("rowsum", signature="x")
setGeneric("colsum", signature="x")

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
    sink <- RealizationSink(ans_dim, ans_dimnames, ans_type)
    on.exit(close(sink))

    y_grid <- colGrid(y)
    ans_spacings <- c(ans_dim[[1L]], ncol(y_grid[[1L]]))
    ans_grid <- RegularArrayGrid(ans_dim, ans_spacings)  # parallel to 'y_grid'
    nblock <- length(y_grid)  # same as 'length(ans_grid)'
    for (b in seq_len(nblock)) {
        if (get_verbose_block_processing())
            message("Processing block ", b, "/", nblock, " ... ",
                    appendLF=FALSE)
        y_viewport <- y_grid[[b]]
        block <- read_block(y, y_viewport)
        block_ans <- x %*% block
        write_block(sink, ans_grid[[b]], block_ans)
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
    sink <- RealizationSink(ans_dim, ans_dimnames, ans_type)
    on.exit(close(sink))
}

setMethod("%*%", c("DelayedMatrix", "DelayedMatrix"), .BLOCK_matrix_mult)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Costing calculations for parallelized multiplication
###
### by Aaron Lun
###
### We make some careful decisions about how to parallelize %*% operations.
### This is done to minimize the redundant data processing, especially when
### chunks of data can only be read in an atomic manner.
###

.define_atomic_chunks <- function(x)
# Defines atomic chunks across rows or columns of 'x'. These are subsets of
# contiguous rows or columns that must be read in their entirety. Wholly
# in-memory matrices have trivial boundaries, i.e., per row/column. Chunks
# for HDF5 matrices are defined by the physical layout.
{
    grid <- chunkGrid(x)
    if (is.null(grid)) {
        row_out <- list(
            bound=seq_len(nrow(x)),
            size=rep(1L, nrow(x))
        )
        row_out$cost <- row_out$size

        col_out <- list(
            bound=seq_len(ncol(x)),
            size=rep(1L, ncol(x))
        )
        col_out$cost <- col_out$size

    } else {
        N <- dim(grid)
        bygrid <- dims(grid)

        row_dims <- bygrid[seq_len(N[1]),1]
        row_out <- list(
            bound=cumsum(row_dims),
            size=row_dims
        )
        # file-backed penalty (relative to 1).
        # TODO: accommodate penalties for combined seeds that are only partially file-backed.
        row_out$cost <- row_out$size * 10L

        col_dims <- bygrid[seq_len(N[2]) * N[1],2]
        col_out <- list(
            bound=cumsum(col_dims),
            size=col_dims
        )
        col_out$cost <- col_out$size * 10L
    }

    list(row=row_out, col=col_out)
}

.grid_by_dimension <- function(x, nworkers)
# Splits a dimension of the matrix into at least 'nworkers' blocks.
# If the block size is too large, it is reduced to obtain the desired
# number of blocks in order for parallelization to be effective.
{
    old <- getAutoBlockSize()
    element_size <- switch(type(x), double=8, 4)

    ideal_size_by_row <- max(1, ceiling(nrow(x)/nworkers) * as.double(ncol(x)) * element_size)
    if (old > ideal_size_by_row) {
        suppressMessages(setAutoBlockSize(ideal_size_by_row)) # TODO: avoid setting the block size just to compute this?
        row_grid <- rowGrid(x)
        suppressMessages(setAutoBlockSize(old))
    } else {
        row_grid <- rowGrid(x)
    }

    ideal_size_by_col <- max(1, ceiling(ncol(x)/nworkers) * as.double(nrow(x)) * 8)
    if (old > ideal_size_by_col) {
        suppressMessages(setAutoBlockSize(ideal_size_by_col))
        col_grid <- colGrid(x)
        suppressMessages(setAutoBlockSize(old))
    } else {
        col_grid <- colGrid(x)
    }

    list(row=row_grid, col=col_grid)
}

.evaluate_split_costs <- function(split, dim_info)
# Given a split of a matrix dimension into putative blocks, 
# evaluate the costs with respect to chunk atomicity.
{
    boundaries <- dim_info$bound
    if (length(boundaries)==0L) { return(0) } # avoid NAs from zero-length dim.
    chunk_cost <- dim_info$cost

    right_chunk_index <- findInterval(split, boundaries, left.open=TRUE) + 1L
    block_starts <- c(1L, head(split, -1L) + 1L)
    left_chunk_index <- findInterval(block_starts, boundaries, left.open=TRUE) + 1L

    cum_cost <- c(0L, cumsum(chunk_cost))
    cum_cost[right_chunk_index+1L] - cum_cost[left_chunk_index]
}

.dispatch_mult <- function(x, y, nworkers, transposed.x=FALSE, transposed.y=FALSE)
# Decides how to split jobs across multiple workers for x %*% y
# or equivalent operations involving their transposes.
{
    atomic_x <- .define_atomic_chunks(x)
    nr_x <- as.double(nrow(x)) # use double to avoid overflow.
    nc_x <- as.double(ncol(x))

    atomic_y <- .define_atomic_chunks(y)
    nr_y <- as.double(nrow(y))
    nc_y <- as.double(ncol(y))

    x_grids <- .grid_by_dimension(x, nworkers)
    y_grids <- .grid_by_dimension(y, nworkers)

    # Mimicking the effect of supplying t(x) or t(y) (without actually computing them).
    if (transposed.x) {
        atomic_x <- list(row=atomic_x$col, col=atomic_x$row)
        x_grids <- list(row=t(x_grids$col), col=t(x_grids$row))
        tmp <- nr_x
        nr_x <- nc_x
        nc_x <- tmp
    }

    if (transposed.y) {
        atomic_y <- list(row=atomic_y$col, col=atomic_y$row)
        y_grids <- list(row=t(y_grids$col), col=t(y_grids$row))
        tmp <- nr_y
        nr_y <- nc_y
        nc_y <- tmp
    }

    # Determining costs for each job splitting scheme.
    # Assessment is based on the maximum cost for a given worker under each scheme,
    # with the aim being to minimize the maximum time to complete all jobs.
    x_by_row <- cumsum(dims(x_grids$row)[,1])
    x_by_row_cost <- .evaluate_split_costs(x_by_row, atomic_x$row)
    x_by_col <- cumsum(dims(x_grids$col)[,2])

    y_by_row <- cumsum(dims(y_grids$row)[,1])
    y_by_col <- cumsum(dims(y_grids$col)[,2])
    y_by_col_cost <- .evaluate_split_costs(y_by_col, atomic_y$col)

    cost_x_by_row <- max(x_by_row_cost * nc_x) + nr_y * sum(y_by_col_cost)
    cost_y_by_col <- sum(x_by_row_cost) * nc_x + max(nr_y * y_by_col_cost)

    if (getAutoMultParallelAgnostic()) { 
        # Ignoring splitting by the common dimension, as this is not
        # agnostic with respect to the number of workers.
        all_options <- c(x_by_row=cost_x_by_row, y_by_col=cost_y_by_col)
            
    } else {
        cost_common_x <- max(
            nr_x * .evaluate_split_costs(x_by_col, atomic_x$col) +
            .evaluate_split_costs(x_by_col, atomic_y$row) * nc_y
        )
        cost_common_y <- max(
            nr_x * .evaluate_split_costs(y_by_row, atomic_x$col) +
            .evaluate_split_costs(y_by_row, atomic_y$row) * nc_y
        )

        all_options <- c(x_by_row=cost_x_by_row, y_by_col=cost_y_by_col, common_x=cost_common_x, common_y=cost_common_y)
    } 

    choice <- names(all_options)[which.min(all_options)]

    if (choice=="x_by_row" || choice=="y_by_col") {
        # Both are returned though only one is used at any given time;
        # see .super_BLOCK_mult() for why we need both available.
        grid <- list(x=x_grids$row, y=y_grids$col)

    } else if (choice=="common_x") {
        grid <- list(x=x_grids$col,
            y=ArbitraryArrayGrid(list(
                cumsum(dims(x_grids$col)[,2]),
                as.integer(nc_y)
            ))
        )
    } else {
        grid <- list(
            x=ArbitraryArrayGrid(list(
                as.integer(nr_x),
                cumsum(dims(y_grids$row)[,1])
            )),
            y=y_grids$row)
    }

    # Undo the transposition to get viewports for the original matrix.
    if (transposed.x) {
        grid$x <- t(grid$x)
    }
    if (transposed.y) {
        grid$y <- t(grid$y)
    }
    choice <- switch(choice, x_by_row="x", y_by_col="y", common_x="common", common_y="common")
    list(choice=choice, grid=grid)
}

.dispatch_self <- function(x, nworkers, transposed=FALSE)
# Decides how to split jobs across multiple workers for crossprod(x).
{
    atomic <- .define_atomic_chunks(x)
    nr <- as.double(nrow(x)) # use double to avoid overflow.
    nc <- as.double(ncol(x))
    grids <- .grid_by_dimension(x, nworkers)

    if (transposed) {
        atomic <- list(row=atomic$col, col=atomic$row)
        grids <- list(row=t(grids$col), col=t(grids$row))
        tmp <- nr
        nr <- nc
        nc <- tmp
    }

    if (getAutoMultParallelAgnostic()) {
        choice <- "single"
        grid <- grids$col

    } else {
        # Splitting 'x' by column, and taking the crossproduct with all of 'x'.
        # This is the same the other way around, so we don't bother computing that.
        by_col <- cumsum(dims(grids$col)[,2])
        by_col_cost <- .evaluate_split_costs(by_col, atomic$col)
        cost_single <- max(nr * by_col_cost) + nr * sum(by_col_cost)
    
        # Splitting by the common dimension, i.e., along the rows of 'x'
        # and taking the crossproduct of each split block.
        by_row <- cumsum(dims(grids$row)[,1])
        by_row_cost <- .evaluate_split_costs(by_row, atomic$row)
        cost_both <- max(by_row_cost * nc)
    
        all_options <- c(both=cost_both, single=cost_single)
        choice <- names(all_options)[which.min(all_options)]
        grid <- switch(choice, both=grids$row, single=grids$col)
    }

    # Obtaining the grid for the original 'x'.
    if (transposed) {
        grid <- t(grid)
    }
    list(choice=choice, grid=grid)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Parallelized schemes for matrix multiplication.
###
### by Aaron Lun
###
### This splits one or both matrices into blocks according to the
### desired parallelization scheme, and distributes them to workers.
### This also requires care to respect the maximum block size.
###

.single_grid_iterate <- function(grid) {
    force(grid)
    b <- 0L

    # Do NOT put read_block in here, as this will waste time realizing from file in the single worker.
    # TODO: subset in-memory matrices to avoid serializing them in their entirety to the workers.
    function () {
        if (b==length(grid)) return(NULL)
        b <<- b + 1L
        grid[[b]]
    }
}

.left_mult <- function(viewport, x, y, MULT) {
    block <- read_block(x, viewport) # this, and all other calls, had better yield a non-DA, otherwise MULT will recurse endlessly.
    MULT(block, y)
}

.right_mult <- function(viewport, x, y, MULT) {
    block <- read_block(y, viewport)
    MULT(x, block)
}

.dual_grid_iterate <- function(grid.x, grid.y) {
    force(grid.x)
    force(grid.y)
    b <- 0L

    function () {
        if (b==length(grid.x)) return(NULL)
        b <<- b + 1L
        list(x=grid.x[[b]], y=grid.y[[b]])
    }
}

.both_mult <- function(viewports, x, y, MULT) {
    block_x <- read_block(x, viewports$x)
    block_y <- read_block(y, viewports$y)
    MULT(block_x, block_y)
}

.super_BLOCK_mult <- function(x, y, MULT, transposed.x=FALSE, transposed.y=FALSE, BPPARAM=getAutoBPPARAM())
# Controller function that split jobs for a multiplication function "MULT".
# This accommodates %*%, crossprod and tcrossprod for two arguments.
{
    scheme_out <- .dispatch_mult(x, y, bpnworkers(BPPARAM), transposed.x=transposed.x, transposed.y=transposed.y)
    chosen_scheme <- scheme_out$choice
    chosen_grids <- scheme_out$grid

    old <- getAutoBPPARAM()
    on.exit(setAutoBPPARAM(old))
    setAutoBPPARAM(SerialParam()) # avoid re-parallelizing in further calls to 'MULT'.

    if (chosen_scheme=="common") {
        ans <- bpiterate(.dual_grid_iterate(chosen_grids$x, chosen_grids$y), FUN=.both_mult,
            x=x, y=y, MULT=MULT, BPPARAM=BPPARAM, REDUCE=function(x, y) x + y)

    } else {
        # Switch to iteration over the other argument if the chosen one is
        # single-block and non-DA (at which point you might as well iterate
        # over the other argument anyway). This avoids infinite recursion
        # when 'x' or 'y' fail to get realized via read_block().
        if (chosen_scheme=="x" && length(chosen_grids$x)==1L && !is(x, "DelayedMatrix")) {
            chosen_scheme <- "y"
        } else if (chosen_scheme=="y" && length(chosen_grids$y)==1L && !is(y, "DelayedMatrix")) {
            chosen_scheme <- "x"
        }

        if (chosen_scheme=="x") {
            out <- bpiterate(.single_grid_iterate(chosen_grids$x), FUN=.left_mult,
                x=x, y=y, MULT=MULT, BPPARAM=BPPARAM)
            ans <- do.call(rbind, out)
        } else if (chosen_scheme=="y") {
           out <- bpiterate(.single_grid_iterate(chosen_grids$y), FUN=.right_mult,
                x=x, y=y, MULT=MULT, BPPARAM=BPPARAM)
           ans <- do.call(cbind, out)
        }
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

setMethod("crossprod", c("DelayedMatrix", "DelayedMatrix"), function(x, y) .super_BLOCK_mult(x, y, MULT=crossprod, transposed.x=TRUE))

# tcrossprod with vector 'y' doesn't work in base, and so it won't work here either.
setMethod("tcrossprod", c("DelayedMatrix", "ANY"), function(x, y) .super_BLOCK_mult(x, y, MULT=tcrossprod, transposed.y=TRUE))

setMethod("tcrossprod", c("ANY", "DelayedMatrix"), function(x, y) {
    if (is.null(dim(x))) x <- rbind(x)
    .super_BLOCK_mult(x, y, MULT=tcrossprod, transposed.y=TRUE)
})

setMethod("tcrossprod", c("DelayedMatrix", "DelayedMatrix"), function(x, y) .super_BLOCK_mult(x, y, MULT=tcrossprod, transposed.y=TRUE))

.solo_mult <- function(viewport, x, MULT) {
    block <- read_block(x, viewport)
    MULT(block)
}

.super_BLOCK_self <- function(x, MULT, transposed=FALSE, BPPARAM=getAutoBPPARAM())
# Controller function that split jobs for a multiplication function "MULT".
# This accommodates crossprod and tcrossprod for single arguments.
{
    scheme_out <- .dispatch_self(x, bpnworkers(BPPARAM), transposed=transposed)
    chosen_scheme <- scheme_out$choice
    chosen_grids <- scheme_out$grid

    old <- getAutoBPPARAM()
    on.exit(setAutoBPPARAM(old))
    setAutoBPPARAM(SerialParam()) # avoid re-parallelizing in further calls to 'MULT'.

    if (chosen_scheme=="single") {
        out <- bpiterate(.single_grid_iterate(chosen_grids), FUN=.left_mult,
            x=x, y=x, MULT=MULT, BPPARAM=SerialParam())
        ans <- do.call(rbind, out)

    } else if (chosen_scheme=="both") {
        ans <- bpiterate(.single_grid_iterate(chosen_grids), FUN=.solo_mult,
            x=x, MULT=MULT, BPPARAM=BPPARAM, REDUCE=function(x, y) x+y)

    }

    realize(ans)
}

setMethod("crossprod", c("DelayedMatrix", "missing"), function(x, y) .super_BLOCK_self(x, MULT=crossprod))

setMethod("tcrossprod", c("DelayedMatrix", "missing"), function(x, y) .super_BLOCK_self(x, MULT=tcrossprod, transposed=TRUE))

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
