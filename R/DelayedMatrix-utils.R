### =========================================================================
### Common operations on DelayedMatrix objects
### -------------------------------------------------------------------------
###


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

    ideal_size_by_row <- ceiling(nrow(x)/nworkers) * as.double(ncol(x)) * element_size 
    if (old > ideal_size_by_row) {
        setAutoBlockSize(ideal_size_by_row) # TODO: avoid need to set the block size to just compute this.
        row_grid <- rowGrid(x)
        setAutoBlockSize(old)
    } else {
        row_grid <- rowGrid(x)
    }

    ideal_size_by_col <- ceiling(ncol(x)/nworkers) * as.double(nrow(x)) * 8 
    if (old > ideal_size_by_col) {
        setAutoBlockSize(ideal_size_by_col)
        col_grid <- colGrid(x)
        setAutoBlockSize(old)
    } else {
        col_grid <- colGrid(x)
    }

    list(row=row_grid, col=col_grid)
}

.evaluate_split_costs <- function(split, dim_info) 
# Given a split, evaluate the costs with respect to chunk atomicity. 
{
    boundaries <- dim_info$bound
    chunk_cost <- dim_info$cost

    last_chunk <- findInterval(split, boundaries, left.open=TRUE) + 1L
    prev_chunk <- findInterval(c(1L, head(split, -1L) + 1L), boundaries, left.open=TRUE)

    cum_cost <- c(0L, cumsum(chunk_cost))
    cum_cost[last_chunk+1L] - cum_cost[prev_chunk+1L]
}

.dispatch_mult <- function(x, y, nworkers, transposed.x=FALSE, transposed.y=FALSE) 
# Decides how to split jobs across multiple workers for x %*% y.
{
    atomic_x <- .define_atomic_chunks(x)
    nr_x <- as.double(nrow(x)) # use double to avoid overflow.
    nc_x <- as.double(ncol(x))

    atomic_y <- .define_atomic_chunks(y)
    nr_y <- as.double(nrow(y))
    nc_y <- as.double(ncol(y))
    
    x_grids <- .grid_by_dimension(x, nworkers)
    y_grids <- .grid_by_dimension(y, nworkers)

    # Mimicking the effect of supplying t(x) or t(y) (without actually computing that).
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
    x_by_row <- cumsum(dims(x_grids$row)[,1])
    x_by_row_cost <- .evaluate_split_costs(x_by_row, atomic_x$row)
    x_by_col <- cumsum(dims(x_grids$col)[,2])

    y_by_row <- cumsum(dims(y_grids$row)[,1])
    y_by_col <- cumsum(dims(y_grids$col)[,2])
    y_by_col_cost <- .evaluate_split_costs(y_by_col, atomic_y$col)

    cost_x_by_row <- max(x_by_row_cost * nc_x) + nr_y * sum(y_by_col_cost)
    cost_y_by_col <- sum(x_by_row_cost) * nc_x + max(nr_y * y_by_col_cost)
    cost_common_x <- max(
        nr_x * .evaluate_split_costs(x_by_col, atomic_x$col) + 
        .evaluate_split_costs(x_by_col, atomic_y$row) * nc_y
    )
    cost_common_y <- max(
        nr_x * .evaluate_split_costs(y_by_row, atomic_x$col) + 
        .evaluate_split_costs(y_by_row, atomic_y$row) * nc_y
    )

    # Assessment is based on the maximum cost for a given worker under each scheme,
    # with the aim being to minimize the maximum time to complete all jobs.
    all_options <- c(x_by_row=cost_x_by_row, y_by_col=cost_y_by_col, common_x=cost_common_x, common_y=cost_common_y)
    choice <- names(all_options)[which.min(all_options)]

    if (choice=="x_by_row") {
        grid <- list(x=x_grids$row)
    } else if (choice=="y_by_row") {
        grid <- list(y=y_grids$col)
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
    if (transposed.x && !is.null(grid$x)) {
        grid$x <- t(grid$x)
    }
    if (transposed.y && !is.null(grid$y)) {
        grid$y <- t(grid$y)
    }
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

    by_row <- cumsum(dims(grids$row)[,1])
    by_row_cost <- .evaluate_split_costs(by_row, atomic_x$row)
    by_col <- cumsum(dims(grids$col)[,2])
    by_col_cost <- .evaluate_split_costs(by_col, atomic_x$col)

    # Splitting 'x' by column, and taking the crossproduct with all of 'x'.
    # This is the same the other way around, so we don't bother computing that.
    cost_single <- max(nr * by_col_cost) + nr * sum(by_col_cost)

    # Splitting by the common dimension, i.e., along the rows of 'x'
    # and taking the crossproduct of each split block.
    cost_both <- max(by_row_cost * nc) * 2

    all_options <- c(both=cost_both, single=cost_single)
    choice <- names(all_options)[which.min(all_options)]
    grid <- switch(choice, both=grids$row, single=grids$col)
    if (transposed) {
        grid <- t(grid)
    }
    list(choice=choice, grid=grid)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Parallelized schemes for matrix multiplication.
###
### This splits one or both matrices into blocks according to the
### desired parallelization scheme, and distributes them to workers.
### This also requires care to respect the maximum block size.
###

.super_BLOCK_mult <- function(x, y, MULT, transposed.x=FALSE, transposed.y=FALSE, BPPARAM=bpparam()) {
    if (!is(x, "DelayedMatrix") && !is(y, "DelayedMatrix")) {
        return(MULT(x, y))
    }

    scheme_out <- .dispatch_mult(x, y, bpnworkers(BPPARAM), transposed.x=transposed.x, transposed.y=transposed.y)
    chosen_scheme <- scheme_out$choice
    chosen_grids <- scheme_out$grid

    old <- bpparam()
    on.exit(register(old))
    register(SerialParam()) # avoid re-parallelizing in further calls to 'MULT'.

    if (chosen_scheme=="x_by_row") {
        # TODO: figure out best serialization scheme.
        out <- bplapply(chosen_grids[[1]], FUN=function(viewport, x, y) {
            block <- read_block(x, viewport) # this, and all other calls, had better yield a non-DA, otherwise MULT will recurse endlessly.
            MULT(block, y)
        }, x=x, y=y, BPPARAM=BPPARAM)
        ans <- do.call(rbind, out)

    } else if (chosen_scheme=="y_by_col") {
        out <- bplapply(chosen_grids[[1]], FUN=function(viewport, x, y) {
            block <- read_block(y, viewport)
            MULT(x, block)
        }, x=x, y=y, BPPARAM=BPPARAM)
        ans <- do.call(cbind, out)

    } else if (chosen_scheme=="common_x") {
        # TODO: switch to bpiterate to reduce products once they are created.
        out <- bpmapply(viewport_x=chosen_grids[[1]], viewport_y=chosen_grids[[2]], FUN=function(viewport_x, viewport_y, x, y) {
            block_x <- read_block(x, viewport_x)
            block_y <- read_block(y, viewport_y)
            MULT(block_x, block_y)
        }, MoreArgs=list(x=x, y=y), BPPARAM=BPPARAM, SIMPLIFY=FALSE)
        ans <- Reduce("+", out)

    } else  {
        out <- bpmapply(viewport_x=chosen_grids[[1]], viewport_y=chosen_grids[[2]], FUN=function(viewport_x, viewport_y, x, y) {
            block_x <- read_block(x, viewport_x)
            block_y <- read_block(y, viewport_y)
            MULT(block_x, block_y)
        }, MoreArgs=list(x=x, y=y), BPPARAM=BPPARAM, SIMPLIFY=FALSE)
        ans <- Reduce("+", out)

    }

    ans
}

setMethod("%*%", c("DelayedMatrix", "ANY"), function(x, y) .super_BLOCK_mult(x, y, MULT=`%*%`))

setMethod("%*%", c("ANY", "DelayedMatrix"), function(x, y) .super_BLOCK_mult(x, y, MULT=`%*%`))

setMethod("%*%", c("DelayedMatrix", "DelayedMatrix"), function(x, y) .super_BLOCK_mult(x, y, MULT=`%*%`))

setMethod("crossprod", c("DelayedMatrix", "ANY"), function(x, y) .super_BLOCK_mult(x, y, MULT=crossprod, transposed.x=TRUE))

setMethod("crossprod", c("ANY", "DelayedMatrix"), function(x, y) .super_BLOCK_mult(x, y, MULT=crossprod, transposed.x=TRUE))

setMethod("crossprod", c("DelayedMatrix", "DelayedMatrix"), function(x, y) .super_BLOCK_mult(x, y, MULT=crossprod, transposed.x=TRUE))

setMethod("tcrossprod", c("DelayedMatrix", "ANY"), function(x, y) .super_BLOCK_mult(x, y, MULT=tcrossprod, transposed.y=TRUE))

setMethod("tcrossprod", c("ANY", "DelayedMatrix"), function(x, y) .super_BLOCK_mult(x, y, MULT=tcrossprod, transposed.y=TRUE))

setMethod("tcrossprod", c("DelayedMatrix", "DelayedMatrix"), function(x, y) .super_BLOCK_mult(x, y, MULT=tcrossprod, transposed.y=TRUE))
