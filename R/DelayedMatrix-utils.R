### =========================================================================
### Common operations on DelayedMatrix objects
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### rowSums(), colSums(), rowMeans(), colMeans()
###

.normarg_dims <- function(dims, method)
{
    if (!identical(dims, 1))
        stop("\"", method, "\" method for DelayedMatrix objects ",
             "does not support the 'dims' argument yet")
}

.DelayedMatrix_block_rowSums <- function(x, na.rm=FALSE, dims=1)
{
    .normarg_dims(dims, "rowSums")
    if (is(x, "DelayedArray") && x@is_transposed)
        return(.DelayedMatrix_block_colSums(t(x), na.rm=na.rm, dims=dims))

    REDUCE <- function(submatrix) rowSums(submatrix, na.rm=na.rm)
    COMBINE <- function(i, subarray, init, reduced) { init + reduced }
    init <- numeric(nrow(x))
    ans <- colblock_REDUCE_and_COMBINE(x, REDUCE, COMBINE, init)
    setNames(ans, rownames(x))
}

.DelayedMatrix_block_colSums <- function(x, na.rm=FALSE, dims=1)
{
    .normarg_dims(dims, "colSums")
    if (is(x, "DelayedArray") && x@is_transposed)
        return(.DelayedMatrix_block_rowSums(t(x), na.rm=na.rm, dims=dims))

    colsums_list <- colblock_APPLY(x, colSums, na.rm=na.rm,
                                   if_empty=numeric(0))
    unlist(colsums_list, recursive=FALSE)
}

setMethod("rowSums", "DelayedMatrix", .DelayedMatrix_block_rowSums)
setMethod("colSums", "DelayedMatrix", .DelayedMatrix_block_colSums)

.DelayedMatrix_block_rowMeans <- function(x, na.rm=FALSE, dims=1)
{
    .normarg_dims(dims, "rowMeans")
    if (is(x, "DelayedMatrix") && x@is_transposed)
        return(.DelayedMatrix_block_colMeans(t(x), na.rm=na.rm, dims=dims))

    REDUCE <- function(submatrix) {
        submatrix_sums <- rowSums(submatrix, na.rm=na.rm)
        submatrix_nvals <- ncol(submatrix)
        if (na.rm)
            submatrix_nvals <- submatrix_nvals - rowSums(is.na(submatrix))
        cbind(submatrix_sums, submatrix_nvals)
    }
    COMBINE <- function(i, subarray, init, reduced) { init + reduced }
    init <- cbind(
        numeric(nrow(x)),  # sums
        numeric(nrow(x))   # nvals
    )
    ans <- colblock_REDUCE_and_COMBINE(x, REDUCE, COMBINE, init)
    setNames(ans[ , 1L] / ans[ , 2L], rownames(x))
}

.DelayedMatrix_block_colMeans <- function(x, na.rm=FALSE, dims=1)
{
    .normarg_dims(dims, "colMeans")
    if (is(x, "DelayedArray") && x@is_transposed)
        return(.DelayedMatrix_block_rowMeans(t(x), na.rm=na.rm, dims=dims))

    colmeans_list <- colblock_APPLY(x, colMeans, na.rm=na.rm,
                                    if_empty=numeric(0))
    unlist(colmeans_list, recursive=FALSE)
}

setMethod("rowMeans", "DelayedMatrix", .DelayedMatrix_block_rowMeans)
setMethod("colMeans", "DelayedMatrix", .DelayedMatrix_block_colMeans)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Matrix multiplication
###
### We only support multiplication of an ordinary matrix (typically
### small) by a DelayedMatrix object (typically big). Multiplication of 2
### DelayedMatrix objects is not supported.
###

### Write a new HDF5 dataset to disk. Return an HDF5Matrix object that points
### to this new dataset.
.DelayedMatrix_block_mult_by_left_matrix <- function(x, y)
{
    stopifnot(is.matrix(x),
              is(y, "DelayedMatrix") || is.matrix(y),
              ncol(x) == nrow(y))

    require_HDF5Array()
    out_file <- HDF5Array::getHDF5DumpFile()
    out_name <- HDF5Array::getHDF5DumpName()

    ans_type <- typeof(match.fun(type(x))(1) * match.fun(type(y))(1))
    HDF5Array:::h5createDataset2(out_file, out_name, c(nrow(x), ncol(y)),
                                 storage.mode=ans_type)
    on.exit(HDF5Array::setHDF5DumpName())

    colblock_APPLY(y,
        function(submatrix) x %*% submatrix,
        out_file=out_file,
        out_name=out_name
    )

    ## TODO: Investigate the possiblity to store the dimnames in the HDF5 file
    ## so the HDF5Array() constructor can bring them back. Then we wouldn't
    ## need to explicitely set them on 'ans' like we do below.
    ans <- HDF5Array::HDF5Array(out_file, out_name, type=ans_type)
    ans_rownames <- rownames(x)
    ans_colnames <- colnames(y)
    if (!(is.null(ans_rownames) && is.null(ans_colnames)))
        dimnames(ans) <- list(ans_rownames, ans_colnames)
    ans
}

setMethod("%*%", c("DelayedMatrix", "matrix"),
    function(x, y) t(t(y) %*% t(x))
)

setMethod("%*%", c("matrix", "DelayedMatrix"),
    .DelayedMatrix_block_mult_by_left_matrix
)

setMethod("%*%", c("DelayedMatrix", "DelayedMatrix"),
    function(x, y)
        stop(wmsg("multiplication of 2 DelayedMatrix objects is not ",
                  "supported, only multiplication of an ordinary matrix by ",
                  "a DelayedMatrix object at the moment"))
)

