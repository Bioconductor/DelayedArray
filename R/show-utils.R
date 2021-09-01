### =========================================================================
### Compact display of an array-like object
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Some wrappers around extract_array_by_Nindex()
###
### All of them return a list of 2 ordinary arrays. The dimnames need to
### propagate so we use extract_array_by_Nindex() instead of extract_array().
###

.extract_two_1Darrays <- function(x, i1, i2)
{
    a <- extract_array_by_Nindex(x, list(c(i1, i2)))
    ii1 <- seq_along(i1)
    ii2 <- seq_along(i2) + length(i1)
    x1 <- a[ii1, drop=FALSE]
    x2 <- a[ii2, drop=FALSE]
    list(x1, x2)
}

.extract_two_matrices_by_row <- function(x, i1, i2)
{
    a <- extract_array_by_Nindex(x, list(c(i1, i2), NULL))
    ii1 <- seq_along(i1)
    ii2 <- seq_along(i2) + length(i1)
    x1 <- a[ii1, , drop=FALSE]
    x2 <- a[ii2, , drop=FALSE]
    list(x1, x2)
}

.extract_two_matrices_by_col <- function(x, j1, j2)
{
    a <- extract_array_by_Nindex(x, list(NULL, c(j1, j2)))
    jj1 <- seq_along(j1)
    jj2 <- seq_along(j2) + length(j1)
    x1 <- a[ , jj1, drop=FALSE]
    x2 <- a[ , jj2, drop=FALSE]
    list(x1, x2)
}

.extract_four_matrices <- function(x, i1, i2, j1, j2)
{
    a <- extract_array_by_Nindex(x, list(c(i1, i2), c(j1, j2)))
    ii1 <- seq_along(i1)
    ii2 <- seq_along(i2) + length(i1)
    jj1 <- seq_along(j1)
    jj2 <- seq_along(j2) + length(j1)
    x11 <- a[ii1, jj1, drop=FALSE]
    x21 <- a[ii2, jj1, drop=FALSE]
    x12 <- a[ii1, jj2, drop=FALSE]
    x22 <- a[ii2, jj2, drop=FALSE]
    list(x11, x21, x12, x22)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Some wrappers around format()
###

### 'x' must be an ordinary vector or matrix of atomic or recursive type.
### 'max.width' takes effect only if 'x' is character or list (i.e. if the
### underlying type inherits from character or list when 'x' is a matrix).
.format2 <- function(x, justify, quote=TRUE, max.width=22L)
{
    if (is.character(x) && length(x) != 0L && quote) {
        na_idx <- which(is.na(x))
        x <- paste0("\"", x, "\"")
        x[na_idx] <- "NA"
    }
    ans <- format(x, justify=justify)
    if ((is.character(x) || is.list(x)) && length(x) != 0L) {
        truncate_me <- nchar(ans) > max.width
        stop <- ifelse(truncate_me, max.width - 2L, nchar(ans))
        ans <- substr(ans, 1L, stop)
        ans <- paste0(ans, ifelse(truncate_me, "..", ""))
    }
    ans
}

.format_as_character_vector <- function(x, justify, quote=TRUE)
{
    x_names <- names(x)
    x <- as.vector(x)
    ans <- .format2(x, justify, quote)
    names(ans) <- x_names
    ans
}

.format_as_character_matrix <- function(x, justify, quote=TRUE)
{
    x_dim <- dim(x)
    x_dimnames <- dimnames(x)
    x <- as.matrix(x)
    ans <- .format2(x, justify, quote)
    ans <- set_dim(ans, x_dim)
    set_dimnames(ans, x_dimnames)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Some low-level utilities to handle array names/dimnames
###

.wrap_subscripts <- function(i, x_names=NULL, left="[", right="]")
{
    if (length(i) == 0L)
        return(character(0))
    if (is.null(x_names))
        return(paste0(left, i, right, sep=""))
    x_names[i]
}

.make_1Darray_names <- function(i1, i2, x_names, justify)
{
    s1 <- .wrap_subscripts(i1, x_names)
    s2 <- .wrap_subscripts(i2, x_names)
    format(c(s1, ".", s2), justify=justify)
}

.make_rownames <- function(i1, i2, x_rownames, justify)
{
    s1 <- .wrap_subscripts(i1, x_rownames, right=",]")
    s2 <- .wrap_subscripts(i2, x_rownames, right=",]")
    max_width <- max(nchar(s1, type="width"), nchar(s2, type="width"))
    if (max_width <= 1L) {
        ellipsis <- "."
    } else if (max_width == 2L) {
        ellipsis <- ".."
    } else {
        ellipsis <- "..."
    }
    format(c(s1, ellipsis, s2), justify=justify)
}

.make_colnames <- function(j1, j2, x_colnames, justify)
{
    s1 <- .wrap_subscripts(j1, x_colnames, left="[,")
    s2 <- .wrap_subscripts(j2, x_colnames, left="[,")
    ans <- format(c(s1, s2), justify=justify)
    c(head(ans, n=length(s1)), "...", tail(ans, n=length(s2)))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The extract_and_stitch_* family of helpers
###

.extract_and_stitch_two_1Darrays <- function(x, i1, i2,
                                             justify, quote=TRUE)
{
    x12 <- .extract_two_1Darrays(x, i1, i2)
    ans1 <- .format_as_character_vector(x12[[1L]], justify, quote=quote)
    ans2 <- .format_as_character_vector(x12[[2L]], justify, quote=quote)

    ans <- c(ans1, ".", ans2)

    names(ans) <- .make_1Darray_names(i1, i2, names(x), justify)
    ans
}

.extract_and_stitch_two_matrices_by_row <- function(x, i1, i2,
                                                    justify, quote=TRUE)
{
    x12 <- .extract_two_matrices_by_row(x, i1, i2)
    ans1 <- .format_as_character_matrix(x12[[1L]], justify, quote=quote)
    ans2 <- .format_as_character_matrix(x12[[2L]], justify, quote=quote)

    hdots <- rep.int(".", ncol(ans1))
    ans <- rbind(ans1, matrix(hdots, nrow=1L), ans2)

    rownames(ans) <- .make_rownames(i1, i2, rownames(x), justify)
    ans
}

.extract_and_stitch_two_matrices_by_col <- function(x, j1, j2,
                                                    justify, quote=TRUE)
{
    x12 <- .extract_two_matrices_by_col(x, j1, j2)
    ans1 <- .format_as_character_matrix(x12[[1L]], justify, quote=quote)
    ans2 <- .format_as_character_matrix(x12[[2L]], justify, quote=quote)

    vdots <- rep.int(".", nrow(ans1))
    ans <- cbind(ans1, matrix(vdots, ncol=1L), ans2)

    colnames(ans) <- .make_colnames(j1, j2, colnames(x), justify)
    ans
}

.extract_and_stitch_four_matrices <- function(x, i1, i2, j1, j2,
                                              justify, quote=TRUE)
{
    x1234 <- .extract_four_matrices(x, i1, i2, j1, j2)
    ans11 <- .format_as_character_matrix(x1234[[1L]], justify, quote=quote)
    ans21 <- .format_as_character_matrix(x1234[[2L]], justify, quote=quote)
    ans12 <- .format_as_character_matrix(x1234[[3L]], justify, quote=quote)
    ans22 <- .format_as_character_matrix(x1234[[4L]], justify, quote=quote)

    hdots <- rep.int(".", ncol(ans11))
    ans1 <- rbind(ans11, matrix(hdots, nrow=1L), ans21)

    hdots <- rep.int(".", ncol(ans12))
    ans2 <- rbind(ans12, matrix(hdots, nrow=1L), ans22)

    vdots <- rep.int(".", nrow(ans1))
    ans <- cbind(ans1, matrix(vdots, ncol=1L), ans2)

    rownames(ans) <- .make_rownames(i1, i2, rownames(x), justify)
    colnames(ans) <- .make_colnames(j1, j2, colnames(x), justify)
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### 1D array
###

.prepare_1Darray_sample <- function(x, n1, n2, justify, quote=TRUE)
{
    x_len <- length(x)
    x_names <- names(x)
    if (x_len <= n1 + n2 + 1L) {
        ans <- .format_as_character_vector(x, justify, quote=quote)
        i1 <- seq_len(x_len)
        i2 <- integer(0)
        names(ans) <- .make_1Darray_names(i1, i2, x_names, justify)[i1]
    } else {
        i1 <- seq_len(n1)
        i2 <- seq(to=x_len, by=1L, length.out=n2)
        ans <- .extract_and_stitch_two_1Darrays(x, i1, i2, justify, quote=quote)
    }
    ans
}

.print_1Darray_data <- function(x, n1, n2, quote=TRUE)
{
    stopifnot(length(dim(x)) == 1L)
    right <- type(x) != "character"
    justify <- if (right) "right" else "left"
    out <- .prepare_1Darray_sample(x, n1, n2, justify, quote=quote)
    print(out, quote=FALSE, right=right, max=length(out))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### 2D array
###

.rsplit_2Darray_data <- function(x, m1, m2, justify, quote=TRUE)
{
    x_nrow <- nrow(x)
    i1 <- seq_len(m1)
    i2 <- seq(to=x_nrow, by=1L, length.out=m2)
    .extract_and_stitch_two_matrices_by_row(x, i1, i2, justify, quote=quote)
}

.csplit_2Darray_data <- function(x, n1, n2, justify, quote=TRUE)
{
    x_ncol <- ncol(x)
    j1 <- seq_len(n1)
    j2 <- seq(to=x_ncol, by=1L, length.out=n2)
    .extract_and_stitch_two_matrices_by_col(x, j1, j2, justify, quote=quote)
}

.split_2Darray_data <- function(x, m1, m2, n1, n2, justify, quote=TRUE)
{
    x_nrow <- nrow(x)
    i1 <- seq_len(m1)
    i2 <- seq(to=x_nrow, by=1L, length.out=m2)
    x_ncol <- ncol(x)
    j1 <- seq_len(n1)
    j2 <- seq(to=x_ncol, by=1L, length.out=n2)
    .extract_and_stitch_four_matrices(x, i1, i2, j1, j2, justify, quote=quote)
}

.prepare_2Darray_sample <- function(x, m1, m2, n1, n2, justify, quote=TRUE)
{
    ## An attempt at reducing the nb of columns to display when 'x' has
    ## dimnames so the object fits in getOption("width"). Won't necessarily
    ## pick up the optimal nb of columns so should be revisited at some point.
    x_rownames <- rownames(x)
    if (is.null(x_rownames)) {
        rownames_width <- 6L
    } else {
        rownames_width <- nchar(x_rownames[[1L]])
    }
    half_width <- (getOption("width") - rownames_width) / 2
    x_colnames <- colnames(x)
    if (!is.null(x_colnames)) {
        colnames1 <- head(x_colnames, n=n1)
        colnames2 <- tail(x_colnames, n=n2)
        n1 <- pmax(sum(cumsum(nchar(colnames1) + 1L) < half_width), 1L)
        n2 <- pmax(sum(cumsum(nchar(colnames2) + 1L) < half_width), 1L)
    }

    x_nrow <- nrow(x)
    x_ncol <- ncol(x)
    if (x_nrow <= m1 + m2 + 1L) {
        if (x_ncol <= n1 + n2 + 1L) {
            ans <- .format_as_character_matrix(x, justify, quote=quote)
            ## Only needed because of this bug in base::print.default:
            ##   https://stat.ethz.ch/pipermail/r-devel/2016-March/072479.html
            ## TODO: Remove when the bug is fixed.
            if (is.null(colnames(ans))) {
                j1 <- seq_len(ncol(ans))
                j2 <- integer(0)
                colnames(ans) <- .make_colnames(j1, j2, NULL, justify)[j1]
            }
        } else {
            ans <- .csplit_2Darray_data(x, n1, n2, justify, quote=quote)
        }
    } else {
        if (x_ncol <= n1 + n2 + 1L) {
            ans <- .rsplit_2Darray_data(x, m1, m2, justify, quote=quote)
            ## Only needed because of this bug in base::print.default:
            ##   https://stat.ethz.ch/pipermail/r-devel/2016-March/072479.html
            ## TODO: Remove when the bug is fixed.
            if (is.null(colnames(ans))) {
                j1 <- seq_len(ncol(ans))
                j2 <- integer(0)
                colnames(ans) <- .make_colnames(j1, j2, NULL, justify)[j1]
            }
        } else {
            ans <- .split_2Darray_data(x, m1, m2, n1, n2, justify, quote=quote)
        }
    }
    ans
}

.print_2Darray_data <- function(x, m1, m2, n1, n2, quote=TRUE)
{
    stopifnot(length(dim(x)) == 2L)
    right <- type(x) != "character"
    justify <- if (right) "right" else "left"
    out <- .prepare_2Darray_sample(x, m1, m2, n1, n2, justify, quote=quote)
    print(out, quote=FALSE, right=right, max=length(out))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Array of arbitrary dimensions
###

.print_2D_slices <- function(x, m1, m2, n1, n2, grid, idx, quote=TRUE)
{
    x_dimnames <- dimnames(x)
    for (i in idx) {
        viewport <- grid[[i]]
        s <- make_string_from_ArrayViewport(viewport, dimnames=x_dimnames,
                                            as.2Dslice=TRUE)
        cat(s, "\n", sep="")
        Nindex <- makeNindexFromArrayViewport(viewport)
        slice <- extract_array_by_Nindex(x, Nindex)
        slice <- set_dim(slice, dim(slice)[1:2])
        .print_2Darray_data(slice, m1, m2, n1, n2, quote=quote)
        cat("\n")
    }
}

.print_nDarray_data <- function(x, n1, n2, quote=TRUE)
{
    x_dim <- dim(x)
    x_nrow <- x_dim[[1L]]
    x_ncol <- x_dim[[2L]]
    if (x_ncol <= 5L) {
        if (x_nrow <= 3L) {
            m1 <- m2 <- 3L  # print all rows of each slice
            z1 <- z2 <- 3L  # print first 3 and last 3 slices
        } else {
            m1 <- m2 <- 2L  # print first 2 and last 2 rows of each slice
            z1 <- z2 <- 1L  # print only first and last slices
        }
    } else {
        if (x_nrow <= 3L) {
            m1 <- m2 <- 2L  # print first 2 and last 2 rows of each slice
            z1 <- z2 <- 2L  # print first 2 and last 2 slices
        } else {
            m1 <- m2 <- 2L  # print first 2 and last 2 rows of each slice
            z1 <- z2 <- 1L  # print only first and last slices
        }
    }
    spacings <- x_dim
    spacings[-(1:2)] <- 1L
    grid <- RegularArrayGrid(x_dim, spacings)
    nblock <- length(grid)  # should be equal to prod(x_dim[-(1:2)])
    if (nblock <= z1 + z2 + 1L) {
        idx <- seq_len(nblock)
        .print_2D_slices(x, m1, m2, n1, n2, grid, idx, quote=quote)
    } else {
        idx1 <- seq_len(z1)
        idx2 <- seq(to=nblock, by=1L, length.out=z2)
        .print_2D_slices(x, m1, m2, n1, n2, grid, idx1, quote=quote)
        cat("...\n\n")
        .print_2D_slices(x, m1, m2, n1, n2, grid, idx2, quote=quote)
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### show_compact_array()
###

.print_array_data <- function(x, n1, n2, quote=TRUE)
{
    x_dim <- dim(x)
    if (length(x_dim) == 1L)
        return(.print_1Darray_data(x, n1, n2, quote=quote))
    if (length(x_dim) == 2L) {
        nhead <- get_showHeadLines()
        ntail <- get_showTailLines()
        return(.print_2Darray_data(x, nhead, ntail, n1, n2, quote=quote))
    }
    .print_nDarray_data(x, n1, n2, quote=quote)
}

### NOT exported but used in the HDF5Array package!
array_as_one_line_summary <- function(x)
{
    x_dim <- dim(x)
    sprintf("<%s>%s %s of class %s and type \"%s\"",
            paste0(x_dim, collapse=" x "),
            if (is_sparse(x)) " sparse" else "",
            if (length(x_dim) == 2L) "matrix" else "array",
            class(x),
            type(x))
}

### Work on any array-like object that complies with the "seed contract" i.e.
### that supports dim(), dimnames(), and extract_array().
### NOT exported.
show_compact_array <- function(object)
{
    cat(array_as_one_line_summary(object))
    if (any(dim(object) == 0L)) {
        cat("\n")
        return()
    }
    cat(":\n")
    if (type(object) == "integer") {
        n1 <- n2 <- 4L
    } else {
        n1 <- 3L
        n2 <- 2L
    }
    .print_array_data(object, n1, n2)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### showAsCell()
###

### A re-implementation of S4Vectors:::.showAsCell_array() that avoids doing
### the "reshaping" of the object. Should work on any array-like object that
### supports extract_array().
.showAsCell_array_like <- function(object)
{
    dim1 <- dim(object)[-1L]
    p1 <- prod(dim1)
    if (p1 == 0L)
        return(rep.int("", nrow(object)))
    first_cols <- lapply(seq_len(min(p1, 3L)),
        function(i1) {
            subscripts <- arrayInd(i1, dim1)
            index <- c(list(NULL), as.list(subscripts))
            as.character(extract_array(object, index))
        }
    )
    ans <- do.call(paste, c(first_cols, list(sep=":")))
    if (p1 > 3L)
        ans <- paste0(ans, ":...")
    ans
}

setMethod("showAsCell", "Array", .showAsCell_array_like)

