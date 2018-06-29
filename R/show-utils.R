### =========================================================================
### Compact display of an array-like object
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


### 'x' must be an ordinary vector or matrix of atomic or recursive type.
### 'max.width' takes effect only if 'x' is character or list (i.e. if the
### underlying type inherits from character or list when 'x' is a matrix).
.format2 <- function(x, justify, quote=TRUE, max.width=22L)
{
    if (is.character(x) && length(x) != 0L && quote)
        x <- paste0("\"", x, "\"")
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
### 1D array
###

.split_1Darray_names <- function(x_names, idx1, idx2, justify)
{
    make_elt_indices <- function(i) {
        if (length(i) == 0L)
            return(character(0))
        paste0("[", i, "]", sep="")
    }
    if (is.null(x_names)) {
        s1 <- make_elt_indices(idx1)
        s2 <- make_elt_indices(idx2)
    } else {
        s1 <- x_names[idx1]
        s2 <- x_names[idx2]
    }
    format(c(s1, ".", s2), justify=justify)
}

.prepare_1Darray_sample <- function(x, n1, n2, justify, quote=TRUE)
{
    x_len <- length(x)
    x_names <- names(x)
    if (x_len <= n1 + n2 + 1L) {
        ans <- .format_as_character_vector(x, justify, quote=quote)
        idx1 <- seq_len(x_len)
        idx2 <- integer(0)
        names(ans) <- .split_1Darray_names(x_names, idx1, idx2, justify)[idx1]
    } else {
        idx1 <- seq_len(n1)
        idx2 <- seq(to=x_len, by=1L, length.out=n2)
        ## We use extract_array_by_Nindex() instead of extract_array() to
        ## propagate the dimnames.
        x1 <- extract_array_by_Nindex(x, list(idx1))
        x2 <- extract_array_by_Nindex(x, list(idx2))
        ans1 <- .format_as_character_vector(x1, justify, quote=quote)
        ans2 <- .format_as_character_vector(x2, justify, quote=quote)
        ans <- c(ans1, ".", ans2)
        names(ans) <- .split_1Darray_names(x_names, idx1, idx2, justify)
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

.split_rownames <- function(x_rownames, idx1, idx2, justify)
{
    make_row_indices <- function(i) {
        if (length(i) == 0L)
            return(character(0))
        paste0("[", i, ",]", sep="")
    }
    if (is.null(x_rownames)) {
       s1 <- make_row_indices(idx1)
       s2 <- make_row_indices(idx2)
    } else {
       s1 <- x_rownames[idx1]
       s2 <- x_rownames[idx2]
    }
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

.split_colnames <- function(x_colnames, idx1, idx2, justify)
{
    make_col_indices <- function(j) {
        if (length(j) == 0L)
            return(character(0))
        paste0("[,", j, "]", sep="")
    }
    if (is.null(x_colnames)) {
        s1 <- make_col_indices(idx1)
        s2 <- make_col_indices(idx2)
    } else {
        s1 <- x_colnames[idx1]
        s2 <- x_colnames[idx2]
    }
    ans <- format(c(s1, s2), justify=justify)
    c(head(ans, n=length(s1)), "...", tail(ans, n=length(s2)))
}

.rsplit_2Darray_data <- function(x, m1, m2, justify, quote=TRUE)
{
    x_nrow <- nrow(x)
    x_rownames <- rownames(x)
    idx1 <- seq_len(m1)
    idx2 <- seq(to=x_nrow, by=1L, length.out=m2)

    ## We use extract_array_by_Nindex() instead of extract_array() to
    ## propagate the dimnames.
    x1 <- extract_array_by_Nindex(x, list(idx1, NULL))
    x2 <- extract_array_by_Nindex(x, list(idx2, NULL))
    ans1 <- .format_as_character_matrix(x1, justify, quote=quote)
    ans2 <- .format_as_character_matrix(x2, justify, quote=quote)
    dots <- rep.int(".", ncol(ans1))
    ans <- rbind(ans1, matrix(dots, nrow=1L), ans2)

    rownames(ans) <- .split_rownames(x_rownames, idx1, idx2, justify)
    ans
}

.csplit_2Darray_data <- function(x, n1, n2, justify, quote=TRUE)
{
    x_ncol <- ncol(x)
    x_colnames <- colnames(x)
    idx1 <- seq_len(n1)
    idx2 <- seq(to=x_ncol, by=1L, length.out=n2)

    ## We use extract_array_by_Nindex() instead of extract_array() to
    ## propagate the dimnames.
    x1 <- extract_array_by_Nindex(x, list(NULL, idx1))
    x2 <- extract_array_by_Nindex(x, list(NULL, idx2))
    ans1 <- .format_as_character_matrix(x1, justify, quote=quote)
    ans2 <- .format_as_character_matrix(x2, justify, quote=quote)
    dots <- rep.int(".", nrow(ans1))
    ans <- cbind(ans1, matrix(dots, ncol=1L), ans2)

    colnames(ans) <- .split_colnames(x_colnames, idx1, idx2, justify)
    ans
}

.split_2Darray_data <- function(x, m1, m2, n1, n2, justify, quote=TRUE)
{
    x_ncol <- ncol(x)
    x_colnames <- colnames(x)
    idx1 <- seq_len(n1)
    idx2 <- seq(to=x_ncol, by=1L, length.out=n2)

    ## We use extract_array_by_Nindex() instead of extract_array() to
    ## propagate the dimnames.
    x1 <- extract_array_by_Nindex(x, list(NULL, idx1))
    x2 <- extract_array_by_Nindex(x, list(NULL, idx2))
    ans1 <- .rsplit_2Darray_data(x1, m1, m2, justify, quote=quote)
    ans2 <- .rsplit_2Darray_data(x2, m1, m2, justify, quote=quote)
    dots <- rep.int(".", nrow(ans1))
    ans <- cbind(ans1, matrix(dots, ncol=1L), ans2)

    colnames(ans) <- .split_colnames(x_colnames, idx1, idx2, justify)
    ans
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
                idx1 <- seq_len(ncol(ans))
                idx2 <- integer(0)
                colnames(ans) <- .split_colnames(NULL, idx1, idx2,
                                                 justify)[idx1]
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
                idx1 <- seq_len(ncol(ans))
                idx2 <- integer(0)
                colnames(ans) <- .split_colnames(NULL, idx1, idx2,
                                                 justify)[idx1]
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
        index <- makeNindexFromArrayViewport(viewport, expand.RangeNSBS=TRUE)
        ## We use extract_array_by_Nindex() instead of extract_array() to
        ## propagate the dimnames.
        slice <- extract_array_by_Nindex(x, index)
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
    grid <- makeRegularArrayGridOfCappedLengthViewports(x_dim,
                                                        prod(x_dim[1:2]),
                                                        "linear")
    nblock <- length(grid)
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

.array_one_line_summary <- function(x)
{
    dim_in1string <- paste0(dim(x), collapse=" x ")
    sprintf("<%s> %s object of type \"%s\"", dim_in1string, class(x), type(x))
}

### Work on any array-like object that complies with the "seed contract" i.e.
### that supports dim(), dimnames(), and extract_array().
show_compact_array <- function(object)
{
    cat(.array_one_line_summary(object))
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

