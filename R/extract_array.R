### =========================================================================
### extract_array()
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level helpers
###

### Return the slice as a list.
.extract_data_frame_slice <- function(x, index)
{
    slice <- subset_by_Nindex(x, index)
    ## Turn into a list and replace factors with character vectors.
    lapply(slice, as.vector)
}
.extract_DataFrame_slice <- function(x, index)
{
    slice <- subset_by_Nindex(x, index)
    slice <- as.data.frame(slice)
    ## Turn into a list and replace factors with character vectors.
    lapply(slice, as.vector)
}

### Return a list with one list element per column in data frame 'x'.
### All the list elements have length 0.
.extract_data_frame_slice0 <- function(x)
{
    slice0 <- x[0L, , drop=FALSE]
    ## Turn into a list and replace factors with character vectors.
    lapply(slice0, as.vector)
}
.extract_DataFrame_slice0 <- function(x)
{
    slice0 <- x[0L, , drop=FALSE]
    slice0 <- as.data.frame(slice0)
    if (ncol(slice0) != ncol(x))
        stop(wmsg("DataFrame object 'x' can be used as the seed of ",
                  "a DelayedArray object only if as.data.frame(x) ",
                  "preserves the number of columns"))
    ## Turn into a list and replace factors with character vectors.
    lapply(slice0, as.vector)
}

### Equivalent to 'typeof(as.matrix(x))' but with an almost-zero
### memory footprint (it avoids the cost of turning 'x' into a matrix).
.get_data_frame_type <- function(x)
{
    if (ncol(x) == 0L)
        return("logical")
    slice0 <- .extract_data_frame_slice0(x)
    typeof(unlist(slice0, use.names=FALSE))
}

### Equivalent to 'typeof(as.matrix(as.data.frame(x)))' but with an
### almost-zero memory footprint (it avoids the cost of turning 'x' first
### into a data frame then into a matrix).
.get_DataFrame_type <- function(x)
{
    if (ncol(x) == 0L)
        return("logical")
    slice0 <- .extract_DataFrame_slice0(x)
    typeof(unlist(slice0, use.names=FALSE))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### extract_array() generic and methods
###

### 'index' is expected to be an unnamed list of subscripts as positive
### integer vectors, one vector per dimension in 'x'. *Missing* list elements
### are allowed and represented by NULLs.
### The "extract_array" methods don't need to support anything else.
### They must return an ordinary array. No need to propagate the dimnames.
setGeneric("extract_array", signature="x",
    function(x, index) standardGeneric("extract_array")
)

setMethod("extract_array", "ANY",
    function(x, index)
    {
        slice <- subset_by_Nindex(x, index)
        as.array(slice)
    }
)

setMethod("extract_array", "array",
    function(x, index) subset_by_Nindex(x, index)
)

### Equivalent to
###
###     subset_by_Nindex(as.matrix(x), index)
###
### but avoids the cost of turning the full data frame 'x' into a matrix so
### memory footprint stays small when 'index' is small.
setMethod("extract_array", "data.frame",
    function(x, index)
    {
        #ans_type <- .get_data_frame_type(x)
        slice0 <- .extract_data_frame_slice0(x)
        slice <- .extract_data_frame_slice(x, index)
        data <- unlist(c(slice0, slice), use.names=FALSE)
        array(data, dim=get_Nindex_lengths(index, dim(x)))
    }
)

### Equivalent to
###
###     subset_by_Nindex(as.matrix(as.data.frame(x)), index)
###
### but avoids the cost of turning the full DataFrame 'x' first into a data
### frame then into a matrix so memory footprint stays small when 'index' is
### small.
setMethod("extract_array", "DataFrame",
    function(x, index)
    {
        #ans_type <- .get_DataFrame_type(x)
        slice0 <- .extract_DataFrame_slice0(x)
        slice <- .extract_DataFrame_slice(x, index)
        data <- unlist(c(slice0, slice), use.names=FALSE)
        array(data, dim=get_Nindex_lengths(index, dim(x)))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### type() generic and default method
###

setGeneric("type", function(x) standardGeneric("type"))

setMethod("type", "array", function(x) typeof(x))

### type() works out-of-the-box on any array-like object for which
### extract_array() works.
setMethod("type", "ANY",
    function(x)
    {
        index <- rep.int(list(integer(0)), length(dim(x)))
        type(extract_array(x, index))
    }
)

