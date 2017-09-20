### =========================================================================
### subset_seed_as_array()
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
### subset_seed_as_array() generic and methods
###

### 'index' is expected to be an unnamed list of subscripts as positive integer
### vectors, one vector per seed dimension. *Missing* list elements are allowed
### and represented by NULLs.
### The "subset_seed_as_array" methods don't need to support anything else.
### They must return an ordinary array. No need to propagate the dimnames.
setGeneric("subset_seed_as_array", signature="seed",
    function(seed, index) standardGeneric("subset_seed_as_array")
)

setMethod("subset_seed_as_array", "ANY",
    function(seed, index)
    {
        slice <- subset_by_Nindex(seed, index)
        as.array(slice)
    }
)

setMethod("subset_seed_as_array", "array",
    function(seed, index)
        subset_by_Nindex(seed, index)
)

### Equivalent to
###
###     subset_by_Nindex(as.matrix(x), index)
###
### but avoids the cost of turning the full data frame 'x' into a matrix so
### memory footprint stays small when 'index' is small.
setMethod("subset_seed_as_array", "data.frame",
    function(seed, index)
    {
        #ans_type <- .get_data_frame_type(seed)
        slice0 <- .extract_data_frame_slice0(seed)
        slice <- .extract_data_frame_slice(seed, index)
        data <- unlist(c(slice0, slice), use.names=FALSE)
        array(data, dim=get_Nindex_lengths(index, dim(seed)))
    }
)

### Equivalent to
###
###     subset_by_Nindex(as.matrix(as.data.frame(x)), index)
###
### but avoids the cost of turning the full DataFrame 'x' first into a data
### frame then into a matrix so memory footprint stays small when 'index' is
### small.
setMethod("subset_seed_as_array", "DataFrame",
    function(seed, index)
    {
        #ans_type <- .get_DataFrame_type(seed)
        slice0 <- .extract_DataFrame_slice0(seed)
        slice <- .extract_DataFrame_slice(seed, index)
        data <- unlist(c(slice0, slice), use.names=FALSE)
        array(data, dim=get_Nindex_lengths(index, dim(seed)))
    }
)

