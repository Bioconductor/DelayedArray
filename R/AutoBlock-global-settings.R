### =========================================================================
### AutoBlock global settings
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### set/getAutoBlockSize()
###
### The automatic block size must be specified in bytes.
###

### We set the automatic block size to 100 Mb by default.
set_auto.block.size <- function(size=1e8)
{
    S4Arrays:::set_user_option("auto.block.size", size)
}

setAutoBlockSize <- function(size=1e8)
{
    if (!isSingleNumber(size) || size < 1)
        stop(wmsg("the block size must be a single number >= 1"))
    prev_size <- S4Arrays:::get_user_option("auto.block.size")
    set_auto.block.size(size)
    message("automatic block size set to ", size, " bytes ",
            "(was ", prev_size, ")")
    invisible(size)
}

getAutoBlockSize <- function()
{
    size <- S4Arrays:::get_user_option("auto.block.size")
    if (!isSingleNumber(size) || size < 1)
        stop(wmsg("DelayedArray user-controlled global option ",
                  "auto.block.size should be a single number >= 1. ",
                  "Fix it with setAutoBlockSize()."))
    size
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### getAutoBlockLength()
###

### The elements of a character vector or a list have a variable size.
### For a character vector: the minimum size of an element is 8 bytes which
### is the overhead of a CHARSXP object. This doesn't account for the string
### data itself.
### For a list: the minimum size of a list element is 8 bytes and is obtained
### when the element is a NULL. However, assuming that a list will typically
### contain more non-NULL than NULL elements and that the non-NULL elements
### will typically be atomic vectors, the average element size is more likely
### to be >= the overhead of an atomic vector which is 56 bytes.
get_type_size <- function(type)
{
    ### Atomic type sizes in bytes.
    TYPE_SIZES <- c(
        logical=4L,
        integer=4L,
        numeric=8L,
        double=8L,
        complex=16L,
        character=8L,  # overhead of a CHARSXP object
        raw=1L,
        list=56L       # overhead of an atomic vector
    )
    if (missing(type))
        return(TYPE_SIZES)
    if (is.factor(type)) {
        type <- as.character(type)
    } else if (!is.character(type)) {
        stop(wmsg("'type' must be a character vector or factor"))
    }
    if (any(type %in% ""))
        stop(wmsg("'type' cannot contain empty strings"))
    idx <- which(!(type %in% c(names(TYPE_SIZES), NA_character_)))
    if (length(idx) != 0L) {
        unsupported_types <- unique(type[idx])
        in1string <- paste0(unsupported_types, collapse=", ")
        stop(wmsg("unsupported type(s): ",  in1string))
    }
    TYPE_SIZES[type]
}

getAutoBlockLength <- function(type)
{
    if (missing(type))
        stop(wmsg("Please specify the type of the array data. ",
                  "See ?getAutoBlockLength"))
    if (!isSingleString(type))
        stop(wmsg("'type' must be a single string"))
    type_size <- get_type_size(type)
    block_size <- getAutoBlockSize()
    ans <- block_size / type_size
    if (ans > .Machine$integer.max)
        stop(wmsg("Automatic block length is too big. Blocks of ",
                  "length > .Machine$integer.max are not supported yet. ",
                  "Please reduce the automatic block length by reducing ",
                  "the automatic block size with setAutoBlockSize()."))
    max(as.integer(ans), 1L)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### set/getAutoBlockShape()
###

SUPPORTED_BLOCK_SHAPES <- c("hypercube",
                            "scale",
                            "first-dim-grows-first",
                            "last-dim-grows-first")

### We set the automatic block shape to "hypercube" by default.
set_auto.block.shape <- function(shape="hypercube")
{
    S4Arrays:::set_user_option("auto.block.shape", shape)
}

setAutoBlockShape <- function(shape=c("hypercube",
                                      "scale",
                                      "first-dim-grows-first",
                                      "last-dim-grows-first"))
{
    shape <- match.arg(shape)
    prev_shape <- S4Arrays:::get_user_option("auto.block.shape")
    set_auto.block.shape(shape)
    message("automatic block shape set to \"", shape, "\" ",
             "(was \"", prev_shape, "\")")
    invisible(shape)
}

getAutoBlockShape <- function()
{
    shape <- S4Arrays:::get_user_option("auto.block.shape")
    if (!(isSingleString(shape) && shape %in% SUPPORTED_BLOCK_SHAPES)) {
        in1string <- paste(paste0("\"", SUPPORTED_BLOCK_SHAPES, "\""),
                           collapse=", ")
        stop(wmsg("DelayedArray user-controlled global option ",
                  "auto.block.shape should be one of: ", in1string, ". ",
                  "Fix it with setAutoBlockShape()."))
    }
    shape
}

