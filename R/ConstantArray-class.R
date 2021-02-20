### =========================================================================
### ConstantArraySeed and ConstantArray objects
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### ConstantArraySeed objects
###

setClass("ConstantArraySeed",
    contains="Array",
    representation(
        dim="integer",  # This gives us dim() for free!
        value="vector"
    ),
    prototype(
        dim=0L,
        value=NA
    )
)

setValidity2("ConstantArraySeed",
    function(object)
    {
        msg <- validate_dim_slot(object, "dim")
        if (!isTRUE(msg))
            return(msg)
        if (length(object@value) != 1L)
            return("'value' must be a vector (atomic or list) of length 1")
        TRUE
    }
)

setMethod("extract_array", "ConstantArraySeed",
    function(x, index) array(x@value, get_Nindex_lengths(index, dim(x)))
)

setMethod("extract_sparse_array", "ConstantArraySeed",
    function(x, index)
    {
        ans_dim <- get_Nindex_lengths(index, dim(x))
        ans_nzdata <- rep.int(x@value, 0L)
        SparseArraySeed(ans_dim, nzdata=ans_nzdata, check=FALSE)
    }
)

setMethod("is_sparse", "ConstantArraySeed",
    function(x)
    {
        zero <- vector(type(x), length=1L)
        identical(x@value, zero)
    }
)

ConstantArraySeed <- function(dim, value=NA)
{
    new2("ConstantArraySeed", dim=as.integer(dim), value=value)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### ConstantArray and ConstantMatrix objects
###

setClass("ConstantArray",
    contains="DelayedArray",
    representation(seed="ConstantArraySeed")
)

setClass("ConstantMatrix", contains=c("ConstantArray", "DelayedMatrix"))

setMethod("matrixClass", "ConstantArray", function(x) "ConstantMatrix")

setMethod("DelayedArray", "ConstantArraySeed",
    function(seed) new_DelayedArray(seed, Class="ConstantArray")
)

ConstantArray <- function(dim, value=NA)
{
    DelayedArray(ConstantArraySeed(dim, value=value))
}

### Automatic coercion method from ConstantArray to ConstantMatrix silently
### returns a broken object (unfortunately these dummy automatic coercion
### methods don't bother to validate the object they return). So we overwrite
### it.
setAs("ConstantArray", "ConstantMatrix",
    function(from) new2("ConstantMatrix", from)
)

### The user should not be able to degrade a ConstantMatrix object to
### a ConstantArray object so 'as(x, "ConstantArray", strict=TRUE)' should
### fail or be a no-op when 'x' is ConstantMatrix object. Making this
### coercion a no-op seems to be the easiest (and safest) way to go.
setAs("ConstantMatrix", "ConstantArray", function(from) from)  # no-op

