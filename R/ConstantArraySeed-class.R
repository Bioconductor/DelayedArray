setClass("ConstantArraySeed", slots=c(dim='integer', value="vector"))

setValidity2("ConstantArraySeed", function(object) {
    msg <- validate_dim_slot(object, "dim")
    if (!isTRUE(msg)) {
        return(msg)
    }
    if (length(object@value) != 1) {
        return("'value' must be a vector of length 1")
    }
    TRUE
})

setMethod("extract_array", "ConstantArraySeed", function(x, index) {
    array(x@value, get_Nindex_lengths(index, dim(x)))
})

setMethod("extract_sparse_array", "ConstantArraySeed", function(x, index) {
    SparseArraySeed(get_Nindex_lengths(index, dim(x)), nzdata=rep(x@value, 0L))
})

setMethod("is_sparse", "ConstantArraySeed", function(x) {
    zero <- vector(type(x), length=1L)
    identical(x@value, zero)
})

ConstantArraySeed <- function(dim, value) {
    new("ConstantArraySeed", dim=as.integer(dim), value=value)
}
