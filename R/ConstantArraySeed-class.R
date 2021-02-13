setClass("ConstantArraySeed", slots=c(dim='integer', value="ANY"))

setValidity("ConstantArraySeed", function(object) {
    msg <- character(0)

    if (length(dim(object)) < 1 || any(dim(object) < 0)) {
        msg <- c(msg, "'dim' must be an integer vector of length > 1 with non-negative values")
    }

    if (length(object@value) != 1 || !is.atomic(object@value)) {
        msg <- c(msg, "'value' must be a single atomic value")
    }

    if (length(msg)) {
        return(msg)
    }
    TRUE
})

.get_constant_dim <- function(x, index) {
    odims <- dim(x)
    for (i in seq_along(odims)) {
        if (!is.null(index[[i]])) {
            odims[i] <- length(index[[i]])
        }
    }
    odims
}

setMethod("extract_array", "ConstantArraySeed", function(x, index) {
    array(x@value, .get_constant_dim(x, index))
})

setMethod("extract_sparse_array", "ConstantArraySeed", function(x, index) {
    SparseArraySeed(.get_constant_dim(x, index), nzdata=rep(x@value, 0))
})

setMethod("is_sparse", "ConstantArraySeed", function(x) isTRUE(x@value==0))

ConstantArraySeed <- function(dim, value) {
    new("ConstantArraySeed", dim=as.integer(dim), value=value)
}
