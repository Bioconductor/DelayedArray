### =========================================================================
### Bind DelayedArray objects along their rows or columns
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### arbind() and acbind()
###

.DelayedArray_arbind <- function(...)
{
    objects <- unname(list(...))
    dims <- IRanges:::get_dims_to_bind(objects, 1L)
    if (is.character(dims))
        stop(wmsg(dims))
    DelayedArray(new_SeedBinder(objects, 1L))
}

.DelayedArray_acbind <- function(...)
{
    objects <- unname(list(...))
    dims <- IRanges:::get_dims_to_bind(objects, 2L)
    if (is.character(dims))
        stop(wmsg(dims))
    DelayedArray(new_SeedBinder(objects, 2L))
}

setMethod("arbind", "DelayedArray", .DelayedArray_arbind)
setMethod("acbind", "DelayedArray", .DelayedArray_acbind)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### rbind() and cbind()
###

setMethod("rbind", "DelayedMatrix", .DelayedArray_arbind)
setMethod("cbind", "DelayedMatrix", .DelayedArray_acbind)

.as_DelayedMatrix_objects <- function(objects)
{
    lapply(objects,
        function(object) {
            if (length(dim(object)) != 2L)
                stop(wmsg("cbind() and rbind() are not supported on ",
                          "DelayedArray objects that don't have exactly ",
                          "2 dimensions. Please use acbind() or arnind() ",
                          "instead."))
            as(object, "DelayedMatrix")
        })
}

.DelayedArray_rbind <- function(...)
{
    objects <- .as_DelayedMatrix_objects(list(...))
    do.call("rbind", objects)
}

.DelayedArray_cbind <- function(...)
{
    objects <- .as_DelayedMatrix_objects(list(...))
    do.call("cbind", objects)
}

setMethod("rbind", "DelayedArray", .DelayedArray_rbind)
setMethod("cbind", "DelayedArray", .DelayedArray_cbind)

