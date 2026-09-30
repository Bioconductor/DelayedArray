### =========================================================================
### Bind DelayedArray objects with abind()
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### abind()
###
### Note that making abind() work on DelayedArray objects also makes arbind(),
### acbind(), rbind(), cbind(), and bindROWS() work on these objects.
###

.abind_DelayedArray_objects <- function(..., along=NULL, rev.along=NULL)
{
    objects <- S4Arrays:::homogenize_abind_input_objects(list(...),
                                           along=along, rev.along=rev.along)
    if (length(objects) == 0L)
        return(NULL)

    x <- objects[[1L]]
    ## Note that unary abind() is not necessarily a no-op because
    ## homogenize_abind_input_objects() can add dimensions to the input object.
    if (length(objects) == 1L)
        return(x)

    stash_DelayedAbind(x, objects[-1L], along=attr(objects, "along"))
}

setMethod("abind", "DelayedArray", .abind_DelayedArray_objects)

