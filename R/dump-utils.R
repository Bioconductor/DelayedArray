### =========================================================================
### Dump management utilities
### -------------------------------------------------------------------------
###

### Virtual class with no slots. Intended to be extended by implementations
### of DelayedArray on-disk backends. Concrete subclasses must implement:
###   1) A constructor function that takes argument 'dim', 'dimnames', and
###      'type'.
###   2) A "write_to_dump" method that works on an ordinary array.
###   3) A "close" method.
###   4) Coercion to DelayedArray.
### See HDF5ArrayDump class in the HDF5Array package for an example.
setClass("OnDiskArrayDump", representation("VIRTUAL"))

setGeneric("write_to_dump", signature=c("x", "dump"),
    function(x, dump, offsets=NULL) standardGeneric("write_to_dump")
)

setGeneric("close")

### Default "close" method for OnDiskArrayDump objects. A default
### "write_to_dump" method is defined in DelayedArray-class.R.
setMethod("close", "OnDiskArrayDump",
    function(con)
        stop(wmsg("don't know how to close a ", class(con), " object"))
)

