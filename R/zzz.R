.onLoad <- function(libname, pkgname)
{
    options(DelayedArray.block.size=DEFAULT_BLOCK_SIZE)
}

.test <- function() BiocGenerics:::testPackage("DelayedArray")

