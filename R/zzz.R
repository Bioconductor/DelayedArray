.onLoad <- function(libname, pkgname)
{
    options(DelayedArray.simplify=TRUE)
    options(DelayedArray.block.size=DEFAULT_BLOCK_SIZE)
}

.test <- function()
{
    ## Skip this on Windows to avoid 'R CMD check' TIMEOUT on the Windows
    ## build machines.
    if (.Platform$OS.type != "windows") {
        setRealizationBackend("RleArray")
        BiocGenerics:::testPackage("DelayedArray")
    }
}

