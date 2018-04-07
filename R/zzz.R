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
    cat("------------------------------------------------------------------\n")
    cat("Running tests with realization backend set to \"RleArray\" ...\n")
    setRealizationBackend("RleArray")
    BiocGenerics:::testPackage("DelayedArray")
    }

    cat("\n")
    cat("------------------------------------------------------------------\n")
    cat("Running tests with realization backend set to \"HDF5Array\" ...\n")
    setRealizationBackend("HDF5Array")
    BiocGenerics:::testPackage("DelayedArray")
}

