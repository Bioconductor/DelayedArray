.onLoad <- function(libname, pkgname)
{
    options(DelayedArray.block.size=DEFAULT_BLOCK_SIZE)
}

.test <- function()
{
    cat("------------------------------------------------------------------\n")
    cat("Running tests with realization backend set to \"RleArray\" ...\n")
    setRealizationBackend("RleArray")
    BiocGenerics:::testPackage("DelayedArray")

    cat("\n")
    cat("------------------------------------------------------------------\n")
    cat("Running tests with realization backend set to \"HDF5Array\" ...\n")
    setRealizationBackend("HDF5Array")
    BiocGenerics:::testPackage("DelayedArray")
}

