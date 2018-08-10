.onLoad <- function(libname, pkgname)
{
    options(DelayedArray.simplify=TRUE)
    setDefaultGridMaker()
    suppressMessages(setDefaultBlockSize())
    suppressMessages(setDefaultBlockShape())
    ## By default bplapply() uses the SnowParam() backend on Windows, which
    ## introduces **a lot** of overhead in the context of block processing.
    ## See https://github.com/Bioconductor/BiocParallel/issues/78
    ## So we disable parallel evaluation by default on this platform, by
    ## setting the default backend to SerialParam().
    if (.Platform$OS.type == "windows")
        BiocParallel::register(SerialParam())
}

.onUnload <- function(libpath)
{
    library.dynam.unload("DelayedArray", libpath)
}

.test <- function()
{
    setRealizationBackend("RleArray")
    BiocGenerics:::testPackage("DelayedArray")
}

