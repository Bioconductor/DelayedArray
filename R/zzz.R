.onLoad <- function(libname, pkgname)
{
    options(DelayedArray.simplify=TRUE)
    setAutoGridMaker()
    suppressMessages(setAutoBlockSize())
    suppressMessages(setAutoBlockShape())
    setAutoBPPARAM()
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

