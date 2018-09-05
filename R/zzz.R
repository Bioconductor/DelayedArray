.onLoad <- function(libname, pkgname)
{
    options(DelayedArray.simplify=TRUE)
    setDefaultGridMaker()
    suppressMessages(setDefaultBlockSize())
    suppressMessages(setDefaultBlockShape())
    setDefaultBPPARAM()
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

