.onLoad <- function(libname, pkgname)
{
    options(DelayedArray.simplify=TRUE)
    #if (!user_options_file_exists()) {
        ## Initialize DelayedArray user controlled global options (setting
        ## the 1st option creates the file where the options are stored).
        setAutoGridMaker()
        set_auto.block.size()
        set_auto.block.shape()
        setAutoBPPARAM()
    #}
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

