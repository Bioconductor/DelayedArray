.onLoad <- function(libname, pkgname)
{
    options(DelayedArray.simplify=TRUE)
    #user_options_file_exists <- S4Arrays:::user_options_file_exists
    #if (!user_options_file_exists()) {
    #    ## Initialize DelayedArray user controlled global options (setting
    #    ## the 1st option creates the file where the options are stored).
    #    setAutoGridMaker()
    #    set_auto.block.size()
    #    set_auto.block.shape()
    #    setAutoBPPARAM()
    #}
    user_option_is_set <- S4Arrays:::user_option_is_set
    if (!user_option_is_set("auto.grid.maker"))
        setAutoGridMaker()
    if (!user_option_is_set("auto.block.size"))
        set_auto.block.size()
    if (!user_option_is_set("auto.block.shape"))
        set_auto.block.shape()
    if (!user_option_is_set("auto.BPPARAM"))
        setAutoBPPARAM()
    if (!user_option_is_set("auto.mult.parallel.agnostic"))
        setAutoMultParallelAgnostic()
}

.onUnload <- function(libpath)
{
    library.dynam.unload("DelayedArray", libpath)
}

.test <- function()
{
    ## Unit tests temporarily disabled on merida1 and kjohnson3 ...
    slow_build_machine <- function() {
        isTRUE(as.logical(Sys.getenv("IS_BIOC_BUILD_MACHINE"))) &&
            (tolower(Sys.info()[["nodename"]]) %in% c("merida1", "kjohnson3"))
    }
    if (!slow_build_machine()) {
        setAutoRealizationBackend("RleArray")
        on.exit(setAutoRealizationBackend())
        BiocGenerics:::testPackage("DelayedArray")
    }
}

