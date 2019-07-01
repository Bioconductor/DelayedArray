### =========================================================================
### Some low-level utilities
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


### TODO: This should probably go to S4Vectors (but maybe find a better name
### for it first).
seq2 <- function(to, by)
{
    stopifnot(isSingleNumber(to), isSingleNumber(by))
    ans <- seq_len(to %/% by) * by
    if (to %% by != 0L)
        ans <- c(ans, to)
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Greatest common divisor and least common multiple of 2 integers
###
### Currently not used.
###
### TODO: Maybe should go to S4Vectors (and be implemented in C).
###

if (FALSE) {  # ------------ BEGIN DISABLED CODE ------------

### Greatest common divisor of 2 integers.
### Can be about 10x faster than this recursive (and vectorized) solution
### from SO:
###   https://stackoverflow.com/questions/21502181/finding-the-gcd-without-looping-r
### Also check user input and handle edge cases.
GCD <- function(x, y)
{
    stopifnot(isSingleNumber(x), isSingleNumber(y))
    if (!is.integer(x))
        x <- as.integer(x)
    if (!is.integer(y))
        y <- as.integer(y)
    if (y == 0L) {
        if (x == 0L)
            return(NA_integer_)
        return(x)
    }
    while (TRUE) {
        r <- x %% y
        if (r == 0L)
            return(y)
        x <- y
        y <- r
    }
}

### Least common multiple of 2 integers.
LCM <- function(x, y)
{
    stopifnot(isSingleNumber(x), isSingleNumber(y))
    if (!is.integer(x))
        x <- as.integer(x)
    if (!is.integer(y))
        y <- as.integer(y)
    gcd <- GCD(x, y)
    (x %/% gcd) * y  # could overflow
}
}             # ------------- END DISABLED CODE -------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### normarg_dim() and normarg_dimnames()
###

normarg_dim <- function(dim)
{
    if (!is.numeric(dim))
        stop(wmsg("'dim' must be an integer vector"))
    if (!is.integer(dim))
        dim <- as.integer(dim)
    if (length(dim) == 0L)
        stop(wmsg("'dim' cannot be an empty vector"))
    if (S4Vectors:::anyMissingOrOutside(dim, 0L))
        stop(wmsg("'dim' cannot contain negative or NA values"))
    dim
}

normarg_dimnames <- function(dimnames, dim)
{
    ndim <- length(dim)
    if (missing(dimnames) || is.null(dimnames))
        return(vector("list", length=ndim))
    if (!is.list(dimnames))
        stop(wmsg("'dimnames' must be NULL or a list"))
    if (length(dimnames) != ndim)
        stop(wmsg("'dimnames' must have one list element per dimension"))
    dimnames
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### 2 wrappers to dim<- and dimnames<- that try to avoid unnecessary copies
### of 'x'
###

set_dim <- function(x, value)
{
    if (!identical(dim(x), value))
        dim(x) <- value
    x
}

set_dimnames <- function(x, value)
{
    if (!identical(dimnames(x), value))
        dimnames(x) <- value
    x
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### get_first_non_NULL_dimnames()
###

### Implement (and extend to the N-ary case) propagation of the dimnames
### following the "first array with dimnames wins" rule that seems to be
### used by all the binary operations from the Arith (e.g. "+", "^", etc...),
### Compare (e.g. "==", ">=", etc...) and Logic ("&", "|") groups.
### Only makes sense to use if all the array-like objects in 'objects'
### have the same dimensions.
### Note that rbind() and cbind() use a different rule for propagation of
### the dimnames. This rule is implemented (and extended to the N-ary case)
### by combine_dimnames_along() defined in bind-arrays.R
get_first_non_NULL_dimnames <- function(objects)
{
    for (object in objects) {
        object_dimnames <- dimnames(object)
        if (!is.null(object_dimnames))
            return(object_dimnames)
    }
    NULL
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### simplify_NULL_dimnames()
###

simplify_NULL_dimnames <- function(dimnames)
{
    if (all(S4Vectors:::sapply_isNULL(dimnames)))
        return(NULL)
    dimnames
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### combine_array_objects()
###
### NOT exported but used in unit tests.
###

### 'objects' must be a list of array-like objects that support as.vector().
combine_array_objects <- function(objects)
{
    if (!is.list(objects))
        stop("'objects' must be a list")
    NULL_idx <- which(S4Vectors:::sapply_isNULL(objects))
    if (length(NULL_idx) != 0L)
        objects <- objects[-NULL_idx]
    if (length(objects) == 0L)
        return(NULL)
    unlist(lapply(objects, as.vector), recursive=FALSE, use.names=FALSE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Used in validity methods
###

### NOT exported but used in the HDF5Array package!
validate_dim_slot <- function(x, slotname="dim")
{
    x_dim <- slot(x, slotname)
    if (!is.integer(x_dim))
        return(paste0("'", slotname, "' slot must be an integer vector"))
    if (length(x_dim) == 0L)
        return(paste0("'", slotname, "' slot cannot be empty"))
    if (S4Vectors:::anyMissingOrOutside(x_dim, 0L))
        return(paste0("'", slotname, "' slot cannot contain negative ",
                      "or NA values"))
    TRUE
}

validate_dimnames_slot <- function(x, dim, slotname="dimnames")
{
    x_dimnames <- slot(x, slotname)
    if (!is.list(x_dimnames))
        return(paste0("'", slotname, "' slot must be a list"))
    if (length(x_dimnames) != length(dim))
        return(paste0("'", slotname, "' slot must have ",
                      "one list element per dimension in the object"))
    ok <- vapply(seq_along(dim),
                 function(along) {
                   dn <- x_dimnames[[along]]
                   if (is.null(dn))
                       return(TRUE)
                   is.character(dn) && length(dn) == dim[[along]]
                 },
                 logical(1),
                 USE.NAMES=FALSE)
    if (!all(ok))
        return(paste0("each list element in '", slotname, "' slot ",
                      "must be NULL or a character vector along ",
                      "the corresponding dimension in the object"))
    TRUE
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### outer2()
###
### Outer product of an arbitrary number of objects. Typically used on
### objects that are vectors and/or arrays.
###

outer2 <- function(objects, FUN="*", ...)
{
    stopifnot(is.list(objects), length(objects) != 0L)
    nobject <- length(objects)
    ans <- objects[[1L]]
    if (nobject >= 2L) {
        for (i in 2:nobject) {
            idx <- objects[[i]]
            ans <- outer(ans, idx, FUN=FUN, ...)
        }
    }
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Translate an index into the whole to an index into the parts
###
### This is .rowidx2rowkeys() from BSgenome/R/OnDiskLongTable-class.R, copied
### here and renamed get_part_index() !
### TODO: Put it somewhere else where it can be shared.
###

.breakpoints2offsets <- function(breakpoints)
{
    breakpoints_len <- length(breakpoints)
    if (breakpoints_len == 0L)
        return(integer(0))
    c(0L, breakpoints[-breakpoints_len])
}

### breakpoints: integer vector of break points that define the parts.
### idx:         index into the whole as an integer vector.
### Return a list of 2 integer vectors parallel to 'idx'. The 1st vector
### contains part numbers and the 2nd vector indices into the parts.
### In addition, if 'breakpoints' has names (part names) then they are
### propagated to the 1st vector.
get_part_index <- function(idx, breakpoints)
{
    part_idx <- findInterval(idx, breakpoints + 1L) + 1L
    names(part_idx) <- names(breakpoints)[part_idx]
    rel_idx <- idx - .breakpoints2offsets(unname(breakpoints))[part_idx]
    list(part_idx, rel_idx)
}

split_part_index <- function(part_index, npart)
{
    ans <- rep.int(list(integer(0)), npart)
    tmp <- split(unname(part_index[[2L]]), part_index[[1L]])
    ans[as.integer(names(tmp))] <- tmp
    ans
}

get_rev_index <- function(part_index)
{
    f <- part_index[[1L]]
    idx <- split(seq_along(f), f)
    idx <- unlist(idx, use.names=FALSE)
    rev_idx <- integer(length(idx))
    rev_idx[idx] <- seq_along(idx)
    rev_idx
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### set_user_option() / get_user_option()
###

if (FALSE) {  # ------------ BEGIN DISABLED CODE ------------

### In the context of BiocParallel::bplapply() and family, we want the workers
### to inherit the user-controlled options defined on the master. Workers can
### also modify the inherited options or define new user-controlled options
### but this should not affect the options defined on the master (one-way
### communication). This mechanism should also work with nested levels of
### parallelization i.e. if workers start subworkers and so on.
### The mechanism implemented below achieves that. It relies on 2 important
### assumptions:
### 1st assumption: Environment variables defined in the environment of a
###                 master R session are **always** inherited by the workers
###                 started by BiocParallel::bplapply() and family. "Always"
###                 here means "whatever the BPPARAM backend is".
### 2nd assumption: The workers can access the tempdir() of their master.
###                 WARNING: It has been reported that with some cluster
###                 configurations the nodes don't have access to the
###                 tempdir() of the head. So this is NOT a safe assumption!

.make_new_filepath_for_my_own_file <- function(my_pid)
{
    filename <- paste0("DelayeArray-options.", my_pid)
    file.path(tempdir(), filename)
}

.filepath_is_mine <- function(filepath, my_pid)
{
    filename <- basename(filepath)
    if (!grepl(".", filename, fixed=TRUE))
        return(FALSE)
    ## We extract the PID embedded in the filename by grabbing everything
    ## that is after the last occurence of '.'.
    pid_from_filename <- sub("(.*)\\.([^.]*)$", "\\2", filename)
    pid_from_filename == my_pid
}

### Create a file with no options in it.
.create_root_file <- function(my_pid)
{
    new_filepath <- .make_new_filepath_for_my_own_file(my_pid)
    saveRDS(list(), new_filepath)
    new_filepath
}

### If the 2nd assumption (see above) turns out to be incorrect then we want
### to fail early and graciously. So perform the check below before trying
### to read or copy the file located at 'filepath'.
.check_2nd_assumption <- function(filepath)
{
    if (!file.exists(filepath))
        stop(wmsg("worker cannot access file ", filepath,
                  " created by master"))
}

.copy_file_if_not_mine <- function(filepath, my_pid)
{
    is_mine <- .filepath_is_mine(filepath, my_pid)
    if (is_mine)
        return(filepath)
    ## If the file is not mine then the current process is a worker that
    ## was started by a master, so I'm not allowed to write to the file.
    ## I need to make my own copy it first.
    .check_2nd_assumption(filepath)
    new_filepath <- .make_new_filepath_for_my_own_file(my_pid)
    ok <- file.copy(filepath, new_filepath)
    if (!ok)
        stop(wmsg("worker was not able to copy file ", filepath,
                  " created by master"))
    new_filepath
}

### Name of the environment variable that contains the path to the RDS file
### containing the serialized options.
.env_var_name <- "DelayedArray_OPTIONS_FILEPATH"

.get_user_options_filepath <- function() Sys.getenv(.env_var_name)

.set_user_options_filepath <- function(filepath)
{
    do.call(Sys.setenv, setNames(list(filepath), .env_var_name))
}

user_options_file_exists <- function() .get_user_options_filepath() != ""

### Implements a copy-on-write mechanism in the context of multiple workers.
set_user_option <- function(name, value)
{
    stopifnot(isSingleString(name))
    filepath <- .get_user_options_filepath()
    my_pid <- as.character(Sys.getpid())
    if (filepath == "") {
        ## Options file doesn't exist yet. This means that the current
        ## process is the root of the (possibly nested) master/workers tree,
        ## that is, the R session explicitely started by the user at the
        ## command line (or via RStudio).
        writable_filepath <- .create_root_file(my_pid)
    } else {
        ## The current process is a worker (which in turn could itself
        ## become the master of some subworkers at some point).
        writable_filepath <- .copy_file_if_not_mine(filepath, my_pid)
    }
    if (writable_filepath != filepath)
        .set_user_options_filepath(writable_filepath)
    user_options <- readRDS(writable_filepath)
    ## We use 'x[name] <- list(value)' instead of 'x[[name]] <- value'
    ## because the latter removes the list element if value is NULL and
    ## we don't want that.
    user_options[name] <- list(value)
    saveRDS(user_options, writable_filepath)
    invisible(value)
}

get_user_options <- function()
{
    ## Only 2 possibilities: either the current process "owns" 'filepath'
    ## or not. In the latter case it means that the current process is a
    ## worker that was started by a master.
    filepath <- .get_user_options_filepath()
    .check_2nd_assumption(filepath)
    user_options <- try(readRDS(filepath), silent=TRUE)
    if (inherits(user_options, "try-error"))
        stop(wmsg("worker was not able to read file ", filepath,
                  " created by master"))
    user_options
}

get_user_option <- function(name)
{
    stopifnot(isSingleString(name))
    user_options <- get_user_options()
    ## 'user_options' should be a list.
    idx <- match(name, names(user_options))
    if (is.na(idx))
        stop(wmsg("Unkown DelayedArray user-controlled global option: ", name))
    user_options[[idx]]
}
}             # ------------- END DISABLED CODE -------------

set_user_option <- function(name, value)
{
    stopifnot(isSingleString(name))
    name <- paste0("DelayedArray.", name)
    options(setNames(list(value), name))
    invisible(value)
}

get_user_option <- function(name)
{
    stopifnot(isSingleString(name))
    name <- paste0("DelayedArray.", name)
    getOption(name)
}

user_option_is_set <- function(name)
{
    stopifnot(isSingleString(name))
    name <- paste0("DelayedArray.", name)
    name %in% names(.Options)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### bplapply2()
###
### A simple wrapper to BiocParallel::bplapply() that propagates
### base::.Options to the workers.
### Should probably go in BiocParallel, or even better, replace bplapply().
###

bplapply2 <- function(X, FUN, ..., BPREDO=list(), BPPARAM=bpparam())
{
    FUN2 <- function(x, opts)
    {
        force(opts)
        base::options(opts)
        FUN(x, ...)
    }
    bplapply(X, FUN2, opts=base::.Options, BPREDO=BPREDO, BPPARAM=BPPARAM)
}

