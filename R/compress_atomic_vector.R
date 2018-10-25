### =========================================================================
### compress_atomic_vector()
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###

.RAW_OVERHEAD <- object.size(raw(0L))  # 48 bytes on a 64-bit system

### Word size should be 8 on a 64-bit system and 4 on a 32-bit system.
.WORD_SIZE <- object.size(raw(1L)) - .RAW_OVERHEAD

### Predict the size of 'raw(length)' without creating the vector (which can
### be expensive).
### Note that the predicted size is almost always accurate (i.e. almost
### always equal to 'object.size(raw(length))') except for some vector
### lengths for which the raw() constructor strangely tends to "over allocate"
### memory. For example, on my laptop (64-bit Ubuntu 16.04.5 LTS, R 3.5.1)
### 'object.size(raw(64))' is 112, as predicted, but 'object.size(raw(65))'
### is 176 so more than the predicted 120 bytes.
.predict_raw_size <- function(length)
    .RAW_OVERHEAD + ((length - 1L) %/% .WORD_SIZE + 1L) * .WORD_SIZE

.LIST2_OVERHEAD <- object.size(list(NULL, NULL))  # 64 bytes on a 64-bit system

.predict_compressed_size <- function(x_len, ux)
    .LIST2_OVERHEAD + .predict_raw_size(x_len) + object.size(ux)

### If encoding is successful, the compressed vector is returned as a list
### with 2 elements:
###   1) A raw vector of the same length as input vector 'x'.
###   2) A vector of levels (obtained with 'unique(x)').
### Will only compress if encoding with a raw vector is possible (i.e. if
### 'x' contains <= 256 unique values) and worth it (i.e. if compressing
### reduces the object size). Otherwise the input vector is returned
### unmodified (no-op).
compress_atomic_vector <- function(x)
{
    stopifnot(is.atomic(x))
    ## Should we compress?
    if (is.raw(x))
        return(x)  # nothing to gain
    x_len <- length(x)
    if (x_len <= 2L)
        return(x)  # not worth it
    ux <- unique(x)
    ux_len <- length(ux)
    if (ux_len > 256L)
        return(x)  # not possible
    if (.predict_compressed_size(x_len, ux) >= object.size(x))
        return(x)  # not worth it
    ## Encoding is possible and worth it (based on predicted size).
    ans <- list(as.raw(match(x, ux) - 1L), ux)
    if (object.size(ans) >= object.size(x))
        return(x)  # not worth it
    ans
}

### 'decompress_atomic_vector(compress_atomic_vector(x))' should be a no-op
### for any unnamed atomic vector 'x'.
decompress_atomic_vector <- function(x)
{
    if (is.atomic(x))
        return(x)
    stopifnot(is.list(x), length(x) == 2L)
    x[[2L]][as.integer(x[[1L]]) + 1L]
}

