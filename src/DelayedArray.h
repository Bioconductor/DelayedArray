#include <Rdefines.h>


/* compress_atomic_vector.c */

SEXP C_simple_object_size(SEXP x);
SEXP C_encode_atomic_vector(SEXP x);
SEXP C_decode_atomic_vector(SEXP x);
SEXP C_compress_atomic_vector(SEXP x);
SEXP C_decompress_atomic_vector(SEXP x);

/* abind.c */

SEXP C_abind(SEXP objects, SEXP nblock, SEXP ans_dim);

