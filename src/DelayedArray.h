#include <Rdefines.h>


/* compress_atomic_vector.c */

SEXP C_simple_object_size(SEXP x);
SEXP C_encode_atomic_vector(SEXP x);
SEXP C_decode_atomic_vector(SEXP x);
SEXP C_compress_atomic_vector(SEXP x);
SEXP C_decompress_atomic_vector(SEXP x);

/* array_selection.c */

SEXP C_Lindex2Mindex(SEXP Lindex, SEXP dim, SEXP use_names);
SEXP C_Mindex2Lindex(SEXP Mindex, SEXP dim, SEXP use_names, SEXP as_integer);

/* abind.c */

SEXP C_abind(SEXP objects, SEXP nblock, SEXP ans_dim);

/* sparseMatrix_utils.c */

SEXP C_rowsum_dgCMatrix(SEXP x, SEXP group, SEXP ngroup, SEXP na_rm);
SEXP C_colMins_dgCMatrix(SEXP x, SEXP na_rm);
SEXP C_colMaxs_dgCMatrix(SEXP x, SEXP na_rm);
SEXP C_colRanges_dgCMatrix(SEXP x, SEXP na_rm);
SEXP C_colVars_dgCMatrix(SEXP x, SEXP na_rm);

