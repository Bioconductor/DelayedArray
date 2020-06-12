#include "DelayedArray.h"
#include <R_ext/Rdynload.h>

#define CALLMETHOD_DEF(fun, numArgs) {#fun, (DL_FUNC) &fun, numArgs}

static const R_CallMethodDef callMethods[] = {

/* compress_atomic_vector.c */
	CALLMETHOD_DEF(C_simple_object_size, 1),
	CALLMETHOD_DEF(C_encode_atomic_vector, 1),
	CALLMETHOD_DEF(C_decode_atomic_vector, 1),
	CALLMETHOD_DEF(C_compress_atomic_vector, 1),
	CALLMETHOD_DEF(C_decompress_atomic_vector, 1),

/* array_selection.c */
	CALLMETHOD_DEF(C_Lindex2Mindex, 3),
	CALLMETHOD_DEF(C_Mindex2Lindex, 4),

/* abind.c */
	CALLMETHOD_DEF(C_abind, 3),

/* dgCMatrix_stats.c */
	CALLMETHOD_DEF(C_rowsum_dgCMatrix, 4),
	CALLMETHOD_DEF(C_colMins_dgCMatrix, 2),
	CALLMETHOD_DEF(C_colMaxs_dgCMatrix, 2),
	CALLMETHOD_DEF(C_colRanges_dgCMatrix, 2),
	CALLMETHOD_DEF(C_colVars_dgCMatrix, 2),

	{NULL, NULL, 0}
};


void R_init_DelayedArray(DllInfo *info)
{
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);

	return;
}

