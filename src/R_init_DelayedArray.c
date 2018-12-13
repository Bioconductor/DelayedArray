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

/* abind.c */
	CALLMETHOD_DEF(C_abind, 3),

	{NULL, NULL, 0}
};


void R_init_DelayedArray(DllInfo *info)
{
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);

	return;
}

