#include "DelayedArray.h"
#include <R_ext/Rdynload.h>

#define CALLMETHOD_DEF(fun, numArgs) {#fun, (DL_FUNC) &fun, numArgs}

static const R_CallMethodDef callMethods[] = {

/* abind.c */
	CALLMETHOD_DEF(abind, 3),

	{NULL, NULL, 0}
};


void R_init_DelayedArray(DllInfo *info)
{
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);

	return;
}

