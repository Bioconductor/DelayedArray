/****************************************************************************
 *                 A Nested Containment List implementation                 *
 *                            Author: H. Pag\`es                            *
 ****************************************************************************/
#include "DelayedArray.h"
#include "S4Vectors_interface.h"


/****************************************************************************
 * 2 low-level helpers to get the length and values of an "extended numeric"
 * vector i.e. of an integer, numeric, or LLint vector
 *
 * TODO: Maybe move them to S4Vectors.
 */

static R_xlen_t get_xnum_length(SEXP x)
{
	if (IS_INTEGER(x) || IS_NUMERIC(x))
		return XLENGTH(x);
	if (!is_LLint(x))
		error("error in get_xnum_length(): "
		      "'x' must be an \"extended numeric\" vector");
	return get_LLint_length(x);
}

static long long int get_xnum_val(SEXP x, R_xlen_t i)
{
	double x_elt;

	if (IS_INTEGER(x))
		return (long long int) INTEGER(x)[i];
	if (is_LLint(x))
		return get_LLint_dataptr(x)[i];
	if (!IS_NUMERIC(x))
		error("error in get_xnum_val(): "
		      "'x' must be an \"extended numeric\" vector");
	x_elt = REAL(x)[i];
	if (x_elt > (double) LLONG_MAX || x_elt < (double) -LLONG_MAX)
		error("error in get_xnum_val(): "
		      "'x[i]' not in the long long int range");
	return (long long int) x_elt;
}


/****************************************************************************
 * abind()
 */

SEXP abind(SEXP objects, SEXP nblock, SEXP ans_dim)
{
	int nobject;
	long long int nblock0, i, j, ans_offset, ans_block_nelt, block_nelt;
	R_xlen_t object_len, ans_len;
	SEXPTYPE ans_type;
	SEXP object, ans, dim;

	if (!isVectorList(objects))  // IS_LIST() is broken
		error("'objects' must be a list");
	nobject = LENGTH(objects);
	if (nobject == 0)
		error("'objects' must contain at least one object");
        if (get_xnum_length(nblock) != 1)
		error("'nblock' must be a single number");
	nblock0 = get_xnum_val(nblock, 0);
	if (nblock0 <= 0)
		error("'nblock' must be > 0");

	/* Determine 'ans_len' and 'ans_type'. */
	ans_type = NILSXP;  // avoid gcc maybe-uninitialized warning
	ans_len = 0;
	for (i = 0; i < nobject; i++) {
		object = VECTOR_ELT(objects, i);
		if (i == 0) {
			ans_type = TYPEOF(object);
		} else if (TYPEOF(object) != ans_type) {
			error("the arrays to bind must have the same type");
		}
		object_len = XLENGTH(object);
		if (object_len % nblock0 != 0)
			error("the arrays to bind must have a length that "
			      "is a multiple of 'nblock'");
		ans_len += object_len;
	}
	ans_block_nelt = ans_len / nblock0;

	/* Alloc and fill 'ans'. */
	ans = PROTECT(allocVector(ans_type, ans_len));
	ans_offset = 0;
	for (i = 0; i < nobject; i++) {
		object = VECTOR_ELT(objects, i);
		object_len = XLENGTH(object);
		block_nelt = object_len / nblock0;
		for (j = 0; j < nblock0; j++) {
			copy_vector_block(ans, ans_offset + j * ans_block_nelt,
				object, j * block_nelt,
				block_nelt);
		}
		ans_offset += block_nelt;
	}

	/* Set "dim" attribute on 'ans'. */
	dim = PROTECT(duplicate(ans_dim));
	SET_DIM(ans, dim);
	UNPROTECT(2);
	return ans;
}

