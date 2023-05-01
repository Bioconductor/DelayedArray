/****************************************************************************
 *              Compression/decompression of an atomic vector               *
 *                            Author: H. Pag\`es                            *
 ****************************************************************************/
#include "DelayedArray.h"


/* TODO: This stuff has nothing to do here! Get rid of it or move it to
   S4Vectors. */

/****************************************************************************
 * C_simple_object_size()
 */

#define FIXED_SIZE	48  /* object.size(raw(0)) or object.size(list()) */

/* Only supports naked atomic vectors, possibly wrapped in a naked list. */
static size_t simple_object_size(SEXP x)
{
	R_xlen_t x_len;
	size_t s;

	if (ATTRIB(x) != R_NilValue)
		error("attributes not supported by simple_object_size()");
	x_len = XLENGTH(x);
	s = (size_t) x_len;
	switch (TYPEOF(x)) {
	    case LGLSXP: case INTSXP:
		s *= sizeof(int);
	    break;
	    case REALSXP:
		s *= sizeof(double);
	    break;
	    case CPLXSXP:
		s *= sizeof(Rcomplex);
	    break;
	    case RAWSXP:
	    break;
	    case VECSXP:
		s *= sizeof(SEXP);
		for (R_xlen_t i = 0; i < x_len; i++)
			s += simple_object_size(VECTOR_ELT(x, i));
	    break;
	    default:
		error("object of type %s not supported "
		      "by simple_object_size()",
		      CHAR(type2str(TYPEOF(x))));
	}
	return s + FIXED_SIZE;
}

SEXP C_simple_object_size(SEXP x)
{
	size_t s;

	s = simple_object_size(x);
	return s <= INT_MAX ? ScalarInteger((int) s) : ScalarReal((double) s);
}


/****************************************************************************
 * C_encode_atomic_vector() / C_decode_atomic_vector()
 */

SEXP C_encode_atomic_vector(SEXP x)
{
	return R_NilValue;
}

SEXP C_decode_atomic_vector(SEXP x)
{
	return R_NilValue;
}


/****************************************************************************
 * C_compress_atomic_vector() / C_decompress_atomic_vector()
 */

SEXP C_compress_atomic_vector(SEXP x)
{
	return R_NilValue;
}

SEXP C_decompress_atomic_vector(SEXP x)
{
	return R_NilValue;
}

