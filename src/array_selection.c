/****************************************************************************
 *                     Manipulation of array selections                     *
 *                            Author: H. Pag\`es                            *
 ****************************************************************************/
#include "DelayedArray.h"
#include "S4Vectors_interface.h"

#include <limits.h>  /* for INT_MAX, LLONG_MAX, LLONG_MIN */

/*
  An array selection is just an index into an array-like object that defines
  a subset of the array elements. This index can take various forms but 3
  special forms are particularly useful and extensively used thoughout
  the DelayedArray framework:

  1. Linear index (also called "Lindex"): A numeric vector with no NAs
     where each value is >= 1 and <= the length of the array-like object.
     When using an Lindex to subset an array-like object, the returned
     value is a vector-like object (i.e. no dimensions) of the same length
     as the Lindex.

       Example:
         a <- array(101:124, 4:2)
         Lindex <- c(7, 2, 24, 2)
         a[Lindex]

  2. Matrix index (also called "Mindex"): An integer matrix with one column
     per dimension in the array-like object and one row per array element
     in the selection. No NAs. The values in each column must be >= 1 and
     <= the extent of the array-like object along the corresponding dimension.
     When using an Mindex to subset an array-like object, the returned
     value is a vector-like object (i.e. no dimensions) of length the
     number of rows in the Mindex.

       Example:
         a <- array(101:124, 4:2)
         Mindex <- rbind(c(3, 2, 1),
                         c(2, 1, 1),
                         c(4, 3, 2),
                         c(2, 1, 1))
         a[Mindex]

  3. N-dimensional index (also called "Nindex"): A list with one list element
     per dimension in the array-like object. Each list element must be a
     subscript describing the selection along the corresponding dimension
     of the array-like object.
     IMPORTANT: A NULL subscript is interpreted as a missing subscript,
     that is, as a subscript that runs along the full extend of the
     corresponding dimension of the array-like object. This means that
     before an Nindex can be used in a call to `[`, `[<-`, `[[` or `[[<-`,
     the NULL list elements in it must be replaced with objects of
     class "name".
     When using an Nindex to subset an array-like object, the returned value
     is an array-like object of dimensions the lengths of the selections
     along each dimensions.

       Examples:
         a <- array(101:124, 4:2)

         ## Normalized Nindex:

         Nindex <- list(c(1, 4, 1), NULL, 1)
         ## Same as a[c(1, 4, 1), , 1, drop=FALSE]:
         DelayedArray:::subset_by_Nindex(a, Nindex)

         Nindex <- list(integer(0), NULL, 1)
         ## Same as a[integer(0), , 1, drop=FALSE]:
         DelayedArray:::subset_by_Nindex(a, Nindex)

         ## Non-normalized Nindex:

         Nindex <- list(-3, NULL, 1)
         Nindex <- DelayedArray:::normalizeNindex(Nindex, a)
         ## Same as a[-3, , 1, drop=FALSE]:
         DelayedArray:::subset_by_Nindex(a, Nindex)

         Nindex <- list(IRanges(2, 4), NULL, 1)
         Nindex <- DelayedArray:::normalizeNindex(Nindex, a)
         ## Same as a[2:4, , 1, drop=FALSE]:
         DelayedArray:::subset_by_Nindex(a, Nindex)

         dimnames(a)[[1]] <- LETTERS[1:4]
         Nindex <- list(c("D", "B"), NULL, 1)
         Nindex <- DelayedArray:::normalizeNindex(Nindex, a)
         ## Same as a[c("D", "B"), , 1, drop=FALSE]:
         DelayedArray:::subset_by_Nindex(a, Nindex)
*/

#define ERRMSG_BUF_LENGTH 256

static char *errmsg_buf()
{
	static char buf[ERRMSG_BUF_LENGTH];

	return buf;
}

#define PRINT_TO_ERRMSG_BUF(...) \
	snprintf(errmsg_buf(), ERRMSG_BUF_LENGTH, __VA_ARGS__)

#define NOT_A_FINITE_NUMBER(x) \
	(R_IsNA(x) || R_IsNaN(x) || (x) == R_PosInf || (x) == R_NegInf)

static inline int get_untrusted_elt(SEXP x, int i, long long int *val,
				    const char *what)
{
	int tmp1;
	double tmp2;

	if (IS_INTEGER(x)) {
		tmp1 = INTEGER(x)[i];
		if (tmp1 == NA_INTEGER) {
			PRINT_TO_ERRMSG_BUF("%s[%d] is NA", what, i + 1);
			return -1;
		}
		*val = (long long int) tmp1;
	} else {
		tmp2 = REAL(x)[i];
		if (NOT_A_FINITE_NUMBER(tmp2)) {
			PRINT_TO_ERRMSG_BUF("%s[%d] is NA or NaN "
					    "or not a finite number",
					    what, i + 1);
			return -1;
		}
		if (tmp2 > (double) LLONG_MAX || tmp2 < (double) LLONG_MIN) {
			PRINT_TO_ERRMSG_BUF("%s[%d] is too large (= %e)",
					    what, i + 1, tmp2);
			return -1;
		}
		*val = (long long int) tmp2;
	}
	return 0;
}

static int get_matrix_nrow_ncol(SEXP m, int *nrow, int *ncol)
{
	SEXP m_dim;
	R_xlen_t len;

	if (!IS_INTEGER(m))
		return -1;
	m_dim = GET_DIM(m);
	if (m_dim == R_NilValue) {
		*nrow = 1;
		len = XLENGTH(m);
		if (len > INT_MAX)
			error("too many dimensions");
		*ncol = (int) len;
		return 0;
	}
	if (LENGTH(m_dim) == 2) {
		*nrow = INTEGER(m_dim)[0];
		*ncol = INTEGER(m_dim)[1];
		return 0;
	}
	return -1;
}

static long long int safe_dim_prod(const int *dim, int ndim)
{
	long long int prod;
	int i, d;

	prod = 1;
	reset_ovflow_flag();
	for (i = 0; i < ndim; i++) {
		d = dim[i];
		if (d == NA_INTEGER || d < 0)
			error("'dim' cannot contain NAs or negative values");
		prod = safe_llint_mult(prod, (long long int) d);
	}
	if (get_ovflow_flag())
		error("dimensions are too big");
	return prod;
}

static void set_names(SEXP x, SEXP names)
{
	SEXP x_names;

	if (names == R_NilValue)
		return;
	x_names = PROTECT(duplicate(names));
	SET_NAMES(x, x_names);
	UNPROTECT(1);
	return;
}

static void set_rownames(SEXP m, SEXP rownames)
{
	SEXP m_dimnames, m_rownames;

	if (rownames == R_NilValue)
		return;
	m_dimnames = PROTECT(NEW_LIST(2));
	SET_DIMNAMES(m, m_dimnames);
	UNPROTECT(1);

	m_rownames = PROTECT(duplicate(rownames));
	SET_VECTOR_ELT(m_dimnames, 0, m_rownames);
	UNPROTECT(1);
	return;
}


/****************************************************************************
 * Convert back and forth between Lindex and Mindex
 */

static int L2M(const int *dim, int ndim, int dim_nrow, SEXP L, int *M)
{
	int M_nrow, i, ret, along, d;
	long long int x;
	R_xlen_t dim_off, M_off;

	M_nrow = LENGTH(L);
	for (i = 0; i < M_nrow; i++) {
		ret = get_untrusted_elt(L, i, &x, "Lindex");
		if (ret < 0)
			return -1;
		if (x < 1) {
			PRINT_TO_ERRMSG_BUF("Lindex[%d] is < 1", i + 1);
			return -1;
		}
		x--;
		dim_off = 0;
		if (dim_nrow != 1)
			dim_off += i;
		M_off = i;
		for (along = 0; along < ndim; along++) {
			d = dim[dim_off];
			if (d == NA_INTEGER || d < 0) {
				PRINT_TO_ERRMSG_BUF("'dim' cannot contain NAs "
						    "or negative values");
				return -1;
			}
			if (d == 0) {
				PRINT_TO_ERRMSG_BUF("'dim' cannot contain "
						    "zeroes (unless 'Lindex' "
						    "is empty)");
				return -1;
			}
			M[M_off] = x % d + 1;
			x /= d;
			dim_off += dim_nrow;
			M_off += M_nrow;
		}
		if (x != 0) {
			PRINT_TO_ERRMSG_BUF("Lindex[%d] is > prod(dim)", i + 1);
			return -1;
		}
	}
	return 0;
}

static int M2L(const int *dim, int ndim, int dim_nrow,
	       const int *M, int M_nrow, SEXP L)
{
	R_xlen_t dim_len, M_len, dim_off, M_off;
	int i, along, m, d;
	long long int x;
	double val;

	dim_len = dim_nrow * ndim;
	M_len = M_nrow * ndim;
	if (TYPEOF(L) != INTSXP)
		reset_ovflow_flag();
	for (i = 0; i < M_nrow; i++) {
		x = 0;
		dim_off = dim_len;
		if (dim_nrow != 1)
			dim_off += i;
		M_off = M_len + i;
		for (along = ndim - 1; along >= 0; along--) {
			dim_off -= dim_nrow;
			M_off -= M_nrow;
			d = dim[dim_off];
			if (d == NA_INTEGER || d < 0) {
				PRINT_TO_ERRMSG_BUF("'dim' cannot contain NAs "
						    "or negative values");
				return -1;
			}
			if (d == 0) {
				PRINT_TO_ERRMSG_BUF("'dim' cannot contain "
						    "zeroes (unless 'Mindex' "
						    "is empty)");
				return -1;
			}
			m = M[M_off];
			if (m == NA_INTEGER || m < 1 || m > d) {
				PRINT_TO_ERRMSG_BUF("Mindex[%d, %d] is NA "
						    "or < 1 or > dim[%d]",
						    i + 1, along + 1,
						    along + 1);
				return -1;
			}
			if (TYPEOF(L) == INTSXP) {
				x *= d;
				x += m - 1;
			} else {
				x = safe_llint_mult(x, (long long int) d);
				x = safe_llint_add(x, (long long int) m - 1);
			}
		}
		if (TYPEOF(L) == INTSXP) {
			x++;
			INTEGER(L)[i] = (int) x;
		} else {
			x = safe_llint_add(x, 1);
			val = (double) x;
			if (get_ovflow_flag() || (long long int) val != x) {
				PRINT_TO_ERRMSG_BUF("dimensions in dim[%d, ] "
						    "are too big", i + 1);
				return -1;
			}
			REAL(L)[i] = val;
		}
	}
	return 0;
}

/* --- .Call ENTRY POINT --- */
SEXP C_Lindex2Mindex(SEXP Lindex, SEXP dim, SEXP use_names)
{
	int ret, dim_nrow, dim_ncol;
	R_xlen_t Lindex_len;
	SEXP ans;

	/* Check 'dim'. */
	ret = get_matrix_nrow_ncol(dim, &dim_nrow, &dim_ncol);
	if (ret < 0)
		error("'dim' must be an integer vector (or matrix)");

	/* Check 'Lindex'. */
	if (!(IS_INTEGER(Lindex) || IS_NUMERIC(Lindex)))
		error("'Lindex' must be an integer (or numeric) vector");
	Lindex_len = XLENGTH(Lindex);
	/* R does not support matrices with dimensions > INT_MAX yet (as
	   of April 2020) which means that we cannot support a long linear
	   index. */
	if (Lindex_len > INT_MAX)
		error("'Lindex' is too long");
	if (dim_nrow != 1 && dim_nrow != Lindex_len)
		error("'dim' must have a single row or "
		      "one row per element in 'Lindex'");

	ans = PROTECT(allocMatrix(INTSXP, (int) Lindex_len, dim_ncol));

	ret = L2M(INTEGER(dim), dim_ncol, dim_nrow, Lindex, INTEGER(ans));
	if (ret < 0) {
		UNPROTECT(1);
		error(errmsg_buf());
	}

	if (LOGICAL(use_names)[0])
		set_rownames(ans, GET_NAMES(Lindex));

	UNPROTECT(1);
	return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP C_Mindex2Lindex(SEXP Mindex, SEXP dim, SEXP use_names, SEXP as_integer)
{
	int ret, dim_nrow, dim_ncol, Mindex_nrow, Mindex_ncol;
	long long int dim_prod;
	SEXPTYPE ans_Rtype;
	SEXP ans;

	/* Check 'dim'. */
	ret = get_matrix_nrow_ncol(dim, &dim_nrow, &dim_ncol);
	if (ret < 0)
		error("'dim' must be an integer vector (or matrix)");

	/* Check 'Mindex'. */
	ret = get_matrix_nrow_ncol(Mindex, &Mindex_nrow, &Mindex_ncol);
	if (ret < 0)
		error("'Mindex' must be an integer matrix (or vector)");
	if (Mindex_ncol != dim_ncol)
		error("'Mindex' must have one %s per dimension",
		      GET_DIM(Mindex) != R_NilValue ? "column" : "element");
	if (dim_nrow != 1 && dim_nrow != Mindex_nrow)
		error("'dim' must have a single row or "
		      "the same number of rows as 'Mindex'");

	/* Determine the type of the linear index (integer or numeric). */
	if (LOGICAL(as_integer)[0]) {
		/* Risky! We will return garbage if some L-index values
		   don't fit in the integer type. */
		ans_Rtype = INTSXP;
	} else {
		if (dim_nrow == 1) {
			/* All the L-index values we're going to compute are
			   guaranteed to be <= prod(dim). */
			dim_prod = safe_dim_prod(INTEGER(dim), dim_ncol);
			ans_Rtype = dim_prod <= INT_MAX ? INTSXP : REALSXP;
		} else {
			/* Computing an upper bound for the L-index values
			   would require walking on the entire 'dim' matrix
			   so would be too expensive. */
			ans_Rtype = REALSXP;
		}
	}
	ans = PROTECT(allocVector(ans_Rtype, Mindex_nrow));

	ret = M2L(INTEGER(dim), dim_ncol, dim_nrow,
		  INTEGER(Mindex), Mindex_nrow, ans);
	if (ret < 0) {
		UNPROTECT(1);
		error(errmsg_buf());
	}

	if (LOGICAL(use_names)[0])
		set_names(ans, GET_ROWNAMES(Mindex));

	UNPROTECT(1);
	return ans;
}

