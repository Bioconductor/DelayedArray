/****************************************************************************
 *                           C_dgCMatrix_rowsum()                           *
 *                            Author: H. Pag\`es                            *
 ****************************************************************************/
#include "DelayedArray.h"
#include "S4Vectors_interface.h"


static void check_group(SEXP group, int x_nrow, int ngroup)
{
	int i, g;

	if (!IS_INTEGER(group))
		error("the grouping vector must be "
		      "an integer vector or factor");
	if (LENGTH(group) != x_nrow)
		error("the grouping vector must have "
		      "one element per row in 'x'");
	for (i = 0; i < x_nrow; i++) {
		g = INTEGER(group)[i];
		if (g == NA_INTEGER) {
			if (ngroup < 1)
				error("'ngroup' must be >= 1 when 'group' "
				      "contains missing values");
		} else {
			if (g < 1 || g > ngroup)
				error("all non-NA values in 'group' must "
				      "be >= 1 and <= 'ngroup'");
		}
	}
	return;
}

static void compute_rowsum(const double *x, const int *row, int x_len,
			   const int *group, double *out, int out_len,
			   int narm)
{
	int i, g;
	double x_elt;

	for (i = 0; i < out_len; i++)
		out[i] = 0.0;
	for (i = 0; i < x_len; i++) {
		g = group[row[i]];
		if (g == NA_INTEGER)
			g = out_len;
		g--;  // from 1-base to 0-base
		x_elt = x[i];
		if (narm && (R_IsNA(x_elt) || R_IsNaN(x_elt)))
			continue;
		out[g] += x_elt;
	}
	return;
}

/* --- .Call ENTRY POINT --- */
SEXP C_dgCMatrix_rowsum(SEXP x, SEXP group, SEXP ngroup, SEXP na_rm)
{
	SEXP x_Dim, x_x, x_p, x_i, ans, ans_dim;
	int x_nrow, x_ncol, narm, ans_nrow, ans_len, j, offset, count;
	double *ans_p;

	x_Dim = GET_SLOT(x, install("Dim"));
	x_nrow = INTEGER(x_Dim)[0];
	x_ncol = INTEGER(x_Dim)[1];
	x_x = GET_SLOT(x, install("x"));
	x_p = GET_SLOT(x, install("p"));
	x_i = GET_SLOT(x, install("i"));
	narm = LOGICAL(na_rm)[0];

	ans_nrow = INTEGER(ngroup)[0];
	check_group(group, x_nrow, ans_nrow);

	reset_ovflow_flag();
	ans_len = safe_int_mult(ans_nrow, x_ncol);
	if (get_ovflow_flag())
		error("too many groups (matrix of sums will be too big)");

	ans = PROTECT(allocVector(REALSXP, ans_len));
	ans_p = REAL(ans);
	for (j = 0; j < x_ncol; j++) {
		offset = INTEGER(x_p)[j];
		count = INTEGER(x_p)[j + 1] - offset;
		compute_rowsum(REAL(x_x) + offset, INTEGER(x_i) + offset, count,
			       INTEGER(group), ans_p, ans_nrow, narm);
		ans_p += ans_nrow;
	}

	ans_dim = PROTECT(NEW_INTEGER(2));
	INTEGER(ans_dim)[0] = ans_nrow;
	INTEGER(ans_dim)[1] = x_ncol;
	SET_DIM(ans, ans_dim);

	UNPROTECT(2);
	return ans;
}

