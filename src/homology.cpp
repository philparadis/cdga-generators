#include "homology.h"

#include <assert.h>
#include <NTL/LLL.h>

NTL_CLIENT

void FindCocycleBasis(int **differential_matrix, int rows_size, int cols_size, OrderedLCBasis &cocycleBasis, const OrderedBasis &source)
{
	mat_ZZ D, U;
	long rank;
	ZZ d;
	
	// Empty the parameter "cocycleBasis" that was passed
	cocycleBasis.clear();

	// If the differential is zero, then a basis of cocycles is just the basis of the space that was passed to us
	if (differential_matrix == 0 || rows_size == 0 || cols_size == 0) {
		OrderedBasis::const_iterator iter;
		for (iter = source.begin(); iter != source.end(); iter++) {
			LinearCombination lc;
			lc.AddTerm(1, *iter);
			cocycleBasis.push_back(lc);
		}
		return;
	}
	
	// Note: We need D to be the transpose of "differential_matrix", because of the "reverse" convention used by NTL.
	D.SetDims(cols_size, rows_size);
	for (int i=0; i<rows_size; i++) {
		for (int j=0; j<cols_size; j++) {
			D[j][i] = differential_matrix[i][j];
		}
	}
	
	cerr << D << endl;

	// Use the LLL algorithm to compute the rank of the matrix D.
	// The first cols_size - rank rows of U will be a basis of the kernel of D.
	rank = LLL(d, D, U);

	int dim_ker = cols_size - rank;

	// Loop through the first dim_ker rows of U and store the result into cocycleBasis
	for (int i=0; i<dim_ker; i++) {
		LinearCombination lc;
		for (int j=0; j<cols_size; j++) {
			assert(NumBits(U[i][j]) < 32);
			int coeff = to_int(U[i][j]);
			if (coeff != 0) {
				lc.AddTerm(coeff, source[j]);
			}
		}
		cocycleBasis.push_back(lc);
	}
}
void GetCoordinates(const LinearCombination &lc, vec_ZZ &coordinates, const OrderedBasis &basis)
{
	if (basis.empty()) {
		return;
	}
	int size = basis.size();
	int *coords_array = new int[size];
	lc.GetCoordinates(coords_array, basis);
	for (int i=0; i<size; i++) {
		coordinates[i] = coords_array[i];
	}
	delete [] coords_array;
}

void FindHomologyBasis(int **d1, int d1_rows_size, int d1_cols_size, int **d2, int d2_rows_size, int d2_cols_size, OrderedLCBasis &homologyBasis, OrderedLCBasis &imageBasis, const OrderedBasis &source)
{
	OrderedLCBasis cocyclesBasis;

	// First, find a basis of cocycles
	FindCocycleBasis(d2, d2_rows_size, d2_cols_size, cocyclesBasis, source);

	// If the differential d1 is zero, we are done here, we can return homologyBasis = cocyclesBasis
	if (d1 == 0 || d1_rows_size == 0 || d1_cols_size == 0) {
		homologyBasis = cocyclesBasis;
		return;
	}

	mat_ZZ D1;
	long rank;
	ZZ d;

	// Note: We need D1 to be the transpose of "d1", because of the "reverse" convention used by NTL.
	D1.SetDims(d1_cols_size, d1_rows_size);
	for (int i=0; i<d1_rows_size; i++) {
		for (int j=0; j<d1_cols_size; j++) {
			D1[j][i] = d1[i][j];
		}
	}

	cerr << D1 << endl;

	rank = image(d, D1);

	// The dimension of ker(d_n) / im(d_{n-1}) is supposed to be dim_ker - dim_img.
	int dim_ker = (int)cocyclesBasis.size();
	int dim_img = rank;

	// Loop through the last dim_img rows of U and store the result into imageBasis
	imageBasis.clear();
	for (int i=d1_cols_size-dim_img; i<d1_cols_size; i++) {
		LinearCombination lc;
		for (int j=0; j<d1_rows_size; j++) {
			assert(NumBits(D1[i][j]) < 32);
			int coeff = to_int(D1[i][j]);
			if (coeff != 0) {
				lc.AddTerm(coeff, source[j]);
			}
		}
		imageBasis.push_back(lc);
	}

	homologyBasis.clear();
	// Here, we need to extend the basis of the image of im(d_{n-1}) to a basis of ker(d_n)
	// We will do this by appending as many elements of cocycleBasis as possible, making sure our set of vectors
	// is linearly independant at every step.
	// Note: This is very unefficient and I wish I could do it in one step, by reducing a single matrix to it's HNF, but
	// I can't figure out how to do this with NTL. In any case I think NTL is fast enough so that it shouldn't be the bottleneck here.

	// TMP:

	mat_ZZ extended_basis;
	extended_basis.SetDims(dim_ker - dim_img, d1_rows_size); // Each column vector of this matrix will be a basis vector in the extended basis
	int extended_basis_vectors_found = 0;
	int cocycle_index = 0;
	if (dim_img < dim_ker) {
		do {
			mat_ZZ D;
			D.SetDims(dim_img+extended_basis_vectors_found+1, d1_rows_size);
			for (int i=0; i<dim_img; i++) {
				D[i] = D1[d1_cols_size-dim_img+i];
			}
			for (int i=dim_img; i<dim_img+extended_basis_vectors_found; i++) {
				D[i] = extended_basis[i-dim_img];
			}
			// Check for linear independence with the following new vector added
			GetCoordinates(cocyclesBasis[cocycle_index], D[dim_img+extended_basis_vectors_found], source);
			int rank = LLL(d, D);
			if (rank > dim_img+extended_basis_vectors_found) {
				// The rank has increased, so the extended set of vectors is still linearly independant
				GetCoordinates(cocyclesBasis[cocycle_index], extended_basis[extended_basis_vectors_found], source);
				homologyBasis.push_back(cocyclesBasis[cocycle_index]);
				++extended_basis_vectors_found;
			}
			cocycle_index++;
		} while (dim_img+extended_basis_vectors_found < dim_ker && cocycle_index < dim_ker);

		if (dim_img+extended_basis_vectors_found < dim_ker)
			throw logic_error("Fatal error. Failed to compute a quotient space basis.");
	}
}