#ifndef _HOMOLOGY__H
#define _HOMOLOGY__H

#include "cdga.h"

// The parameter "differential_matrix" must be the matrix of a differential d_n : X^{n} ---> X^{n+1},
// with respect to the ordered basis "source". The method will return a basis of ker(d_n).
// The parameter "source" must be an ordered matrix in dimension n (hence it must have cols_size vectors)
void FindCocycleBasis(int **differential_matrix, int rows_size, int cols_size, OrderedLCBasis &cocycleBasis, const OrderedBasis &source);

// The parameter "d1" must be the matrix of the differential d_{n-1} : X^{n-1} ---> X^n.
// The parameter "d2" must be the matrix of the differential d_n : X^n ---> X^{n+1}.
// The parameter "source" must be an ordered basis in degree n (hence must have d2_cols_size vectors)
// The method will find a basis of cocycles which are not boundaries in X^n.
// The method also returns a basis for im(d_{n-1}) as the parameter "imageBasis"

void FindHomologyBasis(int **d1, int d1_rows_size, int d1_cols_size, int **d2, int d2_rows_size, int d2_cols_size, OrderedLCBasis &homologyBasis, OrderedLCBasis &imageBasis, const OrderedBasis &source);

#endif 