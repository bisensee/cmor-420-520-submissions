#ifndef _MATRIX
#define _MATRIX

typedef struct Matrix{
  int m;
  int n;
  double ** data;
} Matrix;

Matrix init_matrix_contiguous(int M, int N);

Matrix init_matrix_non_contiguous(int M, int N);

void free_matrix_contiguous(Matrix * mat);

void free_matrix_non_contiguous(Matrix * mat);

void matvec(double * b, const Matrix * A, const double * x);

void matvec_transpose(double * b, const Matrix * A, const double * x);

#endif
