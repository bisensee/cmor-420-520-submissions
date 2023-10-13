#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"

Matrix init_matrix_contiguous(int M, int N){
  double ** mat = (double **) malloc(sizeof(double *) * M);
  mat[0] = (double *) malloc(sizeof(double) * M * N);
  for (int i = 1; i < M; ++i){
    mat[i] = mat[i-1] + N;
  }
  
  Matrix ma = {.m = M, .n = N, .data = mat};
  return ma;
}
  
Matrix init_matrix_non_contiguous(int M, int N){
  double ** mat = (double **) malloc(sizeof(double *) * M);
  for (int i = 0; i < M; ++i){
    mat[i] = (double *) malloc(sizeof(double) * N);
  }
  
  Matrix ma = {.m = M, .n = N, mat};
  return ma;
}

void free_matrix_contiguous(Matrix * mat){
  free(mat->data[0]);
  free(mat->data);
}

void free_matrix_non_contiguous(Matrix * mat){
  for(int i = 0; i < mat->m; ++i){
    free(mat->data[i]);
  }
  free(mat->data);
}
void matvec(double * b, const Matrix * A, const double * x){
  double * tmp = b;

  for(int M = 0; M < A->m; ++M){
    double sum = 0;
    for(int N = 0; N < A->n; ++N){
      sum += A->data[M][N]*x[N];
    }
    tmp[M] = sum;
  }
  *b = *tmp;
}

void matvec_transpose(double * b, const Matrix * A, const double * x){
  double * tmp = b;

  for(int M = 0; M < A->m; ++M){
    double sum = 0;
    for(int N = 0; N < A->n; ++N){
      sum += A->data[N][M]*x[N];
    }
    tmp[M] = sum;
  }
  *b = *tmp;
}
