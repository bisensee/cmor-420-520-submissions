#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "matrix.h"

double time_matvec(void (*mv_ptr)(double *, const Matrix *, const double *),
		   double * b, const Matrix * A, const double * x,
		   double num_samples){
  
  clock_t tic, toc;
  tic = clock();
  for(int i = 0; i < num_samples; i++){
    (*mv_ptr)(b, A, x);
  }
  toc = clock();
  double elapsed_avg = (double) (toc - tic);
  elapsed_avg /= CLOCKS_PER_SEC;
  elapsed_avg /= num_samples;
  return elapsed_avg;
}
		   

int main(int argc, char *argv[]){
  int M = 1000;
  int N = 1000;
  int num_s = 10;
  if(argc == 4){
    M = atoi(argv[1]);
    N = atoi(argv[2]);
    num_s = atoi(argv[3]);
  }
  
  Matrix A = init_matrix_contiguous(M, N);
  double * x = (double *) malloc(sizeof(double) * M);
  double * b = (double *) malloc(sizeof(double) * M);

  
  for(int i = 0; i < M; i++){
    x[i] = 1;
    b[i] = 0;
    for(int j = 0; j < M; j++){
      A.data[i][j] = 1;
    }
  }
  
  void (*mv_ptr)(double *, const Matrix *, const double *) = matvec;
  void (*mv_trans_ptr)(double *, const Matrix *, const double *) = matvec_transpose;
  double avg_runtime_mv = time_matvec(mv_ptr, b, &A, x, num_s);
  double avg_runtime_mv_trans = time_matvec(mv_trans_ptr, b, &A, x, num_s);
  printf("Contiguous matrix multiplication, m = %d, n = %d, number of samples = %d, average time = %f\n", M, N, num_s, avg_runtime_mv);
  printf("Contiguous transpose matrix multiplication, m = %d, n = %d, number of samples = %d, average time = %f\n", M, N, num_s, avg_runtime_mv_trans);

  free_matrix_contiguous(&A);

  Matrix B = init_matrix_non_contiguous(M, N);

  for(int i = 0; i  < M; ++i){
    b[i] = 0;
    for(int j = 0; j < M; ++j){
      B.data[i][j] = 0;
    }
  }

  avg_runtime_mv = time_matvec(mv_ptr, b, &B, x, num_s);
  avg_runtime_mv_trans = time_matvec(mv_trans_ptr, b, &B, x, num_s);
  printf("Non-contiguous matrix multiplication, m = %d, n = %d, number of samples = %d, average time = %f\n", M, N, num_s, avg_runtime_mv);
  printf("Non-contiguous transpose matrix multiplication, m = %d, n = %d, number of samples = %d, average time = %f\n\n", M, N, num_s, avg_runtime_mv_trans);

  free_matrix_non_contiguous(&B);
}
