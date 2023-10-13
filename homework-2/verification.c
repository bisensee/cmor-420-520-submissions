#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"

int main(){


  ////////////////////////////////////
  //contiguous matrix multiplication//
  ////////////////////////////////////
  
  Matrix cont_mat = init_matrix_contiguous(3,3);

  //set matrix mat to specific values
  double arr[3][3] = {{1,2,3}, {4,5,6}, {7,8,9}};
  for(int i = 0; i < 3; ++i){
    for(int j = 0; j < 3; ++j){
      cont_mat.data[i][j] = arr[i][j];
    }
  }

  double * x = (double *) malloc(sizeof(double) * 3);
  double * b = (double *) malloc(sizeof(double) * 3);
  double bExact[3] = {28, 64, 100};
  double bExactTrans[3] = {60, 72, 84};

  //set x, b vectors
  for(int i = 0; i < 3; ++i){
    x[i] = 2*(i+1);
    b[i] = 0;
  }
  matvec(b, &cont_mat, x);

  //contiguous matrix multiplication error
  double err = 0;
  for(int i = 0; i < 3; i++){
    err += fabs(b[i] - bExact[i]);
  }
  printf("Contiguous matrix multiplication error: %f\n", err);

  //reset b vector
  for(int i = 0; i < 3; ++i){
    b[i] = 0;
  }

  matvec_transpose(b, &cont_mat, x);
  //contiguous transpose matrix multiplication error
  err = 0;
  for(int i = 0; i < 3; i++){
    err += fabs(b[i] - bExactTrans[i]);
  }
  printf("Contiguous transpose matrix multiplication error: %f\n", err);
  free_matrix_contiguous(&cont_mat);

  ////////////////////////////////////////
  //non-contiguous matrix multiplication//
  ////////////////////////////////////////

  Matrix non_cont_mat = init_matrix_non_contiguous(3,3);

  //set matrix mat to specific values
  for(int i = 0; i < 3; ++i){
    for(int j = 0; j < 3; ++j){
      non_cont_mat.data[i][j] = arr[i][j];
    }
  }
  
  //set x, b vectors
  for(int i = 0; i < 3; ++i){
    x[i] = 2*(i+1);
    b[i] = 0;
  }
  
  matvec(b, &non_cont_mat, x);

  //contiguous matrix multiplication error
  err = 0;
  for(int i = 0; i < 3; i++){
    err += fabs(b[i] - bExact[i]);
  }
  printf("Non-contiguous matrix multiplication error: %f\n", err);

  //reset b vector. Before transpose matrix multiplication
  for(int i = 0; i < 3; ++i){
    b[i] = 0;
  }

  matvec_transpose(b, &non_cont_mat, x);

  //non-contiguous transpose matrix multiplication error
  err = 0;
  for(int i = 0; i < 3; i++){
    err += fabs(b[i] - bExactTrans[i]);
  }
  printf("Non-contiguous transpose matrix multiplication error: %f\n", err);
  free_matrix_non_contiguous(&non_cont_mat);
}
