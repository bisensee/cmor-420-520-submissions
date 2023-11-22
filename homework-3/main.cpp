#include <iostream>
#include "vector.hpp"
#include "dense_row_matrix.hpp"
#include "sparse_matrix_csr.hpp"
#include <math.h>
#include <cmath>
#include <chrono>

// To run this file, run "make" and then "./main"
// in the terminal.

// Run "make clean" to delete old build files. 

// both verification and timing should be implemented
// in this main.cpp file. Currently, it is intended only
// to demonstrate the usage of the "Vector" class.
double norm(const Vector a){
  double sum = 0;
  for(int i = 0; i < a.length(); i++){
    sum += a(i)*a(i);
  }
  sum = sqrt(sum);
  return sum;
}

int main(void) {

  int n = 100;
  double h = 1.0/n;
  double a = 0.5*h*h;
  DenseRowMatrix A_dense = DenseRowMatrix(n+1, n+1);
  A_dense(0,0) = 1*n*n;
  for(int i = 0; i < n; i++){
    A_dense(i, i+1) = -1*n*n;
    A_dense(i+1, i) = -1*n*n;
    A_dense(i+1, i+1) = 2*n*n;
  }
  A_dense(n,n) = 1*n*n;

  int num_vals = (n-1)*3 + 4;
  double * vals = new double [num_vals];
  int * row_indices = new int [num_vals];
  int * col_indices = new int [num_vals];

  row_indices[0] = 0;
  row_indices[1] = 0;
  int row_ctr = 1;
  for(int i = 2; i < num_vals; i+=3){
    row_indices[i] = row_ctr;
    row_indices[i+1] = row_ctr;
    row_indices[i+2] = row_ctr;
    row_ctr++;
  }
  row_indices[(n-1)*3 + 2] = n;
  row_indices[(n-1)*3 + 3] = n;
  
  col_indices[0] = 0;
  col_indices[1] = 1;
  int col_ctr = 0;
  for(int i = 2; i < num_vals; i+=3){
    col_indices[i] = col_ctr;
    col_indices[i+1] = col_ctr+1;
    col_indices[i+2] = col_ctr+2;
    col_ctr++;
  }

  col_indices[(n-1)*3 + 2] = n-1;
  col_indices[(n-1)*3 + 3] = n;

  vals[0] = 1*n*n;
  vals[1] = -1*n*n;
  for(int i = 2; i < (n*3-2); i+=3){
    vals[i] = -1*n*n;
    vals[i+1] = 2*n*n;
    vals[i+2] = -1*n*n;
  }
  vals[(n-1)*3 + 2] = -1*n*n;
  vals[(n-1)*3 + 3] = 1*n*n;

  SparseMatrixCSR A_sparse = SparseMatrixCSR(n+1, n+1, row_indices, col_indices, vals, num_vals);

  Vector b = Vector(n+1);

  for(int i = 0; i < n+1; i++){
    b(i) = cos(M_PI * i * h);
  }

  Vector r = b;
  double norm_r = norm(r);
  int iteration = 0;
  Vector u = Vector(n+1);
  for(int i = 0; i < n+1; i++){
    u(i) = 0;
  }

  std::cout << "********Allocating version********" << std::endl;
  std::cout << "DENSE MATRIX" << std::endl;
  std::cout << "Iteration " << iteration << ", norm(r) = " << norm_r << std::endl;

  auto start = std::chrono::system_clock::now();
  while(norm_r > 0.001){
    iteration += 1;
    u += a*r;
    r = b - A_dense*u;
    norm_r = norm(r);
  }
  auto stop  = std::chrono::system_clock::now();
  double elapsed = double(std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count());
  std::cout << "Iteration " << iteration << ", norm(r) = " << norm_r << std::endl;
  std::cout << "Elapsed time: " << elapsed/1000.0 << " seconds" << std::endl;
  std::cout << "Average time per iteration: " << elapsed/(1000.0*iteration) << " seconds" << std::endl;
  std::cout << std::endl;
  
  r = b;
  norm_r = norm(r);
  iteration = 0;
  for(int i = 0; i < n+1; i++){
    u(i) = 0;
  }
  std::cout << "SPARSE MATRIX" << std::endl;
  std::cout << "Iteration " << iteration << ", norm(r) = " << norm_r << std::endl;

  start = std::chrono::system_clock::now();
  while(norm_r > 0.001){
    iteration += 1;
    u += a*r;
    r = b - A_sparse*u;
    norm_r = norm(r);
  }
  stop = std::chrono::system_clock::now();
  elapsed = double(std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count());
  std::cout << "Iteration " << iteration << " , norm(r) = " << norm_r << std::endl;
  std::cout << "Elapsed time: " << elapsed/1000.0 << " seconds" << std::endl;
  std::cout << "Average time per iteration: " << elapsed/(1000.0*iteration) << " seconds" << std::endl;

  
  std::cout << std::endl;
  std::cout << std::endl;


  std::cout << "******Non-allocating version******" << std::endl;
  r = b;
  norm_r = norm(r);
  iteration = 0;
  for(int i = 0; i < n+1; i++){
    u(i) = 0;
  }
  Vector firstOp = r;
  Vector secondOp = u;
  Vector thirdOp = b;
  Vector product = Vector(A_dense.num_rows());
  double sum = 0;
  std::cout << "DENSE MATRIX" << std::endl;
  std::cout << "Iteration " << iteration << ", norm(r) = " << norm_r << std::endl;
  start = std::chrono::system_clock::now();
  while(norm_r > 0.001){
    iteration += 1;
    firstOp *= a;
    u += firstOp;
    secondOp = u;
    for(int i = 0; i < A_dense.num_rows(); i++){
      sum = 0;
      for(int j = 0; j < A_dense.num_columns(); j++){
	sum += A_dense(i,j)*secondOp(j);
      }
      product(i) = sum;
    }
    thirdOp -= product;
    r = thirdOp;
    
    norm_r = norm(r);
    firstOp = r;
    thirdOp = b;
  }
  stop  = std::chrono::system_clock::now();
  elapsed = double(std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count());
  std::cout << "Iteration " << iteration << ", norm(r) = " << norm_r << std::endl;
  std::cout << "Elapsed time: " << elapsed/1000.0 << " seconds" << std::endl;
  std::cout << "Average time per iteration: " << elapsed/(1000.0*iteration) << " seconds" << std::endl;
  std::cout << std::endl;


  
  r = b;
  norm_r = norm(r);
  iteration = 0;
  for(int i = 0; i < n+1; i++){
    u(i) = 0;
  }
  firstOp = r;
  thirdOp = b;
  sum = 0;
  std::cout << "SPARSE MATRIX" << std::endl;
  std::cout << "Iteration " << iteration << ", norm(r) = " << norm_r << std::endl;

  start = std::chrono::system_clock::now();
  while(norm_r > 0.001){
    iteration += 1;
    firstOp *= a;
    u += firstOp;
    secondOp = u;

    for (int i = 0; i < A_sparse.num_rows(); ++i){
      int row_start = A_sparse.get_row_index()[i];
      int row_end = A_sparse.get_row_index()[i+1];
      sum = 0;
      for(int j = row_start; j < row_end; ++j){
	int col = A_sparse.get_col_index()[j];
	sum += A_sparse(i, col) * secondOp(col);
      }
      product(i) = sum;
    }
    
    thirdOp -= product;
    r = thirdOp;
    
    norm_r = norm(r);
    firstOp = r;
    thirdOp = b;
  }
  stop = std::chrono::system_clock::now();
  elapsed = double(std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count());
  std::cout << "Iteration " << iteration << ", norm(r) = " << norm_r << std::endl;
  std::cout << "Elapsed time: " << elapsed/1000.0 << " seconds" << std::endl;
  std::cout << "Average time per iteration: " << elapsed/(1000.0*iteration) << " seconds" << std::endl;
}
