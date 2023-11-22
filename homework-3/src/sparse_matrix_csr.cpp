#include "sparse_matrix_csr.hpp"

double SparseMatrixCSR::operator()(int row, int column) const{
  int start = get_row_index()[row];
  int end = get_row_index()[row+1];
  for(int j = start; j < end; ++j){
    if(_col_index[j] == column){
      return _nzval[j];
    }
  }
  return 0;
}

Vector operator*(SparseMatrixCSR & A, Vector & x){
  Vector fin = Vector(A.num_rows());
  
  for (int i = 0; i < A.num_rows(); ++i){
    int row_start = A.get_row_index()[i];
    int row_end = A.get_row_index()[i+1];
    double sum = 0;
    for(int j = row_start; j < row_end; ++j){
      int col = A.get_col_index()[j];
      sum += A(i, col) * x(col);
    }
    fin(i) = sum;
  }
  return fin;
}
