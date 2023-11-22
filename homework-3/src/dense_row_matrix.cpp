#include "dense_row_matrix.hpp"

double & DenseRowMatrix::operator()(int row, int column){
  return get_data()[row*num_columns() + column];
}

Vector operator*(DenseRowMatrix & A, Vector & x){
  Vector fin = Vector(A.num_rows());
  
  for(int i = 0; i < A.num_rows(); ++i){
    double sum = 0;
    for(int j = 0; j < A.num_columns(); ++j){
      sum += A(i, j)*x(j);
    }
    fin(i) = sum;
  }
  return fin;
}
