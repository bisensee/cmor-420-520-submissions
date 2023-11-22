#include "abstract_matrix.hpp"
#include<iostream>

#ifndef _SPARSEMATRIXCSR
#define _SPARSEMATRIXCSR

class SparseMatrixCSR: public AbstractMatrix{
public:
  SparseMatrixCSR(int m, int n, int * row_indices, int * col_indices, double * values, int num_vals): AbstractMatrix(m, n){
    _nzval = values;
    _col_index = col_indices;
    _row_index = new int [m + 1];
    _row_index[0] = 0;
    int ctr = 1;
    int comp = row_indices[0];
    int idx = 1;
    for(int r = 1; r < num_vals; ++r){
      if(row_indices[r] == comp){
	ctr++;
      }
      else{
	_row_index[idx] = _row_index[idx - 1] + ctr;
	idx++;
	ctr = 1;
	comp = row_indices[r];
      }
    }
    _row_index[m] = _row_index[m - 1] + ctr;
  }
  ~SparseMatrixCSR(){
    delete[] _nzval;
    delete[] _row_index;
    delete[] _col_index;
  }
  double * get_nzval() const{ return _nzval; }

  int * get_row_index() const{ return _row_index; }

  int * get_col_index() const{ return _col_index; }

  virtual double operator()(int row, int column) const;

private:
  double * _nzval;
  int * _col_index;
  int * _row_index;
};

#endif

Vector operator*(SparseMatrixCSR & A, Vector & x);
