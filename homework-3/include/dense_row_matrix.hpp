#include "abstract_matrix.hpp"
#include<iostream>

#ifndef _DENSEROWMATRIX
#define _DENSEROWMATRIX

class DenseRowMatrix: public AbstractMatrix{
public:
  DenseRowMatrix(int m, int n): AbstractMatrix(m,n){
    _data = new double [m*n];
  }
  ~DenseRowMatrix(){
    delete[] _data;
  }
  double * get_data() const{ return _data; }

  virtual double & operator()(int row, int column);

private:
  double * _data;
};

#endif
