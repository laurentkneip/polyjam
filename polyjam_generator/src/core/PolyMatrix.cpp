/*************************************************************************
 *                                                                       *
 * polyjam, a polynomial solver generator for C++                        *
 * Copyright (C) 2015 Laurent Kneip, The Australian National University  *
 *                                                                       *
 * This program is free software: you can redistribute it and/or modify  *
 * it under the terms of the GNU General Public License as published by  *
 * the Free Software Foundation, either version 3 of the License, or     *
 * (at your option) any later version.                                   *
 *                                                                       *
 * This program is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 * GNU General Public License for more details.                          *
 *                                                                       *
 * You should have received a copy of the GNU General Public License     *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>. *
 *                                                                       *
 *************************************************************************/

#include <polyjam/core/PolyMatrix.hpp>
#include <iostream>

using namespace std;

//constructors / destructor

polyjam::core::PolyMatrix::PolyMatrix(
    const Poly & base, size_t rows, size_t cols, bool diag )
{
  if( rows == 0 )
  {
    cout << "Error: cannot construct a poly-matrix with 0 rows" << endl;
    rows = 1;
  }
  if( cols == 0 )
  {
    cout << "Error: cannot construct a poly-matrix with 0 cols" << endl;
    cols = 1;
  }
  
  _matrix.reserve(rows);
  
  for( size_t row = 0; row < rows; row++ )
  {
    std::vector<PolyPtr> newRow;
    newRow.reserve(cols);
    
    for( size_t col = 0; col < cols; col++ )
    {
      if( diag && (row != col) )
        newRow.push_back(PolyPtr( new Poly(base.zero()) ));
      else
        newRow.push_back(PolyPtr( new Poly(base.clone()) ));
    }
    
    _matrix.push_back(newRow);
  }
  
  _degreeLimitation = false;
}

polyjam::core::PolyMatrix::PolyMatrix(
    const Poly & base, size_t rows, size_t cols,
    unsigned int maxDegree, bool diag )
{
  if( rows == 0 )
  {
    cout << "Error: cannot construct a poly-matrix with 0 rows" << endl;
    rows = 1;
  }
  if( cols == 0 )
  {
    cout << "Error: cannot construct a poly-matrix with 0 cols" << endl;
    cols = 1;
  }
  
  _matrix.reserve(rows);
  
  for( size_t row = 0; row < rows; row++ )
  {
    std::vector<PolyPtr> newRow;
    newRow.reserve(cols);
    
    for( size_t col = 0; col < cols; col++ )
    {
      if( diag && (row != col) )
        newRow.push_back(PolyPtr( new Poly(base.zero()) ));
      else
        newRow.push_back(PolyPtr( new Poly(base.clone()) ));
    }
    
    _matrix.push_back(newRow);
  }
  
  _degreeLimitation = true;
  _maxDegree = maxDegree;
}

polyjam::core::PolyMatrix::~PolyMatrix()
{}

// deep stuff

void
polyjam::core::PolyMatrix::copy( const PolyMatrix & copy )
{
  size_t rows = copy.rows();
  size_t cols = copy.cols();
  
  _matrix.clear();
  _matrix.reserve(rows);
  
  for( size_t row = 0; row < rows; row++ )
  {
    std::vector<PolyPtr> newRow;
    newRow.reserve(cols);
    
    for( size_t col = 0; col < cols; col++ )
      newRow.push_back(PolyPtr(new Poly( copy._matrix[row][col]->clone() )));
    
    _matrix.push_back(newRow);
  }
  
  _degreeLimitation = copy._degreeLimitation;
  _maxDegree = copy._maxDegree;
}

polyjam::core::PolyMatrix
polyjam::core::PolyMatrix::clone() const
{
  //create a minimal matrix
  PolyMatrix result( _matrix[0][0]->zero(),1,1,false );
  
  //now use clear its container and fill with elements from this matrix
  size_t _rows = rows();
  size_t _cols = cols();
  
  result._matrix.clear();
  result._matrix.reserve(_rows);
  
  for( size_t row = 0; row < _rows; row++ )
  {
    std::vector<PolyPtr> newRow;
    newRow.reserve(_cols);
    
    for( size_t col = 0; col < _cols; col++ )
      newRow.push_back(PolyPtr(new Poly( _matrix[row][col]->clone() )));
    
    result._matrix.push_back(newRow);
  }
  
  result._degreeLimitation = _degreeLimitation;
  result._maxDegree = _maxDegree;
  
  return result;
}

//output/access

size_t
polyjam::core::PolyMatrix::rows() const
{
  return _matrix.size();
}

size_t
polyjam::core::PolyMatrix::cols() const
{
  if( _matrix.size() == 0 )
    return 0;
  return _matrix[0].size();
}

polyjam::core::Poly &
polyjam::core::PolyMatrix::operator()( size_t row, size_t col )
{
  //check that we attempt to access an existing element (0,0) always does!
  if( row >= _matrix.size() || col >= _matrix[0].size() )
  {
    cout << "Error: attempt to access beyond matrix size" << endl;
    return (*_matrix[0][0]);
  }
  return (*_matrix[row][col]);
}

polyjam::core::Poly &
polyjam::core::PolyMatrix::operator[]( size_t index )
{
  size_t totalNumberElements = _matrix.size() * _matrix[0].size();
  
  if( index >= totalNumberElements )
  {
    cout << "Error: attempt to access beyond matrix size" << endl;
    return (*_matrix[0][0]);
  }
  
  // row major access
  size_t row = index / cols();
  size_t col = index % cols();
  
  return (*_matrix[row][col]);
}

//operations

polyjam::core::PolyMatrix
polyjam::core::PolyMatrix::negation() const
{  
  //create a minimal matrix
  PolyMatrix result( _matrix[0][0]->zero(),1,1,false );
  
  //now use clear its container and fill with elements from this matrix
  size_t _rows = rows();
  size_t _cols = cols();
  
  result._matrix.clear();
  result._matrix.reserve(_rows);
  
  for( size_t row = 0; row < _rows; row++ )
  {
    std::vector<PolyPtr> newRow;
    newRow.reserve(_cols);
    
    for( size_t col = 0; col < _cols; col++ )
      newRow.push_back(PolyPtr(new Poly( _matrix[row][col]->negation() )));
    
    result._matrix.push_back(newRow);
  }
  
  result._degreeLimitation = _degreeLimitation;
  result._maxDegree = _maxDegree;
  
  result.clean();
  return result;
}

polyjam::core::PolyMatrix
polyjam::core::PolyMatrix::transpose() const
{
  //create a minimal matrix
  PolyMatrix result( _matrix[0][0]->zero(),1,1,false );
  
  //now use clear its container and fill with elements from this matrix
  size_t _rows = cols();
  size_t _cols = rows();
  
  result._matrix.clear();
  result._matrix.reserve(_rows);
  
  for( size_t row = 0; row < _rows; row++ )
  {
    std::vector<PolyPtr> newRow;
    newRow.reserve(_cols);
    
    for( size_t col = 0; col < _cols; col++ )
      newRow.push_back(PolyPtr(new Poly( _matrix[col][row]->clone() )));
    
    result._matrix.push_back(newRow);
  }
  
  result._degreeLimitation = _degreeLimitation;
  result._maxDegree = _maxDegree;
  
  result.clean();
  return result;
}

polyjam::core::PolyMatrix
polyjam::core::PolyMatrix::operator+( const PolyMatrix & operant ) const
{  
  PolyMatrix result = clone();
  result += operant;
  return result;
}

polyjam::core::PolyMatrix
polyjam::core::PolyMatrix::operator-( const PolyMatrix & operant ) const
{
  PolyMatrix result = clone();
  result -= operant;
  return result;
}

polyjam::core::PolyMatrix
polyjam::core::PolyMatrix::operator*( const PolyMatrix & operant ) const
{
  //check if the matrices are compatible
  if( cols() != operant.rows() )
  {
    cout << "Error: attempt to multiply incompatible matrices" << endl;
    return PolyMatrix((*this));
  }
  
  //prepare the result matrix
  PolyMatrix result( _matrix[0][0]->zero(), rows(), operant.cols() );
  
  //do actual multiplication here
  for( size_t row = 0; row < rows(); row++ )
  {
    for( size_t col = 0; col < operant.cols(); col++ )
    {
      for( size_t i = 0; i < cols(); i++ )
        result(row,col) += ( *(_matrix[row][i]) * *(operant._matrix[i][col])  );
    }
  }
  
  //transfer other data
  result._degreeLimitation = _degreeLimitation || operant._degreeLimitation;
  result._maxDegree = _maxDegree;
  if(operant._maxDegree < _maxDegree)
    result._maxDegree = operant._maxDegree;
    
  //clean the result
  result.clean();
  return result;
}

polyjam::core::PolyMatrix
polyjam::core::PolyMatrix::operator*( const Poly & operant ) const
{  
  //prepare the result matrix
  PolyMatrix result = clone();
  
  //do actual multiplication here
  for( size_t row = 0; row < rows(); row++ )
  {
    for( size_t col = 0; col < cols(); col++ )
      result(row,col) *= operant;
  }
    
  //clean the result
  result.clean();
  return result;
}

//in-place operations

polyjam::core::PolyMatrix &
polyjam::core::PolyMatrix::negationInPlace()
{
  //prepare the result matrix
  size_t rows = this->rows();
  size_t cols = this->cols();
  
  for( size_t row = 0; row < rows; row++ )
  {    
    for( size_t col = 0; col < cols; col++ )
      _matrix[row][col]->negationInPlace();
  }
  
  //clean the result
  clean();
  return (*this);
}

polyjam::core::PolyMatrix &
polyjam::core::PolyMatrix::transposeInPlace()
{
  //prepare the result matrix
  std::vector< std::vector<PolyPtr> > newMatrix;
  size_t _rows = cols();
  size_t _cols = rows();
  
  newMatrix.reserve(_rows);
  
  for( size_t row = 0; row < _rows; row++ )
  {
    std::vector<PolyPtr> newRow;
    newRow.reserve(_cols);
    
    for( size_t col = 0; col < _cols; col++ )
      newRow.push_back(_matrix[col][row]);
    
    newMatrix.push_back(newRow);
  }
  
  _matrix.clear();
  _matrix = newMatrix;
  
  //clean the result
  clean();
  return (*this);
}

polyjam::core::PolyMatrix &
polyjam::core::PolyMatrix::operator+=( const PolyMatrix & operant )
{
  //check if the matrices are compatible
  if( rows() != operant.rows() || cols() != operant.cols() )
  {
    cout << "Error: attempt to add incompatible matrices" << endl;
    return (*this);
  }
  
  //prepare the result matrix
  size_t rows = this->rows();
  size_t cols = this->cols();
  
  for( size_t row = 0; row < rows; row++ )
  {
    for( size_t col = 0; col < cols; col++ )
      *(_matrix[row][col]) += *(operant._matrix[row][col]);
  }
  
  //transfer other data
  if( operant._degreeLimitation )
    _degreeLimitation = true;
  if(operant._maxDegree < _maxDegree)
    _maxDegree = operant._maxDegree;
  
  //clean the result
  clean();
  return (*this);
}

polyjam::core::PolyMatrix &
polyjam::core::PolyMatrix::operator-=( const PolyMatrix & operant )
{
  //check if the matrices are compatible
  if( rows() != operant.rows() || cols() != operant.cols() )
  {
    cout << "Error: attempt to add incompatible matrices" << endl;
    return (*this);
  }
  
  //prepare the result matrix
  size_t rows = this->rows();
  size_t cols = this->cols();
  
  for( size_t row = 0; row < rows; row++ )
  {
    for( size_t col = 0; col < cols; col++ )
      *(_matrix[row][col]) -= *(operant._matrix[row][col]);
  }
  
  //transfer other data
  if( operant._degreeLimitation )
    _degreeLimitation = true;
  if(operant._maxDegree < _maxDegree)
    _maxDegree = operant._maxDegree;
  
  //clean the result
  clean();
  return (*this);
}

polyjam::core::PolyMatrix &
polyjam::core::PolyMatrix::operator*=( const PolyMatrix & operant )
{
  //check if the matrices are compatible
  if( cols() != operant.rows() )
  {
    cout << "Error: attempt to multiply incompatible matrices" << endl;
    return (*this);
  }
  
  //prepare the result matrix
  std::vector< std::vector<PolyPtr> > newMatrix;
  size_t rows = this->rows();
  size_t cols = operant.cols();
  
  newMatrix.reserve(rows);
  
  for( size_t row = 0; row < rows; row++ )
  {
    std::vector<PolyPtr> newRow;
    newRow.reserve(cols);
    
    for( size_t col = 0; col < cols; col++ )
    {
      newRow.push_back(PolyPtr(
          new Poly( *(_matrix[row][0]) * *(operant._matrix[0][col]) )) );
      
      if( this->cols() > 1 )
      {
        for( size_t i = 0; i < this->cols(); i++ )
          *(newRow[col]) += *(_matrix[row][i]) * *(operant._matrix[i][col]);
      }
    }
    
    newMatrix.push_back(newRow);
  }
  
  _matrix.clear();
  _matrix = newMatrix;
  
  //transfer other data
  if( operant._degreeLimitation )
    _degreeLimitation = true;
  if(operant._maxDegree < _maxDegree)
    _maxDegree = operant._maxDegree;
  
  //clean the result
  clean();
  return (*this);
}

polyjam::core::PolyMatrix &
polyjam::core::PolyMatrix::operator*=( const Poly & operant )
{  
  //prepare the result matrix
  size_t rows = this->rows();
  size_t cols = this->cols();
  
  for( size_t row = 0; row < rows; row++ )
  {    
    for( size_t col = 0; col < cols; col++ )
      *(_matrix[row][0]) *= operant;
  }
  
  //clean the result
  clean();
  return (*this);
}

//special operations

polyjam::core::Poly
polyjam::core::PolyMatrix::dot( const PolyMatrix & operant ) const
{
  //prepare the resulting polynomial
  Poly poly = _matrix[0][0]->zero();
  
  //check if correct input
  if( cols() != 1 || operant.cols() != 1 || rows() != operant.rows() )
  {
    cout << "Error: can only compute dot-product on similar nx1 vectors" << endl;
    return poly;
  }
  
  for( size_t row = 0; row < rows(); row++ )
    poly += ( *(_matrix[row][0]) * *(operant._matrix[row][0]) );
  
  //cleaning
  if( _degreeLimitation || operant._degreeLimitation )
  {
    unsigned int maxDegree = _maxDegree;
    if( operant._maxDegree < maxDegree )
      maxDegree = operant._maxDegree;
    poly.lowerDegreeApproximationInPlace(maxDegree);
  }
  return poly;
}

polyjam::core::PolyMatrix
polyjam::core::PolyMatrix::cross( const PolyMatrix & operant ) const
{
  //check if correct input
  if( cols() != 1 || operant.cols() != 1 || rows() != 3 || operant.rows() != 3 )
  {
    cout << "Error: can only compute cross-product on 3x1 vectors" << endl;
    return PolyMatrix((*this));
  }
  
  // prepare result matrix
  Poly zeroPoly = _matrix[0][0]->zero();
  PolyMatrix result( zeroPoly, 3, 1 );
  
  result[0] =
      (*(_matrix[1][0]) * *(operant._matrix[2][0])) -
      (*(_matrix[2][0]) * *(operant._matrix[1][0]));
  result[1] =
      (*(_matrix[2][0]) * *(operant._matrix[0][0])) -
      (*(_matrix[0][0]) * *(operant._matrix[2][0]));
  result[2] =
      (*(_matrix[0][0]) * *(operant._matrix[1][0])) -
      (*(_matrix[1][0]) * *(operant._matrix[0][0]));
  
  //transfer other data
  result._degreeLimitation = _degreeLimitation || operant._degreeLimitation;
  result._maxDegree = _maxDegree;
  if(operant._maxDegree < _maxDegree)
    result._maxDegree = operant._maxDegree;
  
  result.clean();
  return result;
}

polyjam::core::Poly
polyjam::core::PolyMatrix::determinant() const
{
  //std::cout << "computing degree " << rows() << " determinant" << std::endl;
  
  // check if square
  if( rows() != cols() )
  {
    cout << "Error: Cannot compute determinant of non-suqare matrix" << endl;
    return _matrix[0][0]->zero();
  }
  
  // determinant if 1x1
  if( rows() == 1 )
    return *(_matrix[0][0]);
    
  // determinant if 2x2
  if( rows() == 2 )
  {    
    Poly t1 = *(_matrix[0][0]) * *(_matrix[1][1]);
    Poly b1 = *(_matrix[1][0]) * *(_matrix[0][1]);
    
    return t1 - b1;
  }
  
  //determinant if 3x3
  if( rows() == 3 )
  {
    //std::cout << "part1" << std::endl;
    Poly t1 = *(_matrix[0][0]) * *(_matrix[1][1]);
    if( _degreeLimitation )
      t1.lowerDegreeApproximationInPlace(_maxDegree);
    t1 *= *(_matrix[2][2]);
    if( _degreeLimitation )
      t1.lowerDegreeApproximationInPlace(_maxDegree);
    //std::cout << "part2" << std::endl;
    Poly t2 = *(_matrix[0][1]) * *(_matrix[1][2]);
    if( _degreeLimitation )
      t2.lowerDegreeApproximationInPlace(_maxDegree);
    t2 *= *(_matrix[2][0]);
    if( _degreeLimitation )
      t2.lowerDegreeApproximationInPlace(_maxDegree);
    //std::cout << "part3" << std::endl;
    Poly t3 = *(_matrix[0][2]) * *(_matrix[1][0]);
    if( _degreeLimitation )
      t3.lowerDegreeApproximationInPlace(_maxDegree);
    t3 *= *(_matrix[2][1]);
    if( _degreeLimitation )
      t3.lowerDegreeApproximationInPlace(_maxDegree);
    
    //std::cout << "part4" << std::endl;
    Poly b1 = *(_matrix[2][0]) * *(_matrix[1][1]);
    if( _degreeLimitation )
      b1.lowerDegreeApproximationInPlace(_maxDegree);
    b1 *=  *(_matrix[0][2]);
    if( _degreeLimitation )
      b1.lowerDegreeApproximationInPlace(_maxDegree);
    //std::cout << "part5" << std::endl;
    Poly b2 = *(_matrix[2][1]) * *(_matrix[1][2]);
    if( _degreeLimitation )
      b2.lowerDegreeApproximationInPlace(_maxDegree);
    b2 *= *(_matrix[0][0]);
    if( _degreeLimitation )
      b2.lowerDegreeApproximationInPlace(_maxDegree);
    //std::cout << "part6" << std::endl;
    Poly b3 = *(_matrix[2][2]) * *(_matrix[1][0]);
    if( _degreeLimitation )
      b3.lowerDegreeApproximationInPlace(_maxDegree);
    b3 *= *(_matrix[0][1]);
    if( _degreeLimitation )
      b3.lowerDegreeApproximationInPlace(_maxDegree);
    
    Poly result = t1;
    result += t2;
    result += t3;
    result -= b1;
    result -= b2;
    result -= b3;
    
    //cleaning
    if( _degreeLimitation )
      result.lowerDegreeApproximationInPlace(_maxDegree);
    
    return result;
  }
  
  //this is the generic one if the rows are bigger than 3
  //however, it always proceeds by the first column, unless this code is changed.
  
  std::vector<int> subRows;
  std::vector<int> subCols;
  for( size_t i = 1; i < rows(); i++ )
  {
    subRows.push_back(i);
    subCols.push_back(i);
  }
  
  bool positive = true;
  std::vector< Poly* > polys;
  
  for( size_t i = 0; i < rows(); i++ )
  {
    PolyMatrix sub(_matrix[0][0]->zero(),subRows.size(),subCols.size(),_maxDegree);
    
    for( size_t r = 0; r < subRows.size(); r++ )
    {
      for( size_t c = 0; c < subCols.size(); c++ )
        sub(r,c) = _matrix[subRows[r]][subCols[c]]->clone();
    }
    
    Poly temp = sub.determinant();
    if( !positive )
      temp.negationInPlace();
    
    polys.push_back( new Poly(*(_matrix[i][0]) * temp) );
    if( _degreeLimitation )
      polys.back()->lowerDegreeApproximationInPlace(_maxDegree);
  
    if( i != rows() - 1 )
    {
      subRows[i]--;
    
      if( positive )
        positive = false;
      else
        positive = true;
    }
  }
  
  Poly result = *(polys.front());
  for( size_t i = 1; i < polys.size(); i++ )
  {
    result += *(polys[i]);
    
    if( _degreeLimitation )
      result.lowerDegreeApproximationInPlace(_maxDegree);
  }
  
  for( size_t i = 0; i < polys.size(); i++ )
    delete polys[i];
  
  return result;
}

polyjam::core::Poly
polyjam::core::PolyMatrix::trace() const
{
  // check if square
  if( rows() != cols() )
  {
    cout << "Error: Cannot compute determinant of non-suqare matrix" << endl;
    return _matrix[0][0]->zero();
  }
  
  Poly result = _matrix[0][0]->clone();
  for( size_t i = 1; i < rows(); i++ )
    result += *(_matrix[i][i]);
  
  //cleaning
  if( _degreeLimitation )
    result.lowerDegreeApproximationInPlace(_maxDegree);
  
  return result;
}

polyjam::core::PolyMatrix
polyjam::core::PolyMatrix::skew() const
{
  //check if correct input
  if( cols() != 1 || rows() != 3 )
  {
    cout << "Error: can only compute skew-symmetric of 3x1 vectors" << endl;
    return PolyMatrix((*this));
  }
  
  // prepare result matrix
  PolyMatrix result( _matrix[0][0]->zero(), 3, 3 );
      
  result(0,1) -= *_matrix[2][0];
  result(0,2)  = *_matrix[1][0];
  result(1,2) -= *_matrix[0][0];

  result(1,0)  = *_matrix[2][0];
  result(2,0) -= *_matrix[1][0];
  result(2,1)  = *_matrix[0][0];
  
  //transfer other data
  result._degreeLimitation = _degreeLimitation;
  result._maxDegree = _maxDegree;
  
  result.clean();
  return result;
}

polyjam::core::PolyMatrix
polyjam::core::PolyMatrix::quatMult( const PolyMatrix & operant ) const
{
  //check if correct input
  if( cols() != 1 || operant.cols() != 1 || rows() != 4 || operant.rows() != 4 )
  {
    cout << "Error: can only compute conjugate quaternion of 4x1 vectors" << endl;
    return PolyMatrix((*this));
  }
  
  // prepare result matrix
  PolyMatrix result( _matrix[0][0]->zero(), 4, 1 );
      
  result(0,0) = *_matrix[0][0] * *(operant._matrix[0][0]) - *_matrix[1][0] * *(operant._matrix[1][0]) - *_matrix[2][0] * *(operant._matrix[2][0]) - *_matrix[3][0] * *(operant._matrix[3][0]);
  result(1,0) = *_matrix[1][0] * *(operant._matrix[0][0]) + *_matrix[0][0] * *(operant._matrix[1][0]) + *_matrix[2][0] * *(operant._matrix[3][0]) - *_matrix[3][0] * *(operant._matrix[2][0]);
  result(2,0) = *_matrix[2][0] * *(operant._matrix[0][0]) + *_matrix[0][0] * *(operant._matrix[2][0]) - *_matrix[1][0] * *(operant._matrix[3][0]) + *_matrix[3][0] * *(operant._matrix[1][0]);
  result(3,0) = *_matrix[3][0] * *(operant._matrix[0][0]) + *_matrix[0][0] * *(operant._matrix[3][0]) + *_matrix[1][0] * *(operant._matrix[2][0]) - *_matrix[2][0] * *(operant._matrix[1][0]);
  
  //transfer other data
  result._degreeLimitation = _degreeLimitation;
  result._maxDegree = _maxDegree;
  
  result.clean();
  return result;
}

polyjam::core::PolyMatrix
polyjam::core::PolyMatrix::quatConj() const
{
  //check if correct input
  if( cols() != 1 || rows() != 4 )
  {
    cout << "Error: can only compute conjugate quaternion of 4x1 vectors" << endl;
    return PolyMatrix((*this));
  }
  
  // prepare result matrix
  PolyMatrix result( _matrix[0][0]->zero(), 4, 1 );
      
  result(0,0)  = *_matrix[0][0];
  result(1,0) -= *_matrix[1][0];
  result(2,0) -= *_matrix[2][0];
  result(3,0) -= *_matrix[3][0];
  
  //transfer other data
  result._degreeLimitation = _degreeLimitation;
  result._maxDegree = _maxDegree;
  
  result.clean();
  return result;
}

//comparisons

bool
polyjam::core::PolyMatrix::operator==( const PolyMatrix & operant )
{
  //first compare the size of the matrices
  if( rows() != operant.rows() || cols() != operant.cols() )
    return false;
  
  //now compare the matrix
  for( size_t row = 0; row < rows(); row++ )
  {
    for( size_t col = 0; col < cols(); col++ )
    {
      if( *(_matrix[row][col]) != *(operant._matrix[row][col]) )
        return false;
    }
  }
  
  //if we arrive here, all polynomials are same
  return true;
}

bool
polyjam::core::PolyMatrix::operator!=( const PolyMatrix & operant )
{
  return !((*this) == operant);
}

//handy operations

polyjam::core::PolyMatrix &
polyjam::core::PolyMatrix::setToZero()
{
  for( size_t row = 0; row < rows(); row++ )
  {
    for( size_t col = 0; col < cols(); col++)
      _matrix[row][col]->setToZero();
  }
  
  return (*this);
}

polyjam::core::PolyMatrix &
polyjam::core::PolyMatrix::setToIdentity()
{
  for( size_t row = 0; row < rows(); row++ )
  {
    for( size_t col = 0; col < cols(); col++)
    {
      if( row == col )
        _matrix[row][col]->setToOne();
      else
        _matrix[row][col]->setToZero();
    }
  }
  
  return (*this);
}

//handy comparisons

bool
polyjam::core::PolyMatrix::isZero() const
{  
  for( size_t row = 0; row < rows(); row++ )
  {
    for( size_t col = 0; col < cols(); col++)
    {
      if( !_matrix[row][col]->isZero() )
        return false;
    }
  }
  
  return true;
}

bool
polyjam::core::PolyMatrix::isIdentity() const
{
  for( size_t row = 0; row < rows(); row++ )
  {
    for( size_t col = 0; col < cols(); col++)
    {
      if( row == col )
      {
        if( !_matrix[row][col]->isOne() )
          return false;
      }
      else
      {
        if( !_matrix[row][col]->isZero() )
          return false;
      }
    }
  }
  
  return true;
}

//private

void
polyjam::core::PolyMatrix::clean()
{
  //degree limitation
  if( _degreeLimitation )
  {    
    //now do actual degree limitation here
    for( size_t row = 0; row < rows(); row++ )
    {
      for( size_t col = 0; col < cols(); col++ )
        _matrix[row][col]->lowerDegreeApproximationInPlace( _maxDegree );
    }
  }
}
