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

#include "GaussJordan.hpp"


void
polyjam::math::gaussReduction( Eigen::MatrixXd & matrix )
{  
  //create some constant for zero-checking etc.
  double precision = 0.0000000001;

  int rows = matrix.rows();
  int cols = matrix.cols();
  int maxIndentation = cols-1;
  
  //first compute the number of steps to do
  int steps = rows;
  if( rows > cols )
    steps = cols;
  
  //working variables
  int currentIndentation = 0;
  int frontRow = 0;
  
  //first step down
  for( int dummy = 0; dummy < steps; dummy++ )
  {
    // first iterate through the rows and find a row that has a
    // non-zero coefficient in the relevant column
    
    int row;
    
    while(currentIndentation <= maxIndentation)
    {
      double maxValue = -1.0;
      
      for( int tempRow = frontRow; tempRow < rows; tempRow++ )
      {
        double value = fabs(matrix(tempRow,currentIndentation));
        if( value > maxValue )
        {
          row = tempRow;
          maxValue = value;
        }
      }
      
      //if, after checking this column, we found a non-zero element, we are done
      if(maxValue > precision)
        break;
      
      //if not, we have to move on to the next column
      currentIndentation++;
    }
    
    //if we are beyond the maxIndentation, the whole rest is zero. Break!
    if(currentIndentation > maxIndentation)
      break;
    
    //row is now the row that should go in the place of frontRow->swap
    Eigen::MatrixXd rowCopy = matrix.row(row);
    matrix.row(row) = matrix.row(frontRow);
    matrix.row(frontRow) = rowCopy;
    
    //ok, now use frontRow!
    int col = currentIndentation;
    double leadingCoefficient = matrix(frontRow,col);
    matrix.row(frontRow) /= leadingCoefficient;
    
    //iterate through all remaining rows, and subtract correct multiple of 
    //first row (if leading coefficient is non-zero!)
    row = frontRow;
    ++row;
    while( row < rows )
    {
      double leadingCoefficient = matrix(row,currentIndentation);
      if( fabs(leadingCoefficient) > precision )
        matrix.row(row) -= leadingCoefficient * matrix.row(frontRow);
      ++row;
    }
    
    //increment row and currentIndentation
    ++frontRow;
    currentIndentation++;
    if( currentIndentation > maxIndentation )
      break;
  }
  
  //if the following, this is the zero matrix -> nothing to do!
  if( frontRow == 0 )
    return;
  
  //actually, everything starting from frontRow should be zero -> delete
  Eigen::MatrixXd croppedMatrix = matrix.block(0,0,frontRow,cols);
  matrix = croppedMatrix;
  
  //set index to the last non-zero row
  --frontRow;
  
  //Now step up
  while( frontRow != 0 )
  {    
    //indent until we find a non-zero element in frontRow
    int indentations = 0;
    while( fabs(matrix(frontRow,indentations)) < precision && indentations < cols )
      indentations++;
    
    //if the following, there is a problem, maybe zero-matrix!
    if( indentations == cols )
      break;
    
    //get the working row
    int row = frontRow;
    
    do
    {      
      //decrement working row
      --row;
      
      //working column
      
      //now get the leading coefficient
      double leadingCoefficient = matrix(row,indentations);
      
      //Now iterator until the end, and subtract each time the multiplied
      //front-row
      if( fabs(leadingCoefficient) > precision )
        matrix.row(row) -= leadingCoefficient * matrix.row(frontRow);
    }
    while( row != 0 );
    
    --frontRow;
  }
}