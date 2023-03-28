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

#ifndef POLYJAM_MATH_GAUSSJORDAN_HPP_
#define POLYJAM_MATH_GAUSSJORDAN_HPP_

#include <stdlib.h>
#include <stdio.h>
#include <list>

#include <opencv2/opencv.hpp>
#include <polyjam/core/Coefficient.hpp>

namespace polyjam
{

namespace math
{

//some stuff where the signature changes between types
void customPrint( double value );
void customPrint( const core::Coefficient & value );
double customGetZero( double value );
double customGetOne( double value );
double customGetPrecision( double value );
core::Coefficient customGetZero( core::Coefficient value );
core::Coefficient customGetOne( const core::Coefficient & value );
core::Coefficient customGetPrecision( const core::Coefficient & value );

//print the matrix to the console
template<typename COEFFICIENT>
void
printMatrix( std::vector<std::vector<COEFFICIENT>*> matrix )
{  
  for( int row = 0; row < matrix.size(); row++ )
  {
    for( int col = 0; col < matrix[0]->size(); col++ )
      customPrint( (*(matrix[row]))[col] );
    std::cout << std::endl << std::endl;
  }
  
  std::cout << std::endl;
}

template<typename COEFFICIENT>
void
visualizeMatrix(
    std::vector<std::vector<COEFFICIENT>*> & matrix,
    bool createAndDestroy = true,
    bool saveImage = false,
    std::string name = std::string("coefficientMatrix"),
    std::string save_path = std::string("/tmp") )
{
  //define some types to make life easier
  typedef COEFFICIENT coefficient_t;
  
  //get the zero constant
  coefficient_t zero = customGetZero((*(matrix[0]))[0]);
  coefficient_t precision = customGetPrecision( (*(matrix[0]))[0] );
  
  //define some colors
  cv::Scalar Red( 255, 0, 0 );
  cv::Scalar DarkRed( 220, 0, 0 );
  //cv::Scalar Yellow( 255, 255, 0 );
  cv::Scalar Blue( 0, 0, 255 );
  //cv::Scalar Purple( 192, 0, 192 );

  //create a white image with correct size
  int imgWidth = matrix.front()->size();
  int imgHeight = matrix.size();
  
  cv::Mat img( imgHeight, imgWidth, CV_8UC3, cv::Scalar( 0, 0, 0 ) );
  
  //set the whole image with the right color
  for( int row = 0; row < imgHeight; row++ )
  {
    for( int col = 0; col < imgWidth; col++ )
    {
      cv::Scalar color = Blue;
      if( !( (*(matrix[row]))[col] == zero ) )
      //***//if(
      //***//    ( (*(matrix[row]))[col] < zero && (*(matrix[row]))[col] < precision.negation() ) ||
      //***//    (*(matrix[row]))[col] > precision )
      {
        color = Red;
        if( col % 2 == 1 )
          color = DarkRed;
      }
      img.at<cv::Vec3b>(row, col)[0] = color[0];
      img.at<cv::Vec3b>(row, col)[1] = color[1];
      img.at<cv::Vec3b>(row, col)[2] = color[2];
    }
  }
  
  //now visualize the image
  if( saveImage )
  {
    cv::imwrite(save_path + name + std::string(".png"), img);
  }
}

//This is a templated implementation of a Gauss-reduction
//The template parameter COEFFICIENT only has to support standard operations,
//or otherwise implement one of the small helpers at the top of this file
template<typename COEFFICIENT>
void
gaussReduction(
    std::vector<std::vector<COEFFICIENT>*> & matrix,
    bool continuousVisualization = false )
{
  //define some types to make life easier
  typedef COEFFICIENT coefficient_t;

  bool consolePrint = false;
  
  //initialize the continuous visualization
  //this is useful when implementing new improvements
  std::string name("coefficient matrix");
  if( continuousVisualization )
  {
    cv::namedWindow( name, 0 );
    cv::startWindowThread();
  }
  
  //create some constants for zero-checking etc.
  coefficient_t temp = (*(matrix[0]))[0];
  coefficient_t one = customGetOne(temp);
  coefficient_t zero = customGetZero(temp);
  coefficient_t precision = customGetPrecision(temp);
  int rows = matrix.size();
  int cols = matrix[0]->size();
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
    if( consolePrint )
      std::cout << "dummy is " << dummy << "/" << steps << ", indentation is " << currentIndentation << "/" << maxIndentation << std::endl;
    // first iterate through the rows and find a row that has a
    // non-zero coefficient in the relevant column
    
    int row;
    
    bool nonzeroFound = false;
    while(!nonzeroFound)
    {
      for( int tempRow = frontRow; tempRow < rows; tempRow++ )
      {        
        if( (*(matrix[tempRow]))[currentIndentation] != zero )
        {
          row = tempRow;
          nonzeroFound = true;
          break;
        }
      }
      
      //if, after checking this column, we found a non-zero element, we are done
      if(nonzeroFound)
        break;
      
      //if not, we have to move on to the next column
      currentIndentation++;
      if( currentIndentation > maxIndentation )
        break;
    }
    
    //if we are beyond the maxIndentation, the whole rest is zero. Break!
    if(currentIndentation > maxIndentation)
      break;
    
    //rowIter is now the row that should go in the place of frontRow->swap    
    std::swap( matrix[row], matrix[frontRow] );
    
    //ok, now use frontRow!
    
    //first divide all coefficients by the leading coefficient
    int col = currentIndentation;
    coefficient_t leadingCoefficient = (*(matrix[frontRow]))[col] + zero; //the +zero is important! It causes the coefficient to be copied!
    (*(matrix[frontRow]))[col] = one + zero;
    coefficient_t leadingCoefficient_inv = (*(matrix[frontRow]))[col] / leadingCoefficient;
    ++col;
    while( col < cols )
    {
      (*(matrix[frontRow]))[col] *= leadingCoefficient_inv;
      ++col;
    }
    
    //create a vector of bools indicating the cols that need to be manipulated
    col = currentIndentation;
    std::vector<int> nonzeroIdx;
    nonzeroIdx.reserve( cols - currentIndentation );
    while( col < cols )
    {
      if( ((*(matrix[frontRow]))[col]) != zero )
        nonzeroIdx.push_back(col);
      col++;
    }
    
    //iterate through all remaining rows, and subtract correct multiple of 
    //first row (if leading coefficient is non-zero!)
    row = frontRow;
    ++row;
    while( row < rows )
    {
      col = currentIndentation;
      coefficient_t leadingCoefficient = (*(matrix[row]))[col] + zero;
      
      if( leadingCoefficient != zero )
      {
        for( int col = 0; col < (int) nonzeroIdx.size(); col++ )
        {
          (*(matrix[row]))[nonzeroIdx[col]] -= leadingCoefficient * ((*(matrix[frontRow]))[nonzeroIdx[col]]);
          //***//if(
          //***//    ((*(matrix[row]))[nonzeroIdx[col]] < zero && (*(matrix[row]))[nonzeroIdx[col]].negation() < precision ) ||
          //***//    ((*(matrix[row]))[nonzeroIdx[col]] > zero && (*(matrix[row]))[nonzeroIdx[col]] < precision ) )
          //***//  (*(matrix[row]))[nonzeroIdx[col]] = zero + zero;
        }
      }
      
      ++row;
    }
    
    //increment row and currentIndentation
    ++frontRow;
    currentIndentation++;
    if( currentIndentation > maxIndentation )
      break;
    
    //visualize if desired
    if(continuousVisualization)
      visualizeMatrix( matrix, false );
  }
  
  //if the following, this is the zero matrix -> nothing to do!
  if( frontRow == 0 )
    return;
  
  //actually, everything starting from frontRow should be zero -> delete
  if( frontRow < rows )
  {
    for( int tempRow = frontRow; tempRow < rows; tempRow++ )
      delete matrix[tempRow];
    matrix.resize(frontRow);
  }
  
  //set index to the last non-zero row
  --frontRow;
  
  //Now step up
  while( frontRow != 0 )
  {
    if(consolePrint)
      std::cout << "frontRow is " << frontRow << "/" << rows << std::endl;
    
    //indent until we find a non-zero element in frontRow
    int indentations = 0;
    while( (*(matrix[frontRow]))[indentations] == zero && indentations < cols )
      indentations++;
    
    //if the following, there is a problem, maybe zero-matrix!
    if( indentations == cols )
      break;
    
    //create a vector of bools indicating the cols that need to be manipulated
    int col = indentations;
    std::vector<int> nonzeroIdx;
    nonzeroIdx.reserve( cols - indentations );
    while( col < cols )
    {
      if( ((*(matrix[frontRow]))[col]) != zero )
        nonzeroIdx.push_back(col);
      col++;
    }
    
    //get the working row
    int row = frontRow;
    
    do
    {      
      //decrement working row
      --row;
      
      //working column
      
      //now get the leading coefficient
      coefficient_t leadingCoefficient = (*(matrix[row]))[indentations] + zero;
      
      //Now iterator until the end, and subtract each time the multiplied
      //front-row
      if( leadingCoefficient != zero )
      {        
        for( int col = 0; col < (int) nonzeroIdx.size(); col++ )
          (*(matrix[row]))[nonzeroIdx[col]] -= leadingCoefficient * (*(matrix[frontRow]))[nonzeroIdx[col]];
      }
    }
    while( row != 0 );
    
    //visualize if desired
    if(continuousVisualization)
      visualizeMatrix( matrix, false );
    
    --frontRow;
  }
  
  if(continuousVisualization)
    cv::destroyWindow(name);
}

}
}

#endif /* POLYJAM_MATH_GAUSSJORDAN_HPP_ */
