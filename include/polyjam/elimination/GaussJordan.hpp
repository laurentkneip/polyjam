/******************************************************************************
 * Author:   Laurent Kneip                                                    *
 * Contact:  kneip.laurent@gmail.com                                          *
 * License:  Copyright (c) 2013 Laurent Kneip, ANU. All rights reserved.      *
 *                                                                            *
 * Redistribution and use in source and binary forms, with or without         *
 * modification, are permitted provided that the following conditions         *
 * are met:                                                                   *
 * * Redistributions of source code must retain the above copyright           *
 *   notice, this list of conditions and the following disclaimer.            *
 * * Redistributions in binary form must reproduce the above copyright        *
 *   notice, this list of conditions and the following disclaimer in the      *
 *   documentation and/or other materials provided with the distribution.     *
 * * Neither the name of ANU nor the names of its contributors may be         *
 *   used to endorse or promote products derived from this software without   *
 *   specific prior written permission.                                       *
 *                                                                            *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"*
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE  *
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE *
 * ARE DISCLAIMED. IN NO EVENT SHALL ANU OR THE CONTRIBUTORS BE LIABLE        *
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL *
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR *
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER *
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT         *
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY  *
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF     *
 * SUCH DAMAGE.                                                               *
 ******************************************************************************/

#ifndef POLYJAM_ELIMINATION_GAUSSJORDAN_HPP_
#define POLYJAM_ELIMINATION_GAUSSJORDAN_HPP_

#include <stdlib.h>
#include <stdio.h>
#include <list>

#include <opencv2/opencv.hpp>
#include <polyjam/core/Coefficient.hpp>

namespace polyjam
{

namespace elimination
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
    bool saveImage = false )
{
  //define some types to make life easier
  typedef COEFFICIENT coefficient_t;
  
  //get the zero constant
  coefficient_t zero = customGetZero((*(matrix[0]))[0]);
  coefficient_t precision = customGetPrecision( (*(matrix[0]))[0] );
  
  //define some colors
  CvScalar Red = CV_RGB( 255, 0, 0 );
  CvScalar DarkRed = CV_RGB( 220, 0, 0 );
  //CvScalar Yellow = CV_RGB( 255, 255, 0 );
  CvScalar Blue = CV_RGB( 0, 0, 255 );
  //CvScalar Purple = CV_RGB( 192, 0, 192 );

  //create a white image with correct size
  int imgWidth = matrix.front()->size();
  int imgHeight = matrix.size();
  
  IplImage* img =
      cvCreateImage( cvSize( imgWidth, imgHeight ), IPL_DEPTH_8U, 3 );
  
  //set the whole image with the right color
  for( int row = 0; row < imgHeight; row++ )
  {
    for( int col = 0; col < imgWidth; col++ )
    {
      CvScalar color = Blue;
      if( !( (*(matrix[row]))[col] == zero ) )
      //***//if(
      //***//    ( (*(matrix[row]))[col] < zero && (*(matrix[row]))[col] < precision.negation() ) ||
      //***//    (*(matrix[row]))[col] > precision )
      {
        color = Red;
        if( col % 2 == 1 )
          color = DarkRed;
      }
      
      cvSet2D( img, row, col, color );
    }
  }
  
  //now visualize the image
  std::string name("coefficient matrix");
  
  if( createAndDestroy )
  {
    cvNamedWindow( name.c_str(), 0 );
    cvStartWindowThread();
  }
  
  cvResizeWindow( name.c_str(), imgWidth, imgHeight );
  cvShowImage( name.c_str(), img );
  if( saveImage )
  {
    std::string path("/Users/laurent/temp/pe_pattern.png");
    cvSaveImage( path.c_str(), img );
  }
  
  if( createAndDestroy )
  {
    char exit_key_press = 0;
    while ((int) exit_key_press != 27) // or key != ESC
      exit_key_press = cvWaitKey(10);
    //cvDestroyWindow(name.c_str());
  }
  
  cvReleaseImage(&img);
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
    cvNamedWindow( name.c_str(), 0 );
    cvStartWindowThread();
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
    cvDestroyWindow(name.c_str());
}

}
}

#endif /* POLYJAM_ELIMINATION_GAUSSJORDAN_HPP_ */
