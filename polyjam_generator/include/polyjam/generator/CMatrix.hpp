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
 
#ifndef POLYJAM_GENERATOR_CMATRIX_HPP_
#define POLYJAM_GENERATOR_CMATRIX_HPP_

#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <list>

#include <polyjam/core/Poly.hpp>


/**
 * \brief The namespace of this library.
 */
namespace polyjam
{

/**
 * \brief The namespace for the elimination templates.
 */
namespace generator
{

/**
 * CMatrix
 */
class CMatrix
{

public:
  typedef std::list<core::Poly*> polynomials_t;
  typedef std::vector<core::Monomial> monomials_t;
  
  typedef std::vector<core::Coefficient> crow_t;
  typedef std::vector<crow_t*> cmatrix_t;
  
  typedef std::pair<int,core::Monomial> eq_t;
  typedef std::vector<eq_t> eqs_t;
  
  //constructors
  CMatrix( const polynomials_t & polynomials );
  CMatrix( const polynomials_t & polynomials, const eqs_t & equations );
  CMatrix( const polynomials_t & polynomials, const monomials_t & order );
  CMatrix( const polynomials_t & polynomials, const monomials_t & order, const eqs_t & equations );
  
  //destructor
  virtual ~CMatrix();

  //modifiers
  void reduce();
  
  //accessors
  size_t rows();
  size_t cols();
  core::Coefficient operator()( size_t row, size_t col );
  
  CMatrix subMatrix( const std::list<int> & rows );
  monomials_t monomials();
  core::Poly getPolynomial( int row );
  core::Poly getSymbolicPolynomial( int row, const std::string & matrixName ); 
  core::Poly getSymbolicPolynomial2( int row );
  polynomials_t getPolynomials();
  polynomials_t getSymbolicPolynomials( const std::string & matrixName );
  polynomials_t getSymbolicPolynomials2();
  
  //verifiers
  bool contains( const polynomials_t & polynomials );
  void visualize();
  void save( const std::string & name, const std::string & save_path );

private:
  CMatrix( cmatrix_t & matrix, monomials_t & monomials );

  void fillMonomials( const polynomials_t & polynomials );
  void fillMatrix( const polynomials_t & polynomials, bool quickOrdering = true );

  cmatrix_t _matrix;
  monomials_t _monomials;
};

}
}

#endif /* POLYJAM_GENERATOR_CMATRIX_HPP_ */
