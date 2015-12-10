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
 
#ifndef POLYJAM_ELIMINATION_CMATRIX_HPP_
#define POLYJAM_ELIMINATION_CMATRIX_HPP_

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
namespace elimination
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
  void visualize( bool save = false );

private:
  CMatrix( cmatrix_t & matrix, monomials_t & monomials );

  void fillMonomials( const polynomials_t & polynomials );
  void fillMatrix( const polynomials_t & polynomials, bool quickOrdering = true );

  cmatrix_t _matrix;
  monomials_t _monomials;
};

}
}

#endif /* POLYJAM_ELIMINATION_CMATRIX_HPP_ */
