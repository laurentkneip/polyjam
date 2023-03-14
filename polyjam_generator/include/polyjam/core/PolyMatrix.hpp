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

#ifndef POLYJAM_CORE_POLYMATRIX_HPP_
#define POLYJAM_CORE_POLYMATRIX_HPP_

#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <memory>
#include <polyjam/core/Poly.hpp>

//todo:
//-maybe add PolyMatrix zero() const;
//-maybe add PolyMatrix one() const;

/**
 * \brief The namespace of this library.
 */
namespace polyjam
{

/**
 * \brief The namespace of the core objects of polynomials
 */
namespace core
{

/**
 * PolyMatrix defines a matrix-polynomial. It allows to access
 * its elements, and do matrix-operations. Polynomials can only be replaced by
 * different ones via the matrix assignment operator, which is why the matrix
 * stays consistent. Consistency with operants is implicitly checked by the
 * polynomial operations.
 */
class PolyMatrix
{
public:
  /** A pointer to a polynomial */
  typedef std::shared_ptr<Poly> PolyPtr;
  
  //constructors/destructor

  /**
   * \brief Constructor for a nxm matrix.
   * \param[in] base The initial value of the elements.
   * \param[in] rows The number of rows.
   * \param[in] cols The number of columns.
   * \param[in] diag Only fill diagonal (and set rest to zero)?
   */
  PolyMatrix( const Poly & base, size_t rows, size_t cols, bool diag = false );
  /**
   * \brief Constructor for a nxm matrix.
   * \param[in] base The initial value of the elements.
   * \param[in] rows The number of rows.
   * \param[in] cols The number of columns.
   * \param[in] maxDegree Maximum degree of automatic degree-limitation.
   * \param[in] diag Only fill diagonal (and set rest to zero)?
   */
  PolyMatrix(
      const Poly & base,
      size_t rows,
      size_t cols,
      unsigned int maxDegree,
      bool diag = false );

  /**
   * \brief Destructor.
   */
  ~PolyMatrix();
  
  // deep stuff
  void copy( const PolyMatrix & copy );
  PolyMatrix clone() const;
  
  //output/access
  
  /**
   * \brief Retrieve the number of rows in this matrix.
   * \return The number of rows.
   */
  size_t rows() const;
  /**
   * \brief Retrieve the number of columns in this matrix.
   * \return The number of columns.
   */
  size_t cols() const;
  /**
   * \brief Access an element of the matrix (by reference).
   * \param[in] row The row of the element.
   * \param[in] col The column of the element.
   * \return A reference to the element (polynomial).
   */
  Poly & operator()( size_t row, size_t col );
  /**
   * \brief Row-major access to an element (by reference). Useful for vectors.
   * \param[in] index The index of the element.
   * \return A reference to the element (polynomial).
   */
  Poly & operator[]( size_t index );
  
  //operations
  
  /**
   * \brief Compute the negation of this matrix.
   * \return The negation of this matrix.
   */
  PolyMatrix negation() const;
  /**
   * \brief Compute the transpose of this matrix.
   * \return The transpose of this matrix.
   */
  PolyMatrix transpose() const;
  /**
   * \brief Compute the sum of this and another matrix. Dimensions must agree.
   * \param[in] operant The summand.
   * \return The sum of this and operant.
   */
  PolyMatrix operator+( const PolyMatrix & operant ) const;
  /**
   * \brief Compute the different between this and another matrix. Dimensions
   *        must agree.
   * \param[in] operant The other matrix.
   * \return The different between this and operant.
   */
  PolyMatrix operator-( const PolyMatrix & operant ) const;
  /**
   * \brief Compute the matrix-product of this and another matrix. The number of
   *        columns of this matrix must be equal to the number of rows of the
   *        operant.
   * \param[in] operant The multiplicant.
   * \return The matrix-product of this and operant.
   */
  PolyMatrix operator*( const PolyMatrix & operant ) const;
  /**
   * \brief Compute the matrix-product of this and a polynomial. Each element
   *        will be multiplied.
   * \param[in] operant The multiplicant.
   * \return The matrix-product of this and operant.
   */
  PolyMatrix operator*( const Poly & operant ) const;
  
  //in-place operations
  
  /**
   * \brief Negate this matrix.
   * \return A reference to this matrix (updated).
   */
  PolyMatrix & negationInPlace();
  /**
   * \brief Transpose this matrix.
   * \return A reference to this matrix (updated).
   */
  PolyMatrix & transposeInPlace();
  /**
   * \brief Compute the sum of this and another matrix. Dimensions must agree.
   *        In-place operation.
   * \param[in] operant The summand.
   * \return A reference to this matrix (updated).
   */
  PolyMatrix & operator+=( const PolyMatrix & operant );
  /**
   * \brief Compute the different between this and another matrix. Dimensions
   *        must agree. In-place operation.
   * \param[in] operant The other matrix.
   * \return A reference to this matrix (updated).
   */
  PolyMatrix & operator-=( const PolyMatrix & operant );
  /**
   * \brief Compute the matrix-product of this and another matrix. The number of
   *        columns of this matrix must be equal to the number of rows of the
   *        operant. In-place operation.
   * \param[in] operant The multiplicant.
   * \return A reference to this matrix (updated).
   */
  PolyMatrix & operator*=( const PolyMatrix & operant );
  /**
   * \brief Compute the matrix-product of this and a polynomial. Each element
   *        will be multiplied.
   * \param[in] operant The multiplicant.
   * \return A reference to this matrix (updated).
   */
  PolyMatrix & operator*=( const Poly & operant );
  
  //special operations
  
  /**
   * \brief Compute the dot-product of this and another matrix. The matrices
   *        must be column-vectors of same dimension.
   * \param[in] operant The operant.
   * \return The dot-product of this and another matrix.
   */
  Poly dot( const PolyMatrix & operant ) const;
  /**
   * \brief Compute the cross-product of this and another matrix. The matrices
   *        must be 3-vectors.
   * \param[in] operant The operant.
   * \return 3-vector representing the cross-product.
   */
  PolyMatrix cross( const PolyMatrix & operant ) const;
  /**
   * \brief Compute the determinant of this matrix. The matrix must be either
   *        1x1, 2x2, or 3x3.
   * \return The determinant polynomial.
   */
  Poly determinant() const;
  /**
   * \brief Compute the trace of this matrix.
   * \return The trace polynomial of this matrix.
   */
  Poly trace() const;
  /**
   * \brief Take the skew-symmetric of a 3-vector.
   * \return The skew-symmetric.
   */
  PolyMatrix skew() const;
  /**
   * \brief Treat 4-vectors like quaternions and perform multiplication.
   * \param[in] operant The operant.
   * \return The product of the two quaternions.
   */
  PolyMatrix quatMult( const PolyMatrix & operant ) const;
  /**
   * \brief Get the conjugate complex of a quaternion.
   * \return The conjugate complex of this quaternion.
   */
  PolyMatrix quatConj() const;
  

  //comparisons
  
  /**
   * \brief Compare this matrix to another one.
   * \param[in] operant The comparison matrix.
   * \return All members plus dimensions equal?
   */
  bool operator==( const PolyMatrix & operant );
  /**
   * \brief Compare this matrix to another one.
   * \param[in] operant The comparison matrix.
   * \return Any member or dimensions different?
   */
  bool operator!=( const PolyMatrix & operant );
  
  //handy operations
  
  /**
   * \brief Set all elements of this matrix to zero.
   * \return A reference to this matrix.
   */
  PolyMatrix & setToZero();
  /**
   * \brief Set this matrix to identity.
   * \return A reference to this matrix.
   */
  PolyMatrix & setToIdentity();

  //handy comparisons
  
  /**
   * \brief Are all the members of this matrix zero?
   * \return All members zero?
   */
  bool isZero() const;
  /**
   * \brief Does this matrix equal to identity?
   * \return Is this the identity matrix?
   */
  bool isIdentity() const;
  
private:
  /** The matrix elements (pointers to polynomials) */
  std::vector< std::vector<PolyPtr> > _matrix;
  /** Automatic degree limitation active? */
  bool _degreeLimitation;
  /** Maximum degree for result of any operation */
  unsigned int _maxDegree;
  
  /**
   * \brief Effectuate the degree-limitation on this polynomial.
   */
  void clean();
  
public:
  /** Useful named constructor idioms */
  
  // zero matrices
  static PolyMatrix  zeroR( size_t dimensions, size_t rows, size_t cols, Monomial::Order order = Monomial::GREVLEX )
  { return PolyMatrix( Poly::zeroR(dimensions,order),rows,cols); }
  static PolyMatrix  zeroQ( size_t dimensions, size_t rows, size_t cols, Monomial::Order order = Monomial::GREVLEX )
  { return PolyMatrix( Poly::zeroQ(dimensions,order),rows,cols); }
  static PolyMatrix  zeroZ( size_t dimensions, size_t rows, size_t cols, Monomial::Order order = Monomial::GREVLEX )
  { return PolyMatrix( Poly::zeroZ(dimensions,order),rows,cols); }
  static PolyMatrix  zeroS( size_t dimensions, size_t rows, size_t cols, Monomial::Order order = Monomial::GREVLEX )
  { return PolyMatrix( Poly::zeroS(dimensions,order),rows,cols); }
  static PolyMatrix zeroSZ( size_t dimensions, size_t rows, size_t cols, Monomial::Order order = Monomial::GREVLEX )
  { return PolyMatrix(Poly::zeroSZ(dimensions,order),rows,cols); }
  
  //identity matrices
  static PolyMatrix  idR( size_t dimensions, size_t rows, size_t cols, Monomial::Order order = Monomial::GREVLEX )
  { return PolyMatrix( Poly::oneR(dimensions,order),rows,cols,true); }
  static PolyMatrix  idQ( size_t dimensions, size_t rows, size_t cols, Monomial::Order order = Monomial::GREVLEX )
  { return PolyMatrix( Poly::oneQ(dimensions,order),rows,cols,true); }
  static PolyMatrix  idZ( size_t dimensions, size_t rows, size_t cols, Monomial::Order order = Monomial::GREVLEX )
  { return PolyMatrix( Poly::oneZ(dimensions,order),rows,cols,true); }
  static PolyMatrix  idS( size_t dimensions, size_t rows, size_t cols, Monomial::Order order = Monomial::GREVLEX )
  { return PolyMatrix( Poly::oneS(dimensions,order),rows,cols,true); }
  static PolyMatrix idSZ( size_t dimensions, size_t rows, size_t cols, Monomial::Order order = Monomial::GREVLEX )
  { return PolyMatrix(Poly::oneSZ(dimensions,order),rows,cols,true); }
};

}
}

#endif /* POLYJAM_CORE_POLYMATRIX_HPP_ */
