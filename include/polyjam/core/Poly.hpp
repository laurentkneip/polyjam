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

#ifndef POLYJAM_CORE_POLY_HPP_
#define POLYJAM_CORE_POLY_HPP_

#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <set>
#include <boost/shared_ptr.hpp>
#include <polyjam/core/Term.hpp>

//todo:
//-possibly include the polynomial division here
//-add Poly operator/( const Poly & operant ) const;
//-add Poly operator/( const Term & operant ) const;

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
 * Poly defines a polynomial. The idea is that polynomials are always initialized
 * with exactly one term. This way, they can clearly represent zero, and it
 * allows reusing the handy constructors for the terms. The ordering, field,
 * degree, multiplicity, is all given by the first term in the list of terms.
 * Polynomials therefore need to make sure that the list remains "clean",
 * meaning that no terms of different "character" are mixed. Not having an
 * "empty" polynomial is also exploited in the assignment operator ensuring that
 * the assigned polynomial remains of the same type, which in turn is handy to
 * ensure consistency of PolyMatrix.
 *
 * Careful: constructions and copying are lazy. Use copy/clone if more is required!
 */
class Poly
{

public:

  /**
   * The container for the terms. We use a set, which ensured lookup and
   * insertion in exponential time.
   */
  typedef std::set<Term,std::greater<Term> > terms_t;
  /** A pointer to the terms-container */
  typedef boost::shared_ptr<terms_t> termsPtr;

  // constructors/destructor
  
  /**
   * \brief The main constructor. Constructs a polynomial with one term.
   *        This can be a zero term as well.
   * \param[in] term The initial term.
   */
  Poly( const Term & term );
  /**
   * \brief Destructor.
   */
  virtual ~Poly();
  
  //assignment
  
  /**
   * \brief Ensures that the polynomial remains of the same type,
   *        meaning dimension, field(s), ordering. Sugar is copied!
   * \param[in] rhs The right-hand side of "=".
   * \return A reference to this (updated) polynomial (lazy transfer!).
   */
  Poly & operator=( const Poly & rhs );
  
  //output  
  
  /**
   * \brief Create a clone of this polynomial (deep). Sugar is not copied!
   * \param[in] full Maintain multiplicity?
   * \return A clone with the desired coefficients.
   */
  Poly clone( bool full = true ) const;
  /**
   * \brief Copy the content from another polynomial (deep). Sugar is copied.
   * \param[in] copy The original.
   */
  void copy( const Poly & copy );
  /**
   * \brief Create a polynomial that has similar characteristics, but is one.
   * \return A polynomial with similar characteristics, but one.
   */
  Poly one() const;
  /**
   * \brief Create a polynomial that has similar characteristics, but is zero.
   * \return A polynomial with similar characteristics, but zero.
   */
  Poly zero() const;
  /**
   * \brief Set the dominant coefficient in each term.
   * \param[in] index The index of the dominant coefficient.
   */
  void setDominant( int index );
  /**
   * \brief Print this polynomial on the console.
   */
  void print() const;
  /**
   * \brief Get a string that represents this polynomial.
   * \param[in] c_version Get the string in C++-syntax using pow(...,...)?
   * \return The string representing this polynomial.
   */
  std::string getString( bool c_version ) const;
  /**
   * \brief Get a reference to the leading term. Only a const-reference!
   * \return A reference to the leading term.
   */   
  const Term & leadingTerm() const;
  /**
   * \brief Get the leading coefficient as a new term.
   * \return A term which represents the leading coefficient.
   */
  Term leadingCoefficient() const;
  /**
   * \brief Get the leading monomial as a new term.
   * \return A term which represents the leading monomial.
   */
  Term leadingMonomial() const;
  /**
   * \brief Get an iterator through the polynomial terms pointing at the front.
   * \return An iterator pointing at the front term.
   */
  terms_t::iterator begin() const;
  /**
   * \brief Get an iterator through the polynomial terms pointing at the end.
   * \return An iterator pointing at the end term.
   */
  terms_t::iterator end() const;
  /**
   * \brief Get the number of terms in this polynomial.
   * \return The number of terms in this polynomial.
   */
  size_t size() const;
  /**
   * \brief Access the sugar of this polynomial (used in Buchberger)
   * \return a reference to the sugar of this polynomial.
   */
  unsigned int & sugar();
  /**
   * \brief Evaluate this polynomial by replacing all unknowns with the
   *        provided values.
   * \param[in] values The replacement values for each unknown.
   * \return A Coefficient of R that has the value of the evaluated polynomial.
   */
  Coefficient eval( const std::vector<double> & values ) const;
  /**
   * \brief Evaluate this polynomial by replacing all unknowns with the
   *        provided values.
   * \param[in] calues The replacement values for each unknown.
   * \return A Coefficient that has the value of the evaluated polynomial.
   */
  Coefficient eval( const std::vector<Coefficient> & values ) const;
  
  //operations
  
  /**
   * \brief Extract a version of this polynomial that has different order.
   * \param[in] newOrder The desired monomial order.
   * \return The newly ordered polynomial.
   */
  Poly differentOrderVersion( Monomial::Order newOrder ) const;
  /**
   * \brief Generate a lower-degree approximation of this polynomial.
   * \param[in] maxDegree The maximum degree in the output polynomial.
   * \return The lower-degree approximation.
   */
  Poly lowerDegreeApproximation( unsigned int maxDegree ) const;
  /**
   * \brief Get the negative of this polynomial.
   * \return The negative of this polynomial.
   */
  Poly negation() const;
  /**
   * \brief Compute the sum of this and another polynomial.
   * \param[in] operant The other polynomial.
   * \return The sum of this and another polynomial.
   */
  Poly operator+( const Poly & operant ) const;
  /**
   * \brief Compute the sum of this polynomial and a term.
   * \param[in] operant The term.
   * \return The sum of this polynomial and the term.
   */
  Poly operator+( const Term & operant ) const;
  /**
   * \brief Compute the difference of this and another polynomial.
   * \param[in] operant The other polynomial.
   * \return The difference between this and another polynomial.
   */
  Poly operator-( const Poly & operant ) const;
  /**
   * \brief Compute the different of this polynomial and a term.
   * \param[in] operant The term.
   * \return The difference between this polynomial and the term.
   */
  Poly operator-( const Term & operant ) const;
  /**
   * \brief Compute the product of this and another polynomial.
   * \param[in] operant The other polynomial.
   * \return The product of this and another polynomial.
   */
  Poly operator*( const Poly & operant ) const;
  /**
   * \brief Compute the product of this polynomial and a term.
   * \param[in] operant The term.
   * \return The product of this polynomial and the term.
   */
  Poly operator*( const Term & operant ) const;
  
  //in-place operations
  
  /**
   * \brief Extract a version of this polynomial that has different order.
   *        In-place!
   * \param[in] newOrder The desired monomial order.
   * \return A reference to this polynomial (updated).
   */
  Poly & differentOrderVersionInPlace( Monomial::Order newOrder );
  /**
   * \brief Generate a lower-degree approximation of this polynomial.
   *        In-place!
   * \param[in] maxDegree The maximum degree in the output polynomial.
   * \return A reference to this polynomial (updated).
   */
  Poly & lowerDegreeApproximationInPlace( unsigned int maxDegree );
  /**
   * \brief Get the negative of this polynomial.
   *        In-place!
   * \return A reference to this polynomial.
   */
  Poly & negationInPlace();
  /**
   * \brief Compute the sum of this and another polynomial.
   *        In-place!
   * \param[in] operant The other polynomial.
   * \return A reference to this polynomial.
   */
  Poly & operator+=( const Poly & operant );
  /**
   * \brief Compute the sum of this polynomial and a term.
   *        In-place!
   * \param[in] operant The term.
   * \return A reference to this polynomial.
   */
  Poly & operator+=( const Term & operant );
  /**
   * \brief Compute the difference of this and another polynomial.
   *        In-place!
   * \param[in] operant The other polynomial.
   * \return A reference to this polynomial.
   */
  Poly & operator-=( const Poly & operant );
  /**
   * \brief Compute the different of this polynomial and a term.
   *        In-place!
   * \param[in] operant The term.
   * \return A reference to this polynomial.
   */
  Poly & operator-=( const Term & operant );
  /**
   * \brief Compute the product of this and another polynomial.
   *        In-place!
   * \param[in] operant The other polynomial.
   * \return A reference to this polynomial.
   */
  Poly & operator*=( const Poly & operant );
  /**
   * \brief Compute the product of this polynomial and a term.
   *        In-place!
   * \param[in] operant The term.
   * \return A reference to this polynomial.
   */
  Poly & operator*=( const Term & operant );
  
  //possibly include the polynomial division here
  //Poly & operator/=( const Poly & operant );
  //Poly & operator/=( const Term & operant );
  
  // comparisons
  
  /**
   * \brief Check if this polynomial is equal to another one
   * \param[in] operant The comparison polynomial.
   * \return Is equal?
   */
  bool operator==( const Poly & operant ) const;
  /**
   * \brief Check if this polynomial is different from another one.
   * \param[in] operant The comparison polynomial.
   * \return Is different?
   */
  bool operator!=( const Poly & operant ) const;
  /**
   * \brief Check if this polynomial has identical form to another one. By
   *        identical form we understand here that it contains the same
   *        monomials in the same order.
   * \param[in] operant The comparison polynomial.
   * \return has identical form?
   */
  bool isSimilar( const Poly & operant ) const;

  //handy operations
  
  /**
   * \brief Set this polynomial to one.
   * \return A reference to this polynomial (updated).
   */
  Poly & setToOne();
  /**
   * \brief Set this polynomial to zero.
   * \return A reference to this polynomial (updated).
   */
  Poly & setToZero();

  //handy comparisons
  
  /**
   * \brief Is this polynomial zero?
   * \return Is zero?
   */
  bool isZero() const;
  /**
   * \brief Is this polynomial one?
   * \return Is one?
   */
  bool isOne() const;
  
private:
  /** The terms of this polynomial. */
  termsPtr _terms;
  /** The "sugar" of this polynomial. */
  unsigned int _sugar;
  
public:
  /** Useful named constructor idioms */
  
  // zero poly
  static Poly  zeroR( size_t dimensions, Monomial::Order order = Monomial::GREVLEX )
  { return Poly( Term::zeroR(dimensions,order)); };
  static Poly  zeroQ( size_t dimensions, Monomial::Order order = Monomial::GREVLEX )
  { return Poly( Term::zeroQ(dimensions,order)); };
  static Poly  zeroZ( size_t dimensions, Monomial::Order order = Monomial::GREVLEX )
  { return Poly( Term::zeroZ(dimensions,order)); };
  static Poly  zeroS( size_t dimensions, Monomial::Order order = Monomial::GREVLEX )
  { return Poly( Term::zeroS(dimensions,order)); };
  static Poly zeroSZ( size_t dimensions, Monomial::Order order = Monomial::GREVLEX )
  { return Poly(Term::zeroSZ(dimensions,order)); };
  
  // one poly
  static Poly  oneR( size_t dimensions, Monomial::Order order = Monomial::GREVLEX )
  { return Poly( Term::oneR(dimensions,order)); };
  static Poly  oneQ( size_t dimensions, Monomial::Order order = Monomial::GREVLEX )
  { return Poly( Term::oneQ(dimensions,order)); };
  static Poly  oneZ( size_t dimensions, Monomial::Order order = Monomial::GREVLEX )
  { return Poly( Term::oneZ(dimensions,order)); };
  static Poly  oneS( size_t dimensions, Monomial::Order order = Monomial::GREVLEX )
  { return Poly( Term::oneS(dimensions,order)); };
  static Poly oneSZ( size_t dimensions, Monomial::Order order = Monomial::GREVLEX )
  { return Poly(Term::oneSZ(dimensions,order)); };
  
  // constant poly
  static Poly  constR( int constant, size_t dimensions, Monomial::Order order = Monomial::GREVLEX )
  { return Poly( Term::constR(constant,dimensions,order)); };
  static Poly  constQ( int constant, size_t dimensions, Monomial::Order order = Monomial::GREVLEX )
  { return Poly( Term::constQ(constant,dimensions,order)); };
  static Poly  constZ( int constant, size_t dimensions, Monomial::Order order = Monomial::GREVLEX )
  { return Poly( Term::constZ(constant,dimensions,order)); };
  static Poly  constS( int constant, size_t dimensions, Monomial::Order order = Monomial::GREVLEX )
  { return Poly( Term::constS(constant,dimensions,order)); };
  static Poly constSZ( int constant, size_t dimensions, Monomial::Order order = Monomial::GREVLEX )
  { return Poly(Term::constSZ(constant,dimensions,order)); };
  
  // unknowns
  static Poly  uR( size_t index, size_t dimensions, Monomial::Order order = Monomial::GREVLEX )
  { return Poly( Term::uR(index,dimensions,order)); };
  static Poly  uQ( size_t index, size_t dimensions, Monomial::Order order = Monomial::GREVLEX )
  { return Poly( Term::uQ(index,dimensions,order)); };
  static Poly  uZ( size_t index, size_t dimensions, Monomial::Order order = Monomial::GREVLEX )
  { return Poly( Term::uZ(index,dimensions,order)); };
  static Poly  uS( size_t index, size_t dimensions, Monomial::Order order = Monomial::GREVLEX )
  { return Poly( Term::uS(index,dimensions,order)); };
  static Poly uSZ( size_t index, size_t dimensions, Monomial::Order order = Monomial::GREVLEX )
  { return Poly(Term::uSZ(index,dimensions,order)); };
  
  // unknowns multiplied by constant
  static Poly  cuR( int constant, size_t index, size_t dimensions, Monomial::Order order = Monomial::GREVLEX )
  { return Poly( Term::cuR(constant,index,dimensions,order)); };
  static Poly  cuQ( int constant, size_t index, size_t dimensions, Monomial::Order order = Monomial::GREVLEX )
  { return Poly( Term::cuQ(constant,index,dimensions,order)); };
  static Poly  cuZ( int constant, size_t index, size_t dimensions, Monomial::Order order = Monomial::GREVLEX )
  { return Poly( Term::cuZ(constant,index,dimensions,order)); };
  static Poly  cuS( int constant, size_t index, size_t dimensions, Monomial::Order order = Monomial::GREVLEX )
  { return Poly( Term::cuS(constant,index,dimensions,order)); };
  static Poly cuSZ( int constant, size_t index, size_t dimensions, Monomial::Order order = Monomial::GREVLEX )
  { return Poly(Term::cuSZ(constant,index,dimensions,order)); };
  
  // rand terms
  static Poly randR( size_t dimensions, Monomial::Order order = Monomial::GREVLEX )
  { return Poly( Term::randR(dimensions,order)); };
  static Poly randQ( size_t dimensions, Monomial::Order order = Monomial::GREVLEX )
  { return Poly( Term::randQ(dimensions,order)); };
  static Poly randZ( size_t dimensions, Monomial::Order order = Monomial::GREVLEX )
  { return Poly( Term::randZ(dimensions,order)); };
  
  // symbolic terms
  static Poly S( const std::string & name, size_t dimensions, Monomial::Order order = Monomial::GREVLEX )
  { return Poly(Term::S(name,dimensions,order)); };
  static Poly SrandZ( const std::string & name, size_t dimensions, Monomial::Order order = Monomial::GREVLEX )
  { return Poly(Term::SrandZ(name,dimensions,order)); };
};

}
}

#endif /* POLYJAM_CORE_POLY_HPP_ */
