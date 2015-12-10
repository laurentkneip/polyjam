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

/**
 * \file Term.hpp
 * \brief Definition of a term of a polynomial, including coefficient and
 *        monomial.
 */

#ifndef POLYJAM_CORE_TERM_HPP_
#define POLYJAM_CORE_TERM_HPP_

#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <vector>
#include <polyjam/core/Monomial.hpp>
#include <polyjam/core/Coefficient.hpp>

// todo: The printing might benefit from a function to check if negative

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
 * Term defines a polynomial term with coefficient and monomial. Can have a
 * multiple field representation.
 *
 * Careful: contructing and copying here is lazy. Use clone/copy if more is required.
 */
class Term
{
public:
  
  // constructors
  
  /**
   * \brief Constructor for a normal term.
   * \param[in] coefficient The coefficient.
   * \param[in] monomial The monomial.
   */
  Term( const Coefficient & coefficient, const Monomial & monomial );
  /**
   * \brief Constructor for a dual-field resentation term.
   * \param[in] coeff1 The first coefficient.
   * \param[in] coeff2 The second coefficient.
   * \param[in] monomial The monomial.
   */
  Term(
      const Coefficient & coeff1,
      const Coefficient & coeff2,
      const Monomial & monomial );
  /**
   * \brief Constructor for an arbitrary coefficient multiplicity.
   * \param[in] coefficients A vector with all coefficients.
   * \param[in] monomial The monomial.
   */
  Term(
      const std::vector<Coefficient> & coefficients,
      const Monomial & monomial );
  
  // Destructor
	
  virtual ~Term();

  // output
  /**
   * \brief Create a clone of this term (deep).
   * \param[in] full Maintain multiplicity?
   * \return A clone with the desired coefficients.
   */
  Term clone( bool full = true ) const;
  /**
   * \brief Copy the content from another term (deep).
   * \param[in] copy The original.
   */
  void copy( const Term & copy );
  /**
   * \brief Generate a similar term that is zero.
   * \param[in] full Maintain multiplicity?
   * \return A similar term that is zero.
   */
  Term zero( bool full = true ) const;
  /**
   * \brief Generate a similar term that is one.
   * \param[in] full Maintain multiplicity?
   * \param[in] index The index of the coefficient to add.
   * \return A similar term that is one.
   */
  Term one( bool full = true ) const;
  /**
   * \brief Is this a term with multiple coefficients?
   * \return Is multiple?
   */
  bool isMultiple() const;
  /**
   * \brief Set dominant coefficient.
   * \param[in] index The index of the dominant coefficient.
   */
  void setDominant( int index );
  /**
   * \brief Print on console (dominant one).
   */
  void print() const;
  /**
   * \brief Get a string that represents the term (dominant one).
   * \param[in] c_version Get the string in C++-syntax using pow(...,...)?
   * \return The string representing the monomial.
   */
  std::string getString( bool c_version = true ) const;
  /**
   * \brief Get a reference to the monomial.
   * \return A reference to the monomial.
   */
  const Monomial & monomial() const;
  /**
   * \brief Get a reference to the dominant coefficient.
   * \return A reference to the dominant coefficient.
   */
  const Coefficient & coefficient() const;
  /**
   * \brief Reset the monomial order
   * \param[in] newOrder The new monomial order for this term.
   */
  void setOrder( Monomial::Order newOrder );
  
  // operations
  
  /**
   * \brief Compute the negation of this term.
   * \return The negative version of this term.
   */
  Term negation() const;
  /**
   * \brief Compute the sum of this and another term.
   * \param[in] The operant.
   * \return The sum.
   */
  Term operator+( const Term & operant ) const;
  /**
   * \brief Compute the difference of this and another term.
   * \param[in] The operant.
   * \return The difference.
   */
  Term operator-( const Term & operant ) const;
  /**
   * \brief Compute the product of this and another term.
   * \param[in] The operant.
   * \return The product.
   */
  Term operator*( const Term & operant ) const;
  /**
   * \brief Compute the quotient of this and another term.
   * \param[in] The operant.
   * \return The quotient.
   */
  Term operator/( const Term & operant ) const;
  
  // in-place operations
  
  /**
   * \brief Make this term its negative.
   * \return A reference to this term.
   */
  Term & negationInPlace();
  /**
   * \brief Add a term to this one.
   * \param[in] The operant.
   * \return A reference to this term.
   */
  Term & operator+=( const Term & operant );
  /**
   * \brief Subtract a term from this one.
   * \param[in] The operant.
   * \return A reference to this term.
   */
  Term & operator-=( const Term & operant );
  /**
   * \brief Multiply this term by another one.
   * \param[in] The operant.
   * \return A reference to this term.
   */
  Term & operator*=( const Term & operant );
  /**
   * \brief Divide this term by another one.
   * \param[in] The operant.
   * \return A reference to this term.
   */
  Term & operator/=( const Term & operant );
  
  // comparisons
  
  /**
   * \brief Check if this term equals to another term.
   * \param[in] operant The comparison term.
   * \return Is equal?
   */
  bool operator==( const Term & operant ) const;
  /**
   * \brief Check if this term is different from another term.
   * \param[in] operant The comparison term.
   * \return Is different?
   */
  bool operator!=( const Term & operant ) const;
  /**
   * \brief Check if is smaller (forwarded to monomial check/used for
   *        polynomial ordering).
   * \param[in] operant The comparison term.
   * \return Is smaller?
   */
  bool operator>( const Term & operant ) const;
  /**
   * \brief Check if this term is similar to another one. By similar we
   *        understand here that the fields, dimensions, and order are same.
   * \return Are these terms similar?
   */
  bool isSimilar( const Term & operant ) const;
  
  // handy operations
  
  /**
   * \brief Set this term to one.
   */
  Term & setToOne();
  /**
   * \brief Set this term to zero.
   */
  Term & setToZero();
  
  // handy comparisons
  
  /**
   * \brief Check if this term equals to one.
   * \return Is one?
   */
  bool isOne() const;
  /**
   * \brief Check if this term equals to zero.
   * \return Is zero?
   */
  bool isZero() const;
	
private:

  /** The coefficients of this term */
  mutable std::vector<Coefficient> _coefficients;
  /** The currently dominant coefficient */
  mutable int _dominantIndex;
  /** The monomial of this term */
  Monomial _monomial;

  /**
   * \brief Check if this and another term have opposing multiplicity.
   * \param[in] operant The comparison Term.
   * \return Different multiplicity?
   */
  bool isWrong( const Term & operant ) const;

public:
  /** Useful named constructor idioms */
  
  // zero terms of different fields
  static Term  zeroR( size_t dimensions, Monomial::Order order = Monomial::GREVLEX )
  { return Term( Coefficient(fields::Field::R),   Monomial(dimensions,order) ); };
  static Term  zeroQ( size_t dimensions, Monomial::Order order = Monomial::GREVLEX )
  { return Term( Coefficient(fields::Field::Q),   Monomial(dimensions,order) ); };
  static Term  zeroZ( size_t dimensions, Monomial::Order order = Monomial::GREVLEX )
  { return Term( Coefficient(fields::Field::Zp),  Monomial(dimensions,order) ); };
  static Term  zeroS( size_t dimensions, Monomial::Order order = Monomial::GREVLEX )
  { return Term( Coefficient(fields::Field::Sym), Monomial(dimensions,order) ); };
  static Term zeroSZ( size_t dimensions, Monomial::Order order = Monomial::GREVLEX )
  { return Term( Coefficient(fields::Field::Sym), Coefficient(fields::Field::Zp), Monomial(dimensions,order) ); };
  
  // one terms of different fields
  static Term  oneR( size_t dimensions, Monomial::Order order = Monomial::GREVLEX )
  { return Term( Coefficient(1,fields::Field::R),   Monomial(dimensions,order) ); };
  static Term  oneQ( size_t dimensions, Monomial::Order order = Monomial::GREVLEX )
  { return Term( Coefficient(1,fields::Field::Q),   Monomial(dimensions,order) ); };
  static Term  oneZ( size_t dimensions, Monomial::Order order = Monomial::GREVLEX )
  { return Term( Coefficient(1,fields::Field::Zp),  Monomial(dimensions,order) ); };
  static Term  oneS( size_t dimensions, Monomial::Order order = Monomial::GREVLEX )
  { return Term( Coefficient(1,fields::Field::Sym), Monomial(dimensions,order) ); };
  static Term oneSZ( size_t dimensions, Monomial::Order order = Monomial::GREVLEX )
  { return Term( Coefficient(1,fields::Field::Sym), Coefficient(1,fields::Field::Zp), Monomial(dimensions,order) ); };
  
  // const terms of different fields
  static Term  constR( int constant, size_t dimensions, Monomial::Order order = Monomial::GREVLEX )
  { return Term( Coefficient(constant,fields::Field::R),   Monomial(dimensions,order) ); };
  static Term  constQ( int constant, size_t dimensions, Monomial::Order order = Monomial::GREVLEX )
  { return Term( Coefficient(constant,fields::Field::Q),   Monomial(dimensions,order) ); };
  static Term  constZ( int constant, size_t dimensions, Monomial::Order order = Monomial::GREVLEX )
  { return Term( Coefficient(constant,fields::Field::Zp),  Monomial(dimensions,order) ); };
  static Term  constS( int constant, size_t dimensions, Monomial::Order order = Monomial::GREVLEX )
  { return Term( Coefficient(constant,fields::Field::Sym), Monomial(dimensions,order) ); };
  static Term constSZ( int constant, size_t dimensions, Monomial::Order order = Monomial::GREVLEX )
  { return Term( Coefficient(constant,fields::Field::Sym), Coefficient(constant,fields::Field::Zp), Monomial(dimensions,order) ); };
  
  // unknowns of different fields
  static Term  uR( size_t index, size_t dimensions, Monomial::Order order = Monomial::GREVLEX )
  { return Term( Coefficient(1,fields::Field::R),   Monomial(dimensions,index,order) ); };
  static Term  uQ( size_t index, size_t dimensions, Monomial::Order order = Monomial::GREVLEX )
  { return Term( Coefficient(1,fields::Field::Q),   Monomial(dimensions,index,order) ); };
  static Term  uZ( size_t index, size_t dimensions, Monomial::Order order = Monomial::GREVLEX )
  { return Term( Coefficient(1,fields::Field::Zp),  Monomial(dimensions,index,order) ); };
  static Term  uS( size_t index, size_t dimensions, Monomial::Order order = Monomial::GREVLEX )
  { return Term( Coefficient(1,fields::Field::Sym), Monomial(dimensions,index,order) ); };
  static Term uSZ( size_t index, size_t dimensions, Monomial::Order order = Monomial::GREVLEX )
  { return Term( Coefficient(1,fields::Field::Sym), Coefficient(1,fields::Field::Zp), Monomial(dimensions,index,order) ); };
  
  // unknowns of different fields (multiplied by constant)
  static Term  cuR( int constant, size_t index, size_t dimensions, Monomial::Order order = Monomial::GREVLEX )
  { return Term( Coefficient(constant,fields::Field::R),   Monomial(dimensions,index,order) ); };
  static Term  cuQ( int constant, size_t index, size_t dimensions, Monomial::Order order = Monomial::GREVLEX )
  { return Term( Coefficient(constant,fields::Field::Q),   Monomial(dimensions,index,order) ); };
  static Term  cuZ( int constant, size_t index, size_t dimensions, Monomial::Order order = Monomial::GREVLEX )
  { return Term( Coefficient(constant,fields::Field::Zp),  Monomial(dimensions,index,order) ); };
  static Term  cuS( int constant, size_t index, size_t dimensions, Monomial::Order order = Monomial::GREVLEX )
  { return Term( Coefficient(constant,fields::Field::Sym), Monomial(dimensions,index,order) ); };
  static Term cuSZ( int constant, size_t index, size_t dimensions, Monomial::Order order = Monomial::GREVLEX )
  { return Term( Coefficient(constant,fields::Field::Sym), Coefficient(constant,fields::Field::Zp), Monomial(dimensions,index,order) ); };
  
  // rand terms of different fields
  static Term randR( size_t dimensions, Monomial::Order order = Monomial::GREVLEX )
  { return Term( Coefficient(fields::Field::R, true), Monomial(dimensions,order) ); };
  static Term randQ( size_t dimensions, Monomial::Order order = Monomial::GREVLEX )
  { return Term( Coefficient(fields::Field::Q, true), Monomial(dimensions,order) ); };
  static Term randZ( size_t dimensions, Monomial::Order order = Monomial::GREVLEX )
  { return Term( Coefficient(fields::Field::Zp,true), Monomial(dimensions,order) ); };
  
  // symbolic terms
  static Term S( const std::string & name, size_t dimensions, Monomial::Order order = Monomial::GREVLEX )
  { return Term( Coefficient(name), Monomial(dimensions,order) ); };
  static Term SrandZ( const std::string & name, size_t dimensions, Monomial::Order order = Monomial::GREVLEX )
  { return Term( Coefficient(name), Coefficient( fields::Field::Zp, true ), Monomial(dimensions,order) ); };
};

}
}

#endif /* POLYJAM_CORE_TERM_HPP_ */
