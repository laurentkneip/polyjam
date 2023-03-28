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

/**
 * \file Coefficient.hpp
 * \brief Class for a Coefficient (e.g. of a term of a polynomial).
 */

#ifndef POLYJAM_CORE_COEFFICIENT_HPP_
#define POLYJAM_CORE_COEFFICIENT_HPP_

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string>
#include <memory>

#include <polyjam/fields/Field.hpp>

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
 * This is the class for a coefficient. It holds an abstract member of a
 * field, and allows to do various operations on the coefficients. It implements
 * the factory pattern to create coefficients from different fields, and
 * therefore prevents the use of templates throughout the rest of the library.
 *
 * Careful: assignment and copy is lazy by default!
 */
class Coefficient
{
public:

  /** A pointer to the actual field-object underneath */
  typedef std::shared_ptr<fields::Field> FieldPtr;
  
  // Constructors

  /**
   * \brief Creates a zero member of a certain field.
   * \param[in] kind The kind of Field.
   * \param[in] random Set the coefficient to random (only if not Symbolic).
   */
  Coefficient( fields::Field::Kind kind, bool random = false );
  /**
   * \brief Constructor for a member of R.
   * \param[in] value The value of the coefficient.
   */
  Coefficient( double value );
  /**
   * \brief Constructor for a member of Q.
   * \param[in] numerator The numberator of the coefficient.
   * \param[in] denominator The denominator of the coefficient.
   */
  Coefficient( int numerator, unsigned int denominator );
  /**
   * \brief Constructor for a constant in any field. Interesting for
   *        Zp and Symbolic.
   * \param[in] constant The value of the coefficient.
   * \param[in] kind The kind of field we want this constant to be on?
   */
  Coefficient( int constant, fields::Field::Kind kind );
  /**
   * \brief Constructor for a symbolic coefficient.
   * \param[in] name The name of the coefficient.
   */
  Coefficient( const std::string & name );
  /**
   * \brief Constructor for a symbolic coefficient.
   * \param[in] name The name of the coefficient.
   */
  Coefficient( const char* name );
  
  // Destructor
  
  /**
   * Destructor.
   */
  virtual ~Coefficient();
  
  // deep copy stuff
  
  /**
   * \brief Create a new object that is equal to the present one.
   * \return The clone.
   */
  Coefficient clone() const;
  
  /**
   * \brief Copy the values from another coefficient.
   * \param[in] coefficient The original.
   */
  void copy( const Coefficient & coefficient );
  
  
  // output
  
  /**
   * \brief Print this coefficient on the console.
   */
  void print() const;
  /**
   * \brief Get a C-syntax style string of this coefficient.
   * \return A C-syntax style string of this coefficient.
   */
  std::string getString( bool c_version = true ) const;
  /**
   * \brief Get a C-syntax style string of this coefficient.
   * \return A C-syntax style string of this coefficient.
   * \note: this one is however calling the special one for a symbolic member
   * (without integer coefficients)
   */
  std::string getStringSpecial( bool c_version = true ) const;
  /**
   * \brief Retrieve the kind of this field.
   * \return The kind of this field.
   */
  fields::Field::Kind kind() const;
  /**
   * \brief Generate a similar coefficient that is zero.
   * \return A similar coefficient that is zero.
   */
  Coefficient zero() const;
  /**
   * \brief Generate a similar coefficient that is one.
   * \return A similar coefficient that is one.
   */
  Coefficient one() const;
  /**
   * \brief Get the characteristic (of course only works for Zp)
   * \return Characteristic.
   */
  unsigned int characteristic() const;

  // standard operations
  
  /**
   * \brief Compute the negation of this coefficient.
   * \return The negative version of this coefficient.
   */
  Coefficient negation() const;
  /**
   * \brief Compute the inverse of this coefficient.
   * \return The inverse of this coefficient.
   */
  Coefficient inversion() const;
  /**
   * \brief Compute the sum of this and another coefficient.
   * \param[in] The operant.
   * \return The sum.
   */
  Coefficient operator+( const Coefficient & operant ) const;
  /**
   * \brief Compute the difference of this and another coefficient.
   * \param[in] The operant.
   * \return The difference.
   */
  Coefficient operator-( const Coefficient & operant ) const;
  /**
   * \brief Compute the product of this and another coefficient.
   * \param[in] The operant.
   * \return The product.
   */
  Coefficient operator*( const Coefficient & operant ) const;
  /**
   * \brief Compute the quotient of this and another coefficient.
   * \param[in] The operant.
   * \return The quotient.
   */
  Coefficient operator/( const Coefficient & operant ) const;
  
  // in-place operations
  
  /**
   * \brief Negate this coefficient.
   */
  Coefficient & negationInPlace();
  /**
   * \brief Invert this coefficient.
   */
  Coefficient & inversionInPlace();
  /**
   * \brief Add a coefficient to this one.
   * \param[in] The operant.
   * \return A reference to this coefficient.
   */
  Coefficient & operator+=( const Coefficient & operant );
  /**
   * \brief Subtract a coefficient from this one.
   * \param[in] The operant.
   * \return A reference to this coefficient.
   */
  Coefficient & operator-=( const Coefficient & operant );
  /**
   * \brief Multiply this coefficient by another one.
   * \param[in] The operant.
   * \return A reference to this coefficient.
   */
  Coefficient & operator*=( const Coefficient & operant );
  /**
   * \brief Divide this coefficient by another one.
   * \param[in] The operant.
   * \return A reference to this coefficient.
   */
  Coefficient & operator/=( const Coefficient & operant );
  
  // comparisons
  
  /**
   * \brief Is equal?
   * \param[in] operant The operant.
   * \return Is equal?
   */
  bool operator==( const Coefficient & operant ) const;
  /**
   * \brief Is different?
   * \param[in] operant The operant.
   * \return Is different?
   */
  bool operator!=( const Coefficient & operant ) const;
  /**
   * \brief Is smaller or equal?
   * \param[in] operant The operant.
   * \return Is smaller or equal?
   */
  bool operator<=( const Coefficient & operant ) const;
  /**
   * \brief Is bigger or equal?
   * \param[in] operant The operant.
   * \return Is bigger or equal?
   */
  bool operator>=( const Coefficient & operant ) const;
  /**
   * \brief Is smaller?
   * \param[in] operant The operant.
   * \return Is smaller?
   */
  bool operator<( const Coefficient & operant ) const;
  /**
   * \brief Is bigger?
   * \param[in] operant The operant.
   * \return Is bigger?
   */
  bool operator>( const Coefficient & operant ) const;
  
  // handy operations
  
  /**
   * \brief Set this coefficient to zero.
   */
  Coefficient & setToZero();
  /**
   * \brief Set this coefficient to one.
   */
  Coefficient & setToOne();
  
  // handy comparisons
  
  /**
   * \brief Is this coefficient zero?
   * \return Is zero?
   */
  bool isZero() const;
  /**
   * \brief Is this coefficient one?
   * \return Is one?
   */
  bool isOne() const;
  
private:

  /** The value/field of this coefficient */
  FieldPtr _field;
  
  /**
   * \brief Another constructor (lazy copy, only used internally!)
   * \param[in] field The original.
   */
  Coefficient( fields::Field * field );
};

}
}

#endif /* POLYJAM_CORE_COEFFICIENT_HPP_ */
