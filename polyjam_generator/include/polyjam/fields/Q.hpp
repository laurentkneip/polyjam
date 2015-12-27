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
 * \file Q.hpp
 * \brief Member of the algebraic field Q.
 */

#ifndef POLYJAM_FIELDS_Q_HPP_
#define POLYJAM_FIELDS_Q_HPP_

#include <polyjam/fields/Field.hpp>
#include <stdint.h>

/**
 * \brief The namespace of this library.
 */
namespace polyjam
{

/**
 * \brief The namespace for the fields.
 */
namespace fields
{

/**
 * The class Q defines a member of the algebraic field Q,
 * including addition, subtraction, multiplication, division.
 */
class Q : public Field
{
public:

  /**
   * \brief Constructor for a zero or random number.
   * \param[in] random Set this to a random number?
   */
  Q( bool random = false );
  /**
   * \brief Constructor.
   * \param[in] numerator Initial numerator value.
   * \param[in] denominator Initial denominator value.
   */   
  Q( int64_t numerator, uint64_t denominator );
  /**
   * \brief Constructor.
   * \param[in] value The initial value.
   */
  Q( int value );
  /**
   * \brief Copy-Constructor (deep).
   * \param[in] copy Original.
   */
  Q( const Q & copy );
  /**
   * \brief Copy-Constructor (deep).
   * \param[in] copy Original.
   */
  Q( const Field* copy );
  /**
   * Destructor.
   */
  virtual ~Q();

  // output
  
  /** See base-class documentation */
  std::string getString( bool c_version = true ) const;

  // Get constants from this field
  
  /** See base-class documentation */
  Field* zero() const;
  /** See base-class documentation */
  Field* one() const;

  // standard self-operations

  /** See base-class documentation */
  void negation();
  /** See base-class documentation */
  void inversion();

  // standard operations
    
  /** See base-class documentation */
  void add( const Field* operant );
  /** See base-class documentation */
  void subtract( const Field* operant );
  /** See base-class documentation */
  void multiply( const Field* operant );
  /** See base-class documentation */
  void divide( const Field* operant );
  
  // comparison operations

  /** See base-class documentation */
  bool isEql( const Field* operant ) const;
  /** See base-class documentation */
  int compare( const Field* operant ) const;
  
private:

  /** The numerator */  
  int64_t _numerator;
  /** The denominator */
  uint64_t _denominator;
  
  /**
   * \brief Type converting absolue value.
   * \param[in] value The original value in int64_t
   * \return The return value in uint64_t
   */
  uint64_t specialAbs( int64_t value ) const;
  
  /**
   * \brief Compute the greatest common divisor of two numbers.
   * \param[in] a The first number.
   * \param[in] b The second number.
   * \return The greatest common divisor.
   */
  uint64_t GCD(uint64_t a, uint64_t b) const;
  
  /**
   * \brief Simplify this fraction.
   * \return The simplified fraction.
   */
  void clean();
};

}
}

#endif /* POLYJAM_FIELDS_Q_HPP_ */
