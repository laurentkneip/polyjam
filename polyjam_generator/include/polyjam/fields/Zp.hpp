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
 * \file Zp.hpp
 * \brief Member of the algebraic field Zp (Prime field).
 */

#ifndef POLYJAM_FIELDS_ZP_HPP_
#define POLYJAM_FIELDS_ZP_HPP_

#include <polyjam/fields/Field.hpp>

#define DEFAULT_CHARACTERISTIC 30097

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
 * The class Zp defines a member of the algebraic prime field Zp,
 * including addition, subtraction, multiplication, division.
 */
class Zp : public Field
{
public:

  /**
   * \brief Constructor for zero or random member.
   * \param[in] random Sets this number to a random prime field member.
   * \param[in] prime The characteristic of the prime field.
   */
  Zp( bool random = false, unsigned int characteristic = DEFAULT_CHARACTERISTIC );
  /**
   * \brief Constructor.
   * \param[in] value The initial value of the prime-field member.
   * \param[in] prime The characteristic of the prime field.
   */
  Zp( int value, unsigned int characteristic = DEFAULT_CHARACTERISTIC );
  /**
   * \brief Copy-Constructor (deep).
   * \param[in] copy Original.
   */
  Zp( const Zp & copy );
  /**
   * \brief Copy-Constructor (deep).
   * \param[in] copy Original.
   */
  Zp( const Field* copy );
  /**
   * \brief Destructor.
   */
  virtual ~Zp();
  
  // output
  
  /** See base-class documentation */
  std::string getString( bool c_version = true ) const;

  /**
   * \brief Get the characteristic of this field.
   * \return The characteristic.
   */
  unsigned int characteristic() const;

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
  /** The value of this prime field member */
  unsigned int _value;
  /** The characteristic of this prime field */
  unsigned int _characteristic;

  /**
   * \brief Compute the multiplicative inverse of a prime field member.
   * \param[in] value The value of which we want to compute the multiplicative
   *                  inverse.
   * \return The multiplicative inverse.
   */
  unsigned int getMultiplicativeInverse( unsigned int value ) const;
  /**
   * \brief Take the "modulo-prime" value of a value.
   * \param[in] value The value that we want to run the modulo operation on.
   * \return The value "modulo-prime".
   */
  unsigned int moduloPrime( int value ) const;
  /**
   * \brief Check if the operant is from the wrong field
   * \param[in] operant The operant we want to check.
   * \return Is wrong field?
   */
  bool isWrong( const Field* operant ) const;
};

}
}

#endif /* POLYJAM_FIELDS_ZP_HPP_ */
