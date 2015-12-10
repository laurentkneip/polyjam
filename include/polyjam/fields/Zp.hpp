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
