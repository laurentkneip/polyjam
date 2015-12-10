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
