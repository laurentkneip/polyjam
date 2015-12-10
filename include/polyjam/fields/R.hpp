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
 * \file R.hpp
 * \brief Member of the algebraic field R.
 */

#ifndef POLYJAM_FIELDS_R_HPP_
#define POLYJAM_FIELDS_R_HPP_

#include <polyjam/fields/Field.hpp>

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
 * The class R defines a member of the algebraic field R,
 * including addition, subtraction, multiplication, division.
 */
class R : public Field
{
public:

  /**
   * \brief Constructor for a zero or random number.
   * \param[in] random Set this to a random number?
   */
  R( bool random = false );
  /**
   * \brief Constructor.
   * \param[in] value The initial value.
   */
  R( double value );
  /**
   * \brief Copy-Constructor (deep).
   * \param[in] copy Original.
   */
  R( const R & copy );
  /**
   * \brief Copy-Constructor (deep).
   * \param[in] copy Original.
   */
  R( const Field* copy );
  /**
   * Destructor.
   */
  virtual ~R();

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
  /** The number */
  double _value;
};

}
}

#endif /* POLYJAM_FIELDS_R_HPP_ */
