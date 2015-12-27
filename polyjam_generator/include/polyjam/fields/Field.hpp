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
 * \file Field.hpp
 * \brief Abstract class for a member of an algebraic field.
 */

#ifndef POLYJAM_FIELDS_FIELD_HPP_
#define POLYJAM_FIELDS_FIELD_HPP_

#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <iostream>


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
 * The class Field defines a generic type for an algebraic field,
 * including it's operations addition, subtraction, multiplication, division,
 * ana a member representation in the derivatives.
 */
class Field
{
public:
  
  /** The different kinds of fields that are available */
  enum Kind
  {
    R,
    Q,
    Zp,
    Sym
  };

  /**
   * \brief Default constructor.
   * \param[in] kind The kind of the field.
   */
  Field( Kind kind ) : _kind(kind) {};
  /**
   * \brief Destructor.
   */
  virtual ~Field() {};

  // output

  /**
   * \brief Generic output of the member.
   */
  void print() const
  {
    std::cout << getString(false);
  };
  /**
   * \brief Get a string of the member.
   * \return A string describing the member.
   */
  virtual std::string getString( bool c_version = true ) const = 0;
  /**
   * \brief Get the kind of this field.
   * \return The kind of this field
   */
  Kind kind() const { return _kind; };

  // Get constants from this field
  
  /**
   * \brief Create a field-member that is zero.
   * \return A reference to a zero field-member.
   */
  virtual Field* zero() const = 0;
  /**
   * \brief Create a field-member that is one.
   * \return A reference to a one field-member.
   */
  virtual Field* one() const = 0;

  // standard operations (in-place!)

  /**
   * \brief Compute the negation of this number.
   * \return The negation of this number.
   */
  virtual void negation() = 0;
  /**
   * \brief Compute the inversion of this number.
   * \return The inversion of this number.
   */
  virtual void inversion() = 0;
  
  /**
   * \brief Compute the sum of two numbers.
   * \param[in] operant The operant.
   * \return The sum of the two numbers.
   */
  virtual void add( const Field* operant ) = 0;
  /**
   * \brief Compute the signed difference of two numbers.
   * \param[in] operant The operant.
   * \return The signed difference of the two numbers.
   */
  virtual void subtract( const Field* operant ) = 0;
  /**
   * \brief Compute the multiplication of two numbers.
   * \param[in] operant The operant.
   * \return The multiplication of the two numbers.
   */
  virtual void multiply( const Field* operant ) = 0;
  /**
   * \brief Compute the division of two numbers.
   * \param[in] operant The operant.
   * \return The division of the two numbers.
   */
  virtual void divide( const Field* operant ) = 0;
  
  // comparison operations

  /**
   * \brief Is equal?
   * \param[in] operant The operant.
   * \return Is equal?
   */   
  virtual bool isEql( const Field* operant ) const = 0;
  /**
   * \brief Compare two field members. Is separated from isEql, because isEql
   *        might make sense for all fields, but compare not.
   * \param[in] operant The operant.
   * \return The comparison result (1 = bigger, -1 = smaller, 0 = equal)
   */
  virtual int compare( const Field* operant ) const = 0;
  
  // helper
  
  /**
   * \brief Check if operant is from the wrong field.
   * \param[in] operant The operant we want to check.
   * \return Is wrong field?
   */
  virtual bool isWrong( const Field* operant ) const
  {
    if( operant->kind() != _kind )
    {
      std::cout << "Error: Field not compatible" << std::endl;
      return true;
    }
    return false;
  }
  
private:

  /** The kind of this field */
  Kind _kind;
};

}
}

#endif /* POLYJAM_FIELDS_FIELD_HPP_ */
