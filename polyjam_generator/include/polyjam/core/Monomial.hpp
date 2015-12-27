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
 * \file Monomial.hpp
 * \brief Definition of a monomial of any ordering, degree and dimensions.
 */

#ifndef POLYJAM_CORE_MONOMIAL_HPP_
#define POLYJAM_CORE_MONOMIAL_HPP_

#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <string>


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
 * The monomial defines a monomial with a bunch of handly operations on it.
 */
class Monomial
{
public:

  /**
   * the orders we can have, will be used in polynomial for ordering the
   * terms
   */
  enum Order
  {
    LEX,
    REVLEX,
    GRLEX,
    GREVLEX
  };

  // constructors
  
  /**
   * \brief Zero constructor for a monomial of certain dimensions.
   * \param[in] dimensions The number of unknowns in the monomial.
   * \param[in] order The used monomial-order.
   */
  Monomial( size_t dimensions, Order order = Monomial::GREVLEX );
  /**
   * \brief Constructor for a monomial of certain dimensions.
   * \param[in] dimensions The number of unknowns in the monomial.
   * \param[in] exponents The exponents for each unknown.
   * \param[in] order The used monomial-order.
   */
  Monomial( size_t dimensions, const unsigned int * exponents, Order order = Monomial::GREVLEX );
  /**
   * \brief Constructor given the exponents (for each dimension, also 0).
   * \param[in] exponents The exponents for the unknowns. Needs to contain zeros
   *                      for the unknowns that are not contained.
   * \param[in] order The used monomial-order.
   */
  Monomial( const std::vector<unsigned int> & exponents, Order order = Monomial::GREVLEX );
  /**
   * \brief Constructor for an unknown.
   * \param[in] dimensions The number of unknowns in the monomial.
   * \param[in] indexOne The position of the unknown (starting at 1). The monomial
   *                     will be set to zero if this variable is 0.
   * \param[in] order The used monomial-order.
   */
  Monomial( size_t dimensions, size_t indexOne, Order order = Monomial::GREVLEX );

  // destructors
  
  /**
   * \brief Destuctor.
   */
  virtual ~Monomial();

  // output
  
  /**
   * \brief Print this monomial on the console.
   */
  void print() const;
  /**
   * \brief Get a string that represents the monomial.
   * \param[in] c_version Get the string in C++-syntax using pow(...,...)?
   * \return The string that describes the monomial.
   */
  std::string getString( bool c_version = true ) const;
  /**
   * \brief Get a string with only the exponents (useful for code generator).
   * \return A string with only the exponents.
   */
  std::string getAlpha() const;
  /**
   * \brief Access the exponents.
   * \return A reference to the exponents.
   */
  const std::vector<unsigned int> & exponents() const;
  /**
   * \brief Return a monomial that is similar to this one, but one.
   * \return A monomial that equals to one (all exponents zero).
   */
  Monomial one() const;
  /**
   * \brief Evaluate this monomial by replacing all unknowns with the
   *        provided values.
   * \param[in] values The replacement values for each unknown.
   */
  double eval( const std::vector<double> & values ) const;
  
  // properties
  
  /**
   * \brief The total degree of this monomial.
   * \return The total degree of this monomial.
   */
  unsigned int degree() const;
  /**
   * \brief The number of unknowns in this monomial.
   * \return The number of unknowns in this monomial.
   */
  unsigned int dimensions() const;
  /**
   * \brief The order used for comparing this monomial.
   * \return The order used for comparing this monomial.
   */
  Order order() const;
  /**
   * \brief Reset the monomial order. This function is used when resetting the
   *        order of polynomials.
   * \param[in] newOrder The new monomial order.
   */
  void setOrder( Order newOrder );
  
  // operations
  
  /**
   * \brief Compute the least common multiple of this and another monomial.
   * \param[in] operant The second monomial.
   * \return The least common multiple.
   */
  Monomial leastCommonMultiple( const Monomial & operant ) const;
  /**
   * \brief Multiply this monomial by another one.
   * \param[in] operant The multiplier.
   * \return The multiplied monomial.
   */
  Monomial operator*(const Monomial & operant) const;
  /**
   * \brief Divide this monomial by another one.
   * \param[in] operant The divider.
   * \return The divided monomial.
   */
  Monomial operator/(const Monomial & operant) const;
  /**
   * \brief Multiply this monomial by another one (in place).
   * \param[in] operant The multiplier.
   * \return A reference to this monomial.
   */
  Monomial & operator*=(const Monomial & operant);
  /**
   * \brief Divide this monomial by another one (in place).
   * \param[in] operant The divider.
   * \return A reference to this monomial.
   */
  Monomial & operator/=(const Monomial & operant);
  
  // comparisons
  
  /**
   * \brief Check if is equal to another monomial.
   * \param[in] operant The comparison monomial.
   * \return Is equal?
   */
  bool operator==( const Monomial & operant ) const;
  /**
   * \brief Check if is different from another monomial.
   * \param[in] operant The comparison monomial.
   * \return Is different?
   */
  bool operator!=( const Monomial & operant ) const;
  /**
   * \brief Check if is bigger or equal than another monomial (default ordering).
   * \param[in] operant The comparison monomial.
   * \return Is bigger or equal?
   */
  bool operator>=( const Monomial & operant ) const;
  /**
   * \brief Check if is smaller or equal than another monomial (default ordering).
   * \param[in] operant The comparison monomial.
   * \return Is smaller or equal?
   */
  bool operator<=( const Monomial & operant ) const;
  /**
   * \brief Check if is bigger than another monomial (default ordering).
   * \param[in] operant The comparison monomial.
   * \return Is bigger?
   */
  bool operator>( const Monomial & operant ) const;
  /**
   * \brief Check if is smaller than another monomial (default ordering).
   * \param[in] operant The comparison monomial.
   * \return Is smaller?
   */
  bool operator<( const Monomial & operant ) const;
  /**
   * \brief Check if is bigger than another monomial.
   * \param[in] operant The comparison monomial.
   * \param[in] order The monomial order we are using.
   * \return Is bigger?
   */
  bool isBigger( const Monomial & operant, const Order & order ) const;
  /**
   * \brief Check if is smaller than another monomial.
   * \param[in] operant The comparison monomial.
   * \param[in] order The monomial order we are using.
   * \return Is smaller?
   */
  bool isSmaller( const Monomial & operant, const Order & order ) const;
  /**
   * \brief A general comparison, using the indicated order.
   * \param[in] operant The comparison monomial.
   * \param[in] order The monomial order we are using.
   * \return -1 if this one is smaller, 0 if equal, and 1 if bigger.
   */
  int comparison(const Monomial & operant, const Order & order ) const;
  /**
   * \brief Check if is lex-bigger than another monomial.
   * \param[in] operant The comparison monomial.
   * \return Is lex-bigger?
   */
  bool isLexBigger( const Monomial & operant ) const;
  /**
   * \brief Check if is lex-smaller than another monomial.
   * \param[in] operant The comparison monomial.
   * \return Is lex-smaller?
   */
  bool isLexSmaller( const Monomial & operant ) const;
  /**
   * \brief Check how this monomial lex-compares to another one.
   * \param[in] operant The comparison monomial.
   * \return -1 if this one is smaller, 0 if equal, and 1 if bigger.
   */
  int lexComparison( const Monomial & operant ) const;
  /**
   * \brief Check if is revlex-bigger than another monomial.
   * \param[in] operant The comparison monomial.
   * \return Is revlex-bigger?
   */
  bool isRevlexBigger( const Monomial & operant ) const;
  /**
   * \brief Check if is revlex-smaller than another monomial.
   * \param[in] operant The comparison monomial.
   * \return Is revlex-smaller?
   */
  bool isRevlexSmaller( const Monomial & operant ) const;
  /**
   * \brief Check how this monomial revlex-compares to another one.
   * \param[in] operant The comparison monomial.
   * \return -1 if this one is smaller, 0 if equal, and 1 if bigger.
   */
  int revlexComparison( const Monomial & operant ) const;
  /**
   * \brief Check if is grlex-bigger than another monomial.
   * \param[in] operant The comparison monomial.
   * \return Is grlex-bigger?
   */
  bool isGrlexBigger( const Monomial & operant ) const;
  /**
   * \brief Check if is grlex-smaller than another monomial.
   * \param[in] operant The comparison monomial.
   * \return Is grlex-smaller?
   */
  bool isGrlexSmaller( const Monomial & operant ) const;
  /**
   * \brief Check how this monomial grlex-compares to another one.
   * \param[in] operant The comparison monomial.
   * \return -1 if this one is smaller, 0 if equal, and 1 if bigger.
   */
  int grlexComparison( const Monomial & operant ) const;
  /**
   * \brief Check if is grevlex-bigger than another monomial.
   * \param[in] operant The comparison monomial.
   * \return Is grevlex-bigger?
   */
  bool isGrevlexBigger( const Monomial & operant ) const;
  /**
   * \brief Check if is grevlex-smaller than another monomial.
   * \param[in] operant The comparison monomial.
   * \return Is grevlex-smaller?
   */
  bool isGrevlexSmaller( const Monomial & operant ) const;
  /**
   * \brief Check how this monomial grevlex-compares to another one.
   * \param[in] operant The comparison monomial.
   * \return -1 if this one is smaller, 0 if equal, and 1 if bigger.
   */
  int grevlexComparison( const Monomial & operant ) const;
  /**
   * \brief Check if is dividable by another monomial
   * \param[in] operant The potential divisor.
   * \return Is dividable?
   */
  bool isDividableBy( const Monomial & operant ) const;
  /**
   * \brief Check if is relatively prime to another monomial.
   * \param[in] operant The comparison monomial.
   * \return Is relatively prime?
   */
  bool isRelativelyPrime( const Monomial & operant ) const;
  
  // handy operations
  
  /**
   * \brief Set all exponents to zero
   */
  Monomial & setToOne();
  
  // handy comparisons
  
  /**
   * \brief Check if all exponents are zero.
   * \return Is one?
   */
  bool isOne() const;

private:
  /**
   * The monomial. Each element in the vector represents the exponent for one
   * variable. The length therefore equals the number of dimensions.
   */
  std::vector<unsigned int> _exponents;
  /** The ordering used for the default comparisons */
  Order _order;
  
  /**
   * \brief Check if a monomial has a different number of unknowns, and print
   *        error if so. It does not print a warning in case the order is different!
   * \param[in] operant The comparison monomial.
   * \return Is incompatible?
   */
  bool isIncompatible( const Monomial & operant ) const;
};

}
}

#endif /* POLYJAM_CORE_MONOMIAL_HPP_ */
