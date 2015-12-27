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

#include <polyjam/core/Monomial.hpp>
#include <sstream>
#include <iostream>
#include <math.h>


using namespace std;

// constructors

polyjam::core::Monomial::Monomial( size_t dimensions, Order order ) :
    _order(order)
{
  _exponents.resize(dimensions,0);
}

polyjam::core::Monomial::Monomial(
    size_t dimensions, const unsigned int * exponents, Order order ) :
    _order(order)
{
  _exponents.resize(dimensions,0);
  for( size_t i = 0; i < dimensions; i++ )
    _exponents[i] = exponents[i];
}

polyjam::core::Monomial::Monomial(
    const vector<unsigned int> & exponents, Order order ) :
    _exponents(exponents), _order(order)
{}

polyjam::core::Monomial::Monomial(
    size_t dimensions, size_t indexOne, Order order ) :
    _order(order)
{
  _exponents.clear();
  _exponents.resize(dimensions,0);
  if( indexOne > dimensions )
  {
  	cout << "Error: cannot set dimension " << indexOne << " to one.";
  	cout << "Chose from [1:" << dimensions << "]" << endl;
  }
  else
  {
  	if( indexOne > 0 )
  	  _exponents[indexOne-1] = 1;
  }
}

// destructors

polyjam::core::Monomial::~Monomial()
{}

// output

void
polyjam::core::Monomial::print() const
{
  cout << getString(false);
}

string
polyjam::core::Monomial::getString( bool c_version ) const
{
  stringstream result;
  bool firstPrinted = false;
  
  if( c_version )
  {
    for(size_t i = 0; i < _exponents.size(); i++)
    {
      if( _exponents[i] != 0 )
      {
        if( firstPrinted )
          result << "*";
      
        if( _exponents[i] > 1 )
          result << "pow(x_" << (i+1) << "," << _exponents[i] << ")";
        else
          result << "x_" << (i+1);
        
        firstPrinted = true;
      }
    }
  }
  else
  {
    for(size_t i = 0; i < _exponents.size(); i++)
    {
      if( _exponents[i] != 0 )
      {
        if( firstPrinted )
          result << "*";
      
        result << "x_" << (i+1);
        if( _exponents[i] > 1 )
          result << "^" << _exponents[i];
        firstPrinted = true;
      }
    }
  }
  
  if( !firstPrinted )
    result << "1";
  return result.str();
}

string
polyjam::core::Monomial::getAlpha() const
{
  stringstream alpha;
  for( size_t i = 0; i < _exponents.size(); i++ )
    alpha << _exponents[i];
  return alpha.str();
}

const std::vector<unsigned int> &
polyjam::core::Monomial::exponents() const
{
  return _exponents;
}

polyjam::core::Monomial
polyjam::core::Monomial::one() const
{
  return Monomial(_exponents.size(),_order);
}

double
polyjam::core::Monomial::eval( const std::vector<double> & values ) const
{
  if( _exponents.size() != values.size() )
  {
    cout << "Error: wrong number of values in eval function" << endl;
    return 0.0;
  }
  
  double result = 1.0;
  for( int i = 0; i < (int) values.size(); i++ )
    result *= pow( values[i], _exponents[i] );
  return result;
}

// properties

unsigned int
polyjam::core::Monomial::degree() const
{
  unsigned int totalDegree = 0;

  for( size_t i = 0; i < _exponents.size(); i++ )
    totalDegree += _exponents[i];

  return totalDegree;
}

unsigned int
polyjam::core::Monomial::dimensions() const
{
  return _exponents.size();
}

polyjam::core::Monomial::Order
polyjam::core::Monomial::order() const
{
  return _order;
}

void
polyjam::core::Monomial::setOrder( Order newOrder )
{
  _order = newOrder;
}

// basic operations

polyjam::core::Monomial
polyjam::core::Monomial::leastCommonMultiple(
    const Monomial & operant ) const
{
  if( isIncompatible(operant) )
    return (*this);

  Monomial result(operant);

  for( size_t i = 0; i < _exponents.size(); i++ )
  {
    if( _exponents[i] > operant._exponents[i] )
      result._exponents[i] = _exponents[i];
  }

  return result;
}

polyjam::core::Monomial
polyjam::core::Monomial::operator*(const Monomial & operant) const
{
  if( isIncompatible(operant) )
    return (*this);

  Monomial result(*this);
  for( size_t i = 0; i < result._exponents.size(); i++ )
    result._exponents[i] += operant._exponents[i];

  return result;
}

polyjam::core::Monomial
polyjam::core::Monomial::operator/(const Monomial & operant) const
{
  if( isIncompatible(operant) )
    return (*this);

  if( !isDividableBy(operant) )
  {
    cout << "Error: Monomial not dividable!" << endl;
    return (*this);
  }
  
  Monomial result(*this);
  for( size_t i = 0; i < result._exponents.size(); i++ )
    result._exponents[i] -= operant._exponents[i];

  return result;
}

polyjam::core::Monomial &
polyjam::core::Monomial::operator*=(const Monomial & operant)
{
  if( isIncompatible(operant) )
    return (*this);

  for( size_t i = 0; i < _exponents.size(); i++ )
    _exponents[i] += operant._exponents[i];
  return (*this);
}

polyjam::core::Monomial &
polyjam::core::Monomial::operator/=(const Monomial & operant)
{
  if( isIncompatible(operant) )
    return (*this);

  if( !isDividableBy(operant) )
  {
    cout << "Error: Monomial not dividable!" << endl;
    return (*this);
  }
  
  for( size_t i = 0; i < _exponents.size(); i++ )
    _exponents[i] -= operant._exponents[i];
  return (*this);
}

// comparisons

bool
polyjam::core::Monomial::operator==( const Monomial & operant ) const
{
  if( lexComparison( operant ) == 0 )
    return true;
  return false;
}

bool
polyjam::core::Monomial::operator!=( const Monomial & operant ) const
{
  if( lexComparison( operant ) != 0 )
    return true;
  return false;
}

bool
polyjam::core::Monomial::operator>=( const Monomial & operant ) const
{
  if( comparison(operant,_order) >= 0 )
    return true;
  return false;
}

bool
polyjam::core::Monomial::operator<=( const Monomial & operant ) const
{
  if( comparison(operant,_order) <= 0 )
    return true;
  return false;
}

bool
polyjam::core::Monomial::operator>( const Monomial & operant ) const
{
  if( comparison(operant,_order) > 0 )
    return true;
  return false;
}

bool
polyjam::core::Monomial::operator<( const Monomial & operant ) const
{
  if( comparison(operant,_order) < 0 )
    return true;
  return false;
}

bool
polyjam::core::Monomial::isBigger(
    const Monomial & operant, const Order & order ) const
{
  if( comparison( operant, order ) > 0 )
    return true;
  return false;
}

bool
polyjam::core::Monomial::isSmaller(
    const Monomial & operant, const Order & order ) const
{
  if( comparison( operant, order ) < 0 )
    return true;
  return false;
}

int
polyjam::core::Monomial::comparison(
    const Monomial & operant, const Order & order ) const
{  
  switch(order)
  {
    case Monomial::LEX:
    {
      return lexComparison(operant);
      break;
    }
    case Monomial::REVLEX:
    {
      return revlexComparison(operant);
      break;
    }
    case Monomial::GRLEX:
    {
      return grlexComparison(operant);
      break;
    }
    case Monomial::GREVLEX:
    {
      return grevlexComparison(operant);
      break;
    }
    default:
    {
      //impossible
      return 0;
    }
  }
  
  return 0;
}

bool
polyjam::core::Monomial::isLexBigger( const Monomial & operant ) const
{
  if( lexComparison(operant) > 0 )
    return true;
  return false;
}

bool
polyjam::core::Monomial::isLexSmaller( const Monomial & operant ) const
{
  if( lexComparison(operant) < 0 )
    return true;
  return false;
}

int
polyjam::core::Monomial::lexComparison( const Monomial & operant ) const
{
  if( isIncompatible(operant) )
    return 0;
  
  for( size_t i = 0; i < _exponents.size(); i++ )
  {
    if( _exponents[i] !=  operant._exponents[i] )
    {
      if( _exponents[i] > operant._exponents[i] )
        return 1;
      else
        return -1;
    }
  }
  
  return 0;
}

bool
polyjam::core::Monomial::isRevlexBigger( const Monomial & operant ) const
{
  if( revlexComparison(operant) > 0 )
    return true;
  return false;
}

bool
polyjam::core::Monomial::isRevlexSmaller( const Monomial & operant ) const
{
  if( revlexComparison(operant) < 0 )
    return true;
  return false;
}

int
polyjam::core::Monomial::revlexComparison( const Monomial & operant ) const
{
  if( isIncompatible(operant) )
    return 0;
  
  // This is a bit tricky here, because this definition agrees with RevLex
  // in Macaulay.
  // We reuse this in GrevLex, by taking RevLex only if the grades are same.
  // This does however not correspond with GrevLex in Macaulay!
  // In order to correspond, we have to start from the back!
  
  /*for( size_t i = 0; i < _exponents.size(); i++ )
  {    
    if( _exponents[i] != operant._exponents[i] )
    {
      if( _exponents[i] < operant._exponents[i] )
        return 1;
      else
        return -1;
    }
  }*/
  
  int i = _exponents.size();
  do
  {
    i--;
    if( _exponents[i] != operant._exponents[i] )
    {
      if( _exponents[i] < operant._exponents[i] )
        return 1;
      else
        return -1;
    }
  } while( i != 0);
  
  return 0;
}

bool
polyjam::core::Monomial::isGrlexBigger( const Monomial & operant ) const
{
  if( grlexComparison(operant) > 0 )
    return true;
  return false;
}

bool
polyjam::core::Monomial::isGrlexSmaller( const Monomial & operant ) const
{
  if( grlexComparison(operant) < 0 )
    return true;
  return false;
}

int
polyjam::core::Monomial::grlexComparison(
    const Monomial & operant ) const
{
  if( isIncompatible(operant) )
    return 0;

  unsigned int thisDegree = degree();
  unsigned int thatDegree = operant.degree();

  if( thisDegree == thatDegree )
    return lexComparison(operant);
  else
  {
    if( thisDegree > thatDegree )
      return 1;
    else
      return -1;
  }

  return 0;
}

bool
polyjam::core::Monomial::isGrevlexBigger( const Monomial & operant ) const
{
  if( grevlexComparison(operant) > 0 )
    return true;
  return false;
}

bool
polyjam::core::Monomial::isGrevlexSmaller( const Monomial & operant ) const
{
  if( grevlexComparison(operant) < 0 )
    return true;
  return false;
}

int
polyjam::core::Monomial::grevlexComparison(
    const Monomial & operant ) const
{
  if( isIncompatible(operant) )
    return 0;

  unsigned int thisDegree = degree();
  unsigned int thatDegree = operant.degree();

  if( thisDegree == thatDegree )
    return revlexComparison(operant);
  else
  {
    if( thisDegree > thatDegree )
      return 1;
    else
      return -1;
  }

  return 0;
}

bool
polyjam::core::Monomial::isDividableBy(
    const Monomial & operant ) const
{
  if( isIncompatible(operant) )
    return false;

  for( size_t i = 0; i < operant._exponents.size(); i++ )
  {
    if( operant._exponents[i] > _exponents[i] )
      return false;
  }

  return true;
}

bool
polyjam::core::Monomial::isRelativelyPrime(
    const Monomial & operant ) const
{
  if( isIncompatible(operant) )
    return false;

  Monomial lcm = leastCommonMultiple(operant);
  Monomial product = (*this) * operant;
  return lcm == product;
}

// handy operations

polyjam::core::Monomial &
polyjam::core::Monomial::setToOne()
{
  for( size_t i = 0; i < _exponents.size(); i++ )
    _exponents[i] = 0;
  
  return (*this);
}

// handy comparisons

bool
polyjam::core::Monomial::isOne() const
{
  if( degree() == 0 )
    return true;
  return false;
}

// private

bool
polyjam::core::Monomial::isIncompatible( const Monomial & operant ) const
{
  if( _exponents.size() != operant._exponents.size() )
  {
    cout << "Error: Monomial operation with incompatible dimensions!" << endl;
    return true;
  }
  return false;
}
