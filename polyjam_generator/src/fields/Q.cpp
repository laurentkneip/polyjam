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

#include <polyjam/fields/Q.hpp>
#include <sstream>
#include <iostream>
#include <math.h>
#include <limits>

using namespace std;


polyjam::fields::Q::Q( bool random ) :
    Field(Field::Q), _numerator(0), _denominator(1)
{
  if(random)
  {
    _numerator =
        floor(((double) rand())/ ((double) RAND_MAX) *
        ((double) numeric_limits<int64_t>::max()) );
    _denominator =
        floor(((double) rand())/ ((double) RAND_MAX) *
        ((double) numeric_limits<uint64_t>::max()) );
    clean();
  }
}

polyjam::fields::Q::Q( int64_t numerator, uint64_t denominator ) :
    Field(Field::Q), _numerator(numerator), _denominator(denominator)
{
  clean();
}

polyjam::fields::Q::Q( int value ) :
    Field(Field::Q), _numerator(value), _denominator(1)
{}

polyjam::fields::Q::Q( const Q & copy ) :
    Field(Field::Q), _numerator(copy._numerator),
    _denominator(copy._denominator)
{}

polyjam::fields::Q::Q( const Field* copy ) :
    Field(Field::Q)
{
  if(isWrong(copy))
  {
    _numerator = 0;
    _denominator = 1;
    std::cout << "Error: cannot create Q-copy from non-Q object" << std::endl;
  }
  else
  {
    Q * specCopy = (Q *) copy;
    _numerator = specCopy->_numerator;
    _denominator = specCopy->_denominator;
  }
}

polyjam::fields::Q::~Q()
{}

string
polyjam::fields::Q::getString( bool c_version ) const
{
  stringstream temp;
  temp << "(" << _numerator;
  if(_denominator != 1)
    temp << "/" << _denominator;
  temp << ")";
  return temp.str();
}

polyjam::fields::Field*
polyjam::fields::Q::zero() const
{
  return (new Q(0,1));
}

polyjam::fields::Field*
polyjam::fields::Q::one() const
{
  return (new Q(1,1));
}


void
polyjam::fields::Q::negation()
{  
  _numerator *= -1;
}

void
polyjam::fields::Q::inversion()
{
  if(_denominator == 0)
  {
    cout << "Error: attempt to invert a degenerate Q member!" << endl;
    return;
  }
  
  if(_numerator < 0)
  {
    int64_t temp = _numerator;
    _numerator = -((int64_t) _denominator);
    _denominator = specialAbs(temp);
  }
  else
  {
    if(_numerator > 0)
    {
      int64_t temp = _numerator;
      _numerator = (int64_t) _denominator;
      _denominator = specialAbs(temp);
    }
    else
      cout << "Error: attempt to invert 0 Q member!" << endl;
  }
}

void
polyjam::fields::Q::add( const Field* operant )
{
  if(isWrong(operant))
    return;

  Q * specOperant = (Q *) operant;
  
  if( _denominator == 0 || specOperant->_denominator == 0 )
  {
    cout << "Error: attempt to add degenerate Q member!" << endl;
    return;
  }
    
  int64_t temp =
      _numerator * (int64_t) specOperant->_denominator +
      (int64_t) _denominator * specOperant->_numerator;
  _numerator = temp;
  _denominator *= specOperant->_denominator;
  clean();
}

void
polyjam::fields::Q::subtract( const Field* operant )
{
  if(isWrong(operant))
    return;

  Q * specOperant = (Q *) operant;
  
  if( _denominator == 0 || specOperant->_denominator == 0 )
  {
    cout << "Error: attempt to add degenerate Q member!" << endl;
    return;
  }
    
  int64_t temp =
      _numerator * (int64_t) specOperant->_denominator -
      (int64_t) _denominator * specOperant->_numerator;
  _numerator = temp;
  _denominator *= specOperant->_denominator;
  clean();
}

void
polyjam::fields::Q::multiply( const Field* operant )
{
  if(isWrong(operant))
    return;

  Q * specOperant = (Q *) operant;
  
  if( _denominator == 0 || specOperant->_denominator == 0 )
  {
    cout << "Error: attempt to multiply by degenerate Q member!" << endl;
    return;
  }

  uint64_t gcd1 =
      GCD(specialAbs(specOperant->_numerator),_denominator);
  uint64_t gcd2 =
      GCD(specialAbs(_numerator),specOperant->_denominator);

  int64_t tempOperantNumerator = specOperant->_numerator/((int64_t)gcd1);
  uint64_t tempDenominator = _denominator/gcd1;
  int64_t tempNumerator = _numerator/((int64_t)gcd2);
  uint64_t tempOperantDenominator = specOperant->_denominator/gcd2;

  _numerator = tempNumerator * tempOperantNumerator;
  _denominator = tempDenominator * tempOperantDenominator;

  clean();
}

void
polyjam::fields::Q::divide( const Field* operant )
{
  if(isWrong(operant))
    return;
  
  Q specOperant = *((Q *) operant);
  specOperant.inversion();
  
  this->multiply( (Field*) &specOperant );
}

bool
polyjam::fields::Q::isEql( const Field* operant ) const
{
  if( isWrong(operant) )
    return false;
  
  Q * specOperant = (Q *) operant;
  
  return (_numerator * (int64_t) specOperant->_denominator) ==
      (specOperant->_numerator * (int64_t) _denominator);
}

int
polyjam::fields::Q::compare( const Field* operant ) const
{
  if( isWrong(operant) )
    return 0;
    
  Q * specOperant = (Q *) operant;
  
  int64_t difference =
      (_numerator * (int64_t) specOperant->_denominator) -
      (specOperant->_numerator * (int64_t) _denominator);
  
  if( difference < (int64_t) 0 )
    return -1;
  if( difference > (int64_t) 0 )
    return 1;
  return 0;
}

uint64_t
polyjam::fields::Q::specialAbs( int64_t value ) const
{
  int64_t result;
  if( value < (int64_t) 0 )
    result = -value;
  else
    result = value;

  return (uint64_t) result;
}

uint64_t
polyjam::fields::Q::GCD(uint64_t a, uint64_t b) const
{
  while( 1 )
  {
    a = a % b;
    if( a == (uint64_t) 0 )
      return b;
    b = b % a;
    if( b == (uint64_t) 0 )
      return a;
  }
}

void
polyjam::fields::Q::clean()
{
  uint64_t gcd = GCD(specialAbs(_numerator),_denominator);
  _numerator = _numerator / (int64_t) gcd;
  _denominator = _denominator / gcd;
}
