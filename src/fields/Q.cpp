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
