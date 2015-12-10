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

#include <polyjam/fields/R.hpp>
#include <sstream>
#include <iostream>
#include <limits>

using namespace std;


polyjam::fields::R::R( bool random ) :
    Field(Field::R), _value(0.0)
{
  if(random)
    _value = (((double) rand()) / ((double) RAND_MAX) - 0.5) * 2.0 * 1000.0;
/*        (((double) rand())/ ((double) RAND_MAX)) *
        numeric_limits<double>::max();*/
}

polyjam::fields::R::R( double value) : Field(Field::R), _value(value)
{}

polyjam::fields::R::R( const R & copy ) : Field(Field::R), _value(copy._value)
{}

polyjam::fields::R::R( const Field* copy ) : Field(Field::R)
{
  if(isWrong(copy))
  {
    _value = 0.0;
    std::cout << "Error: cannot create R-copy from non-R object" << std::endl;
  }
  else
  {
    R * specCopy = (R *) copy;
    _value = specCopy->_value;
  }
}

polyjam::fields::R::~R()
{}

string
polyjam::fields::R::getString( bool c_version ) const
{
  stringstream temp;
  temp << _value;
  return temp.str();
}

polyjam::fields::Field*
polyjam::fields::R::zero() const
{
  return (new R(0.0));
}

polyjam::fields::Field*
polyjam::fields::R::one() const
{
  return (new R(1.0));
}

void
polyjam::fields::R::negation()
{  
  _value *= -1.0;
}

void
polyjam::fields::R::inversion()
{
  if( _value == 0.0 )
  {
    cout << "Error: attempt to invert zero!" << endl;
    return;
  }
  _value = 1.0 / _value;
}

void
polyjam::fields::R::add( const Field* operant )
{
  if(isWrong(operant))
    return;

  R * specOperant = (R *) operant;
  _value += specOperant->_value;
}

void
polyjam::fields::R::subtract( const Field* operant )
{
  if(isWrong(operant))
    return;

  R * specOperant = (R *) operant;
  _value -= specOperant->_value;
}

void
polyjam::fields::R::multiply( const Field* operant )
{
  if(isWrong(operant))
    return;

  R * specOperant = (R *) operant;
  _value *= specOperant->_value;
}

void
polyjam::fields::R::divide( const Field* operant )
{
  if(isWrong(operant))
    return;
  
  R * specOperant = (R *) operant;
  _value /= specOperant->_value;
}

bool
polyjam::fields::R::isEql( const Field* operant ) const
{
  if(isWrong(operant))
    return false;

  R * specOperant = (R *) operant;
  return _value == specOperant->_value;
}

int
polyjam::fields::R::compare( const Field* operant ) const
{
  if(isWrong(operant))
    return 0;
    
  R * specOperant = (R *) operant;
  
  double difference = _value - specOperant->_value;
  
  if( difference < 0.0 )
    return -1;
  if( difference > 0.0 )
    return 1;
  return 0;
}
