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
