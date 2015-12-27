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

#include <polyjam/fields/Zp.hpp>
#include <sstream>
#include <iostream>
#include <math.h>

using namespace std;


polyjam::fields::Zp::Zp(
    bool random, unsigned int characteristic ) : Field(Field::Zp),
    _value(0), _characteristic(characteristic)
{
  if( random )
  {
    double dval = ((double) rand())/ ((double) RAND_MAX) * (_characteristic-1);
    _value = floor(dval+0.5);
  }
};

polyjam::fields::Zp::Zp(
    int value, unsigned int characteristic ) : Field(Field::Zp),
    _characteristic(characteristic)
{
  _value = moduloPrime(value);
}

polyjam::fields::Zp::Zp( const Zp & copy ) : Field(Field::Zp),
    _value(copy._value), _characteristic(copy._characteristic)
{}

polyjam::fields::Zp::Zp( const Field* copy ) : Field(Field::Zp)
{
  if(copy->kind() != Field::Zp)
  {
    _value = 0;
    _characteristic = DEFAULT_CHARACTERISTIC;
    std::cout << "Error: cannot create Zp-copy from non-Zp object" << std::endl;
  }
  else
  {
    Zp * specCopy = (Zp *) copy;
    _value = specCopy->_value;
    _characteristic = specCopy->_characteristic;
  }
}

polyjam::fields::Zp::~Zp()
{}

string
polyjam::fields::Zp::getString( bool c_version ) const
{
  stringstream temp;
  temp << _value;
  return temp.str();
}

unsigned int
polyjam::fields::Zp::characteristic() const
{
  return _characteristic;
}

polyjam::fields::Field*
polyjam::fields::Zp::zero() const
{
  return (new Zp(0));
}

polyjam::fields::Field*
polyjam::fields::Zp::one() const
{
  return (new Zp(1));
}

void
polyjam::fields::Zp::negation()
{  
  int temp = (int) _value;
  temp *= -1;
  _value = moduloPrime(temp);
}

void
polyjam::fields::Zp::inversion()
{
  _value = getMultiplicativeInverse(_value);
}

void
polyjam::fields::Zp::add( const Field* operant )
{
  if(isWrong(operant))
    return;
    
  Zp * specOperant = (Zp *) operant;
  _value = moduloPrime(_value + specOperant->_value);
}

void
polyjam::fields::Zp::subtract( const Field* operant )
{
  if(isWrong(operant))
    return;
    
  Zp * specOperant = (Zp *) operant;
  _value = moduloPrime(((int) _value) - ((int) specOperant->_value));
}

void
polyjam::fields::Zp::multiply( const Field* operant )
{
  if(isWrong(operant))
    return;
    
  Zp * specOperant = (Zp *) operant;
  _value = moduloPrime(((int) _value) * ((int) specOperant->_value));
}

void
polyjam::fields::Zp::divide( const Field* operant )
{
  if(isWrong(operant))
    return;
  
  Zp * specOperant = (Zp *) operant;
  int inverse = (int) getMultiplicativeInverse((int) specOperant->_value);
  _value = moduloPrime( ((int) _value) * inverse );
}

bool
polyjam::fields::Zp::isEql( const Field* operant ) const
{
  if(isWrong(operant))
    return false;
  
  Zp * specOperant = (Zp *) operant;
  return _value == specOperant->_value;
}

int
polyjam::fields::Zp::compare( const Field* operant ) const
{
  cout << "Error: a comparison operation on Zp doesn't make sense!" << endl;
  return false;
}

unsigned int
polyjam::fields::Zp::getMultiplicativeInverse( unsigned int value ) const
{
  if( value == 0 )
  {
    cout << "Error: cannot take multiplicative inverse of 0!" << endl;
    return value;
  }
	int l1 = _characteristic; int y1 = 0;
	int l2 = value; int y2 = 1;

	while( l2 > 1 )
	{
		int modulo = l1 % l2;
		if( modulo == 0 )
			break;
		else
		{
			int temp_l2 = l2;
			int temp_y2 = y2;

			int factor = (l1 - modulo) / l2;

			l2 = l1 - factor * l2;
			y2 = y1 - factor * y2;

			l1 = temp_l2;
			y1 = temp_y2;
		}
	}

	if( l2 > 1 )
	{
		cout << "Error: Unable to compute multiplicative inverse!" << endl;
		return value;
	}

	return moduloPrime(y2);
}

unsigned int
polyjam::fields::Zp::moduloPrime( int value ) const
{
	int result = value % (int) _characteristic;

	if( result < 0 )
		result += _characteristic;

	return (unsigned int) result;
}

bool
polyjam::fields::Zp::isWrong( const Field* operant ) const
{
  if( operant->kind() != Field::Zp )
  {
    cout << "Error: Field not compatible with Zp" << endl;
    return true;
  }
  
  Zp * specOperant = (Zp *) operant;
  
  if( specOperant->_characteristic != _characteristic )
  {
    cout << "Error: Combination of Zp field with different characteristic!";
    cout << endl;
    return true;
  }
  
  return false;
}
