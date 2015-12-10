
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

#include <polyjam/core/Coefficient.hpp>

#include <polyjam/fields/R.hpp>
#include <polyjam/fields/Q.hpp>
#include <polyjam/fields/Zp.hpp>
#include <polyjam/fields/Sym.hpp>

#include <iostream>

using namespace std;

// Constructors

polyjam::core::Coefficient::Coefficient( fields::Field::Kind kind, bool random )
{
  fields::Field * newField = NULL;
  
  switch(kind)
  {
    case fields::Field::R:
    {
      newField = new fields::R(random);
      break;
    }
    case fields::Field::Q:
    {
      newField = new fields::Q(random);
      break;
    }
    case fields::Field::Zp:
    {
      newField = new fields::Zp(random);
      break;
    }
    case fields::Field::Sym:
    {
      newField = new fields::Sym();
      break;
    }
    default:
      break;
  }
  
  _field = FieldPtr(newField);
}

polyjam::core::Coefficient::Coefficient(
    double value ) :
    _field(new fields::R(value))
{}

polyjam::core::Coefficient::Coefficient(
    int numerator, unsigned int denominator ) :
    _field(new fields::Q(numerator,denominator))
{}

polyjam::core::Coefficient::Coefficient(
    int constant, fields::Field::Kind kind )
{
  fields::Field * newField = NULL;
  
  switch(kind)
  {
    case fields::Field::R:
    {
      newField = new fields::R((double) constant);
      break;
    }
    case fields::Field::Q:
    {
      newField = new fields::Q(constant);
      break;
    }
    case fields::Field::Zp:
    {
      newField = new fields::Zp(constant);
      break;
    }
    case fields::Field::Sym:
    {
      newField = new fields::Sym(constant);
      break;
    }
    default:
      break;
  }
  
  _field = FieldPtr(newField);
}

polyjam::core::Coefficient::Coefficient( const std::string & name ) :
    _field(new fields::Sym(name))
{}

polyjam::core::Coefficient::Coefficient( const char * name ) :
    _field(new fields::Sym(name))
{}

polyjam::core::Coefficient::Coefficient( fields::Field * field ) :
    _field(field)
{}

// Destructor

polyjam::core::Coefficient::~Coefficient()
{}

// deep copy stuff

polyjam::core::Coefficient
polyjam::core::Coefficient::clone() const
{
  fields::Field * newField = NULL;
  
  switch( kind() )
  {
    case fields::Field::R:
    {
      newField = new fields::R(_field.get());
      break;
    }
    case fields::Field::Q:
    {
      newField = new fields::Q(_field.get());
      break;
    }
    case fields::Field::Zp:
    {
      newField = new fields::Zp(_field.get());
      break;
    }
    case fields::Field::Sym:
    {
      newField = new fields::Sym(_field.get());
      break;
    }
    default:
      break;
  }
  
  return Coefficient(newField);
}

void
polyjam::core::Coefficient::copy( const Coefficient & coefficient )
{
  fields::Field * newField = NULL;
  
  switch( coefficient.kind() )
  {
    case fields::Field::R:
    {
      newField = new fields::R(coefficient._field.get());
      break;
    }
    case fields::Field::Q:
    {
      newField = new fields::Q(coefficient._field.get());
      break;
    }
    case fields::Field::Zp:
    {
      newField = new fields::Zp(coefficient._field.get());
      break;
    }
    case fields::Field::Sym:
    {
      newField = new fields::Sym(coefficient._field.get());
      break;
    }
    default:
      break;
  }
  
  _field = FieldPtr(newField);
}

// output
void
polyjam::core::Coefficient::print() const
{
  _field->print();
}

string
polyjam::core::Coefficient::getString( bool c_version ) const
{
  return _field->getString( c_version );
}

string
polyjam::core::Coefficient::getStringSpecial( bool c_version ) const
{
  if( kind() != fields::Field::Sym )
  {
    cout << "Error: cannot get special string for non symbolic coefficients";
    return 0;
  }

  fields::Sym * sy = (fields::Sym *) _field.get();
  return sy->getStringSpecial( c_version );
}

polyjam::fields::Field::Kind
polyjam::core::Coefficient::kind() const
{
  return _field->kind();
}

polyjam::core::Coefficient
polyjam::core::Coefficient::zero() const
{
  return Coefficient(_field->zero());
}

polyjam::core::Coefficient
polyjam::core::Coefficient::one() const
{
  return Coefficient(_field->one());
}

unsigned int
polyjam::core::Coefficient::characteristic() const
{
  if( kind() != fields::Field::Zp )
  {
    cout << "Error: cannot retrieve the characteristic";
    cout << " from non Zp coefficient." << endl;
    return 0;
  }
  fields::Zp * zp = (fields::Zp *) _field.get();
  return zp->characteristic();
}

// standard operations

polyjam::core::Coefficient
polyjam::core::Coefficient::negation() const
{
  Coefficient newClone = clone();
  newClone.negationInPlace();
  return newClone;
}

polyjam::core::Coefficient
polyjam::core::Coefficient::inversion() const
{
  Coefficient newClone = clone();
  newClone.inversionInPlace();
  return newClone;
}

polyjam::core::Coefficient
polyjam::core::Coefficient::operator+( const Coefficient & operant ) const
{
  Coefficient newClone = clone();
  newClone += operant;
  return newClone;
}

polyjam::core::Coefficient
polyjam::core::Coefficient::operator-( const Coefficient & operant ) const
{
  Coefficient newClone = clone();
  newClone -= operant;
  return newClone;
}

polyjam::core::Coefficient
polyjam::core::Coefficient::operator*( const Coefficient & operant ) const
{
  Coefficient newClone = clone();
  newClone *= operant;
  return newClone;
}

polyjam::core::Coefficient
polyjam::core::Coefficient::operator/( const Coefficient & operant ) const
{
  Coefficient newClone = clone();
  newClone /= operant;
  return newClone;
}

// in-place operations

polyjam::core::Coefficient &
polyjam::core::Coefficient::negationInPlace()
{
  _field->negation();
  return (*this);
}

polyjam::core::Coefficient &
polyjam::core::Coefficient::inversionInPlace()
{
  _field->inversion();
  return (*this);
}

polyjam::core::Coefficient &
polyjam::core::Coefficient::operator+=( const Coefficient & operant )
{
  _field->add(operant._field.get());
  return (*this);
}

polyjam::core::Coefficient &
polyjam::core::Coefficient::operator-=( const Coefficient & operant )
{
  _field->subtract(operant._field.get());
  return (*this);
}

polyjam::core::Coefficient &
polyjam::core::Coefficient::operator*=( const Coefficient & operant )
{
  _field->multiply(operant._field.get());
  return (*this);
}

polyjam::core::Coefficient &
polyjam::core::Coefficient::operator/=( const Coefficient & operant )
{
  _field->divide(operant._field.get());
  return (*this);
}

// comparisons

bool
polyjam::core::Coefficient::operator==( const Coefficient & operant ) const
{
  return _field->isEql(operant._field.get());
}

bool
polyjam::core::Coefficient::operator!=( const Coefficient & operant ) const
{
  return !(_field->isEql(operant._field.get()));
}

bool
polyjam::core::Coefficient::operator<=( const Coefficient & operant ) const
{
  int result = _field->compare(operant._field.get());
  if( result <= 0 )
    return true;
  return false;
}

bool
polyjam::core::Coefficient::operator>=( const Coefficient & operant ) const
{
  int result = _field->compare(operant._field.get());
  if( result >= 0 )
    return true;
  return false;
}

bool
polyjam::core::Coefficient::operator<( const Coefficient & operant ) const
{
  int result = _field->compare(operant._field.get());
  if( result < 0 )
    return true;
  return false;
}

bool
polyjam::core::Coefficient::operator>( const Coefficient & operant ) const
{
  int result = _field->compare(operant._field.get());
  if( result > 0 )
    return true;
  return false;
}

// handy operations

polyjam::core::Coefficient &
polyjam::core::Coefficient::setToZero()
{
  _field = FieldPtr(_field->zero());
  return (*this);
}

polyjam::core::Coefficient &
polyjam::core::Coefficient::setToOne()
{
  _field = FieldPtr(_field->one());
  return (*this);
}

// handy comparisons

bool
polyjam::core::Coefficient::isZero() const
{
  fields::Field* temp = _field->zero();
  bool value = _field->isEql(temp);
  delete temp;
  return value;
  //return _field->isEql(_field->zero());
}

bool
polyjam::core::Coefficient::isOne() const
{
  fields::Field* temp = _field->one();
  bool value = _field->isEql(temp);
  delete temp;
  return value;
  //return _field->isEql(_field->one());
}
