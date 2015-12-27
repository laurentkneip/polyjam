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

#include <polyjam/core/Term.hpp>
#include <iostream>
#include <sstream>

using namespace std;


polyjam::core::Term::Term(
    const Coefficient & coefficient,
    const Monomial & monomial ) :
     _dominantIndex(0), _monomial(monomial)
{
  _coefficients.push_back(coefficient);
}

polyjam::core::Term::Term(
    const Coefficient & coeff1,
    const Coefficient & coeff2,
    const Monomial & monomial ) :
    _dominantIndex(0), _monomial(monomial)
{
  _coefficients.push_back(coeff2);
  _coefficients.push_back(coeff1);
}

polyjam::core::Term::Term(
    const std::vector<Coefficient> & coefficients,
    const Monomial & monomial ) :
    _coefficients(coefficients), _dominantIndex(0), _monomial(monomial)
{}

polyjam::core::Term::~Term()
{}

polyjam::core::Term
polyjam::core::Term::clone( bool full ) const
{
  if(full)
  {
    std::vector<Coefficient> newCoefficients;
    for( size_t i = 0; i < _coefficients.size(); i++ )
      newCoefficients.push_back(_coefficients[i].clone());
    return Term( newCoefficients, _monomial );
  }
  return Term( _coefficients[_dominantIndex].clone(), _monomial );
}

void
polyjam::core::Term::copy( const Term & copy )
{
  _coefficients.clear();
  for( size_t i = 0; i < copy._coefficients.size(); i++ )
    _coefficients.push_back(copy._coefficients[i].clone());
  _dominantIndex = copy._dominantIndex;
  _monomial = copy._monomial;
}

polyjam::core::Term
polyjam::core::Term::zero( bool full ) const
{
  if(full)
  {
    std::vector<Coefficient> newCoefficients;
    for( size_t i = 0; i < _coefficients.size(); i++ )
      newCoefficients.push_back(_coefficients[i].zero());
    return Term( newCoefficients, _monomial.one() );
  }
  return Term( _coefficients[_dominantIndex].zero(), _monomial.one() );
}

polyjam::core::Term
polyjam::core::Term::one( bool full ) const
{
  if(full)
  {
    std::vector<Coefficient> newCoefficients;
    for( size_t i = 0; i < _coefficients.size(); i++ )
      newCoefficients.push_back(_coefficients[i].one());
    return Term( newCoefficients, _monomial.one() );
  }
  return Term( _coefficients[_dominantIndex].one(), _monomial.one() );
}

bool
polyjam::core::Term::isMultiple() const
{
  return (_coefficients.size() > 1);
}

void
polyjam::core::Term::setDominant( int index )
{
  if( !(index < (int) _coefficients.size()) )
  {
    cout << "Error: cannot set dominant index " << index;
    cout << ". Only " << _coefficients.size() << " present" << endl;
  }
  else
    _dominantIndex = index;
}

void
polyjam::core::Term::print() const
{
  cout << getString(false);
}

string
polyjam::core::Term::getString( bool c_version ) const
{
  stringstream result;
  if( _coefficients[_dominantIndex].kind() == fields::Field::Sym )
  {
    result << "(";
    result << _coefficients[_dominantIndex].getString( c_version );
    result << ")";
  }
  else
    result << _coefficients[_dominantIndex].getString( c_version );
  
  if( !_monomial.isOne() )
  {
    result << "*";
    result << _monomial.getString( c_version );
  }
  return result.str();
}

const polyjam::core::Monomial &
polyjam::core::Term::monomial() const
{
  return _monomial;
}

const polyjam::core::Coefficient &
polyjam::core::Term::coefficient() const
{
  return _coefficients[_dominantIndex];
}

void
polyjam::core::Term::setOrder( Monomial::Order newOrder )
{
  _monomial.setOrder(newOrder);
}

polyjam::core::Term
polyjam::core::Term::negation() const
{
  std::vector<Coefficient> newCoefficients;
  for( size_t i = 0; i < _coefficients.size(); i++ )
    newCoefficients.push_back(_coefficients[i].negation());
  return Term( newCoefficients, _monomial );
}

polyjam::core::Term
polyjam::core::Term::operator+( const Term & operant ) const
{
  if( isWrong(operant) )
    return Term(*this);
  
  if( _monomial != operant._monomial )
  {
    cout << "Error: attempt to add incompatible terms!" << endl;
    return Term(*this);
  }
  
  std::vector<Coefficient> newCoefficients;
  for( size_t i = 0; i < _coefficients.size(); i++ )
    newCoefficients.push_back( _coefficients[i] + operant._coefficients[i] );
  
  return Term( newCoefficients, _monomial );
}

polyjam::core::Term
polyjam::core::Term::operator-( const Term & operant ) const
{
  if( isWrong(operant) )
    return Term(*this);
  
  if( _monomial != operant._monomial )
  {
    cout << "Error: attempt to subtract incompatible terms!" << endl;
    return Term(*this);
  }
  
  std::vector<Coefficient> newCoefficients;
  for( size_t i = 0; i < _coefficients.size(); i++ )
    newCoefficients.push_back( _coefficients[i] - operant._coefficients[i] );
  
  return Term( newCoefficients, _monomial );
}

polyjam::core::Term
polyjam::core::Term::operator*( const Term & operant ) const
{
  if( isWrong(operant) )
    return Term(*this);

  std::vector<Coefficient> newCoefficients;
  for( size_t i = 0; i < _coefficients.size(); i++ )
    newCoefficients.push_back( _coefficients[i] * operant._coefficients[i] );
  
  return Term( newCoefficients, _monomial * operant._monomial );
}

polyjam::core::Term
polyjam::core::Term::operator/( const Term & operant ) const
{
  if( isWrong(operant) )
    return Term(*this);

  std::vector<Coefficient> newCoefficients;
  for( size_t i = 0; i < _coefficients.size(); i++ )
    newCoefficients.push_back( _coefficients[i] / operant._coefficients[i] );
  
  return Term( newCoefficients, _monomial / operant._monomial );
}

polyjam::core::Term &
polyjam::core::Term::negationInPlace()
{
  for( size_t i = 0; i < _coefficients.size(); i++ )
    _coefficients[i].negationInPlace();
  
  return *this;
}

polyjam::core::Term &
polyjam::core::Term::operator+=( const Term & operant )
{
  if( isWrong(operant) )
    return *this;
  
  if( _monomial != operant._monomial )
  {
    cout << "Error: attempt to add incompatible terms!" << endl;
    return *this;
  }
  
  for( size_t i = 0; i < _coefficients.size(); i++ )
    _coefficients[i] += operant._coefficients[i];
    
  return *this;
}

polyjam::core::Term &
polyjam::core::Term::operator-=( const Term & operant )
{
  if( isWrong(operant) )
    return *this;
  
  if( _monomial != operant._monomial )
  {
    cout << "Error: attempt to subtract incompatible terms!" << endl;
    return *this;
  }
  
  for( size_t i = 0; i < _coefficients.size(); i++ )
    _coefficients[i] -= operant._coefficients[i];
    
  return *this;
}

polyjam::core::Term &
polyjam::core::Term::operator*=( const Term & operant )
{
  if( isWrong(operant) )
    return *this;
  
  for( size_t i = 0; i < _coefficients.size(); i++ )
    _coefficients[i] *= operant._coefficients[i];
  _monomial *= operant._monomial;
  
  return *this;
}

polyjam::core::Term &
polyjam::core::Term::operator/=( const Term & operant )
{
  if( isWrong(operant) )
    return *this;
  
  for( size_t i = 0; i < _coefficients.size(); i++ )
    _coefficients[i] /= operant._coefficients[i];
  _monomial /= operant._monomial;
  
  return *this;
}

bool
polyjam::core::Term::operator==( const Term & operant ) const
{
  //use of dominant coefficient only!
  
  if( _monomial.order() != operant.monomial().order() )
    cout << "WARNING: Mixing terms with different default-order" << endl;
  
  return (
      _coefficients[_dominantIndex] == operant._coefficients[operant._dominantIndex] &&
      _monomial == operant._monomial );
}

bool
polyjam::core::Term::operator!=( const Term & operant ) const
{
  //use of dominant coefficient only!
  
  if( _monomial.order() != operant.monomial().order() )
    cout << "WARNING: Mixing terms with different default-order" << endl;

  return (
      _coefficients[_dominantIndex] != operant._coefficients[operant._dominantIndex] ||
      _monomial != operant._monomial );
}

bool
polyjam::core::Term::operator>( const Term & operant ) const
{
  if( _monomial.order() != operant.monomial().order() )
    cout << "ERROR: Comparing terms with different default-order" << endl;
  
  return _monomial > operant._monomial;
}

bool
polyjam::core::Term::isSimilar( const Term & operant ) const
{
  if( _coefficients.size() != operant._coefficients.size() )
    return false;
  
  for( size_t i = 0; i < _coefficients.size(); i++ )
  {
    if( _coefficients[i].kind() != operant._coefficients[i].kind() )
      return false;
  }
  
  if( _monomial.dimensions() != operant._monomial.dimensions() )
    return false;
  if( _monomial.order() != operant._monomial.order() )
    return false;
    
  return true;
}

polyjam::core::Term &
polyjam::core::Term::setToOne()
{
  for( size_t i = 0; i < _coefficients.size(); i++ )
    _coefficients[i].setToOne();
  _monomial.setToOne();
  
  return (*this);
}

polyjam::core::Term &
polyjam::core::Term::setToZero()
{
  for( size_t i = 0; i < _coefficients.size(); i++ )
    _coefficients[i].setToZero();
  _monomial.setToOne();
  
  return (*this);
}

bool
polyjam::core::Term::isOne() const
{
  // only compare the dominant coefficient
  return _coefficients[_dominantIndex].isOne() && _monomial.isOne();
}

bool
polyjam::core::Term::isZero() const
{
  // only compare the dominant coefficient
  return _coefficients[_dominantIndex].isZero();
}

bool
polyjam::core::Term::isWrong( const Term & operant ) const
{
  if( _monomial.order() != operant._monomial.order() )
    cout << "WARNING: Mixing terms with different default-order" << endl;
  if( _coefficients.size() != operant._coefficients.size() )
  {
    cout << "Error: attempt to combine incompatible terms (multiplicity)!" << endl;
    return true;
  }
  return false;
}
