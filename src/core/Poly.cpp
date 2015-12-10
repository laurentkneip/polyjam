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

#include <polyjam/core/Poly.hpp>
#include <iostream>
#include <sstream>

using namespace std;

//constructors, destructor

polyjam::core::Poly::Poly( const Term & term ) : _sugar(0)
{
  _terms = termsPtr(new terms_t());
  _terms->insert(term);
}

polyjam::core::Poly::~Poly()
{}

//assignment

polyjam::core::Poly &
polyjam::core::Poly::operator=( const Poly & rhs )
{
  if( !leadingTerm().isSimilar(rhs.leadingTerm()) )
  {
    cout << "Error: cannot assign polynomial that has different character";
    cout << endl;
    return (*this);
  }
  
  _terms = rhs._terms;
  _sugar = rhs._sugar;
  return (*this);
}

//output

polyjam::core::Poly
polyjam::core::Poly::clone( bool full ) const
{
  Poly result(leadingTerm().clone(full));
  
  if( _terms->size() > 1 )
  {
    terms_t::iterator iter = _terms->begin();
    ++iter;
    while( iter != _terms->end() )
    {
      result._terms->insert(result._terms->end(), iter->clone(full));
      ++iter;
    }
  }
  
  return result;
}

void
polyjam::core::Poly::copy( const Poly & copy )
{
  if( !leadingTerm().isSimilar(copy.leadingTerm()) )
  {
    cout << "Error: cannot copy from polynomial that has different character";
    cout << endl;
  }
  
  _terms->clear();
  for(
      terms_t::iterator iter = copy._terms->begin();
      iter != copy._terms->end();
      ++iter )
    _terms->insert(_terms->end(), iter->clone() );
  
  _sugar = copy._sugar;
}

polyjam::core::Poly
polyjam::core::Poly::one() const
{
  return Poly(leadingTerm().one());
}

polyjam::core::Poly
polyjam::core::Poly::zero() const
{
  return Poly(leadingTerm().zero());
}

void
polyjam::core::Poly::setDominant( int index )
{
  for(
      terms_t::iterator iter = _terms->begin();
      iter != _terms->end();
      ++iter )
  {
    Term * currentTerm = const_cast<Term*>(&*iter);
    currentTerm->setDominant( index );
  }
}

void
polyjam::core::Poly::print() const
{
	cout << getString(false) << endl;
}

string
polyjam::core::Poly::getString( bool c_version ) const
{
  stringstream result;
  
  for(
	    terms_t::iterator iter = _terms->begin();
	    iter != _terms->end();
	    ++iter )
	{
	  if( iter != _terms->begin() )
	    result << " + ";
	  result << iter->getString(c_version);
	}
	
	return result.str();
}

const polyjam::core::Term &
polyjam::core::Poly::leadingTerm() const
{
  return *(_terms->begin());
}

polyjam::core::Term
polyjam::core::Poly::leadingCoefficient() const
{
  return Term( leadingTerm().coefficient().clone(), leadingTerm().monomial().one() );
}


polyjam::core::Term
polyjam::core::Poly::leadingMonomial() const
{
  return Term( leadingTerm().coefficient().one(), leadingTerm().monomial() );
}

polyjam::core::Poly::terms_t::iterator
polyjam::core::Poly::begin() const
{
  return _terms->begin();
}

polyjam::core::Poly::terms_t::iterator
polyjam::core::Poly::end() const
{
  return _terms->end();
}

size_t
polyjam::core::Poly::size() const
{
  return _terms->size();
}

unsigned int &
polyjam::core::Poly::sugar()
{
  return _sugar;
}

polyjam::core::Coefficient
polyjam::core::Poly::eval( const std::vector<double> & values ) const
{
  polyjam::core::Coefficient result = leadingTerm().coefficient().zero();
  
  if( leadingTerm().coefficient().kind() != fields::Field::R )
  {
    cout << "The eval function is currently only supported for the field R" << endl;
    return result;
  }
  
  for(
      terms_t::iterator termIter = _terms->begin();
      termIter != _terms->end();
      ++termIter )
    result += termIter->coefficient() * Coefficient(termIter->monomial().eval(values));
  
  return result;
}

polyjam::core::Coefficient
polyjam::core::Poly::eval( const std::vector<Coefficient> & values ) const
{
  polyjam::core::Coefficient result = leadingTerm().coefficient().zero();
  for(
      terms_t::iterator termIter = _terms->begin();
      termIter != _terms->end();
      ++termIter )
  {
    const std::vector<unsigned int> & exponents = termIter->monomial().exponents();
    Coefficient coeff = termIter->coefficient().clone();
    
    for( size_t i = 0; i < exponents.size(); i++ )
    {
      if( exponents[i] > 0 )
      {
        for( size_t j = 0; j < exponents[i]; j++ )
          coeff *= values[i];
      }
    }
    
    result += coeff;
  }
  return result;
}

// operations

polyjam::core::Poly
polyjam::core::Poly::differentOrderVersion( Monomial::Order newOrder ) const
{  
  Term firstTerm(leadingTerm().clone());
  firstTerm.setOrder( newOrder );
  Poly newPoly(firstTerm);
  
  if( _terms->size() > 1 )
  {
    terms_t::iterator iter = _terms->begin();
    ++iter;
    
    while( iter != _terms->end() )
    {
      Term nextTerm(iter->clone());
      nextTerm.setOrder(newOrder);
      newPoly._terms->insert(nextTerm);
      ++iter;
    }
  }
  
  return newPoly;
}

polyjam::core::Poly
polyjam::core::Poly::lowerDegreeApproximation( unsigned int maxDegree ) const
{
  Poly result( leadingTerm().zero() );
  
  for(
      terms_t::iterator iter = _terms->begin();
      iter != _terms->end();
      ++iter )
  {
    if( iter->monomial().degree() <= maxDegree )
      result._terms->insert(result._terms->end(),iter->clone());
  }
  
  //check that the size is bigger than one, because then we have to
  //kill the first term
  if( result._terms->size() > 1 )
    result._terms->erase(result._terms->begin());
  
  return result;
}

polyjam::core::Poly
polyjam::core::Poly::negation() const
{
  Poly result( leadingTerm().negation() );
  
  if( _terms->size() > 1 )
  {
    terms_t::iterator iter = _terms->begin();
    ++iter;
    
    while( iter != _terms->end() )
    {
      result._terms->insert(result._terms->end(),iter->negation());
      ++iter;
    }
  }
  
  return result;
}

polyjam::core::Poly
polyjam::core::Poly::operator+( const Poly & operant ) const
{
  Poly result(this->clone());
  result += operant;
  return result;
}

polyjam::core::Poly
polyjam::core::Poly::operator+( const Term & operant ) const
{
  Poly result(this->clone());
  result += operant;
  return result;
}

polyjam::core::Poly
polyjam::core::Poly::operator-( const Poly & operant ) const
{
  Poly result(this->clone());
  result -= operant;
  return result;
}

polyjam::core::Poly
polyjam::core::Poly::operator-( const Term & operant ) const
{
  Poly result(this->clone());
  result -= operant;
  return result;
}

polyjam::core::Poly
polyjam::core::Poly::operator*( const Poly & operant ) const
{
  //Actually the following would be inefficient
  //Poly result(this->clone());
  //result *= operant;
  //return result;
  
  Poly result(this->zero());
  
  for(
      terms_t::iterator iter = operant._terms->begin();
      iter != operant._terms->end();
      ++iter )
  {
    Poly temp(this->clone());
    temp *= (*iter);
    result += temp;
  }
  
  return result;
}

polyjam::core::Poly
polyjam::core::Poly::operator*( const Term & operant ) const
{
  Poly result(this->clone());
  result *= operant;
  return result;
}


// in-place operations

polyjam::core::Poly &
polyjam::core::Poly::differentOrderVersionInPlace( Monomial::Order newOrder )
{
  termsPtr newTerms = termsPtr(new terms_t);
  
  Term firstTerm(leadingTerm().clone());
  firstTerm.setOrder( newOrder );
  newTerms->insert(firstTerm);
  
  if( _terms->size() > 1 )
  {
    terms_t::iterator iter = _terms->begin();
    ++iter;
    
    while( iter != _terms->end() )
    {
      Term nextTerm(iter->clone());
      nextTerm.setOrder(newOrder);
      newTerms->insert(nextTerm);
      ++iter;
    }
  }
  
  swap(newTerms,_terms);
  return (*this);
}

polyjam::core::Poly &
polyjam::core::Poly::lowerDegreeApproximationInPlace( unsigned int maxDegree )
{
  //get a zero term to push_back in case the list becomes empty
  Term zeroTerm(leadingTerm().zero());
  
  termsPtr newTerms = termsPtr(new terms_t);
  
  for(
      terms_t::iterator iter = _terms->begin();
      iter != _terms->end();
      ++iter )
  {
    if( iter->monomial().degree() <= maxDegree )
      newTerms->insert(newTerms->end(),*iter);
  }
    
  if( newTerms->size() == 0 )
    newTerms->insert(zeroTerm);
  
  swap(newTerms,_terms);
  return (*this);
}

polyjam::core::Poly &
polyjam::core::Poly::negationInPlace()
{
  for(
      terms_t::iterator iter = _terms->begin();
      iter != _terms->end();
      ++iter )
  {
    Term * currentTerm = const_cast<Term*>( &*iter );
    currentTerm->negationInPlace();
  }
  
  return (*this);
}

polyjam::core::Poly &
polyjam::core::Poly::operator+=( const Poly & operant )
{
  //go through the terms and add each one
  
  for(
      terms_t::iterator iter = operant._terms->begin();
      iter != operant._terms->end();
      iter++ )
    (*this) += (*iter);
  
  return (*this);
}

polyjam::core::Poly &
polyjam::core::Poly::operator+=( const Term & operant )
{  
  //first check is the term has same characteristics
  if(!leadingTerm().isSimilar(operant))
  {
    cout << "Error: Attempt to add incompatible term" << endl;
    return (*this);
  }
  
  //also check if the zerm is zero, we won't have to add that :)
  if( operant.isZero() )
    return (*this);
  
  //check if this one is zero
  if( isZero() )
  {
    _terms->clear();
    _terms->insert(operant.clone());
    return (*this);
  }
  
  //now insert
  terms_t::iterator insertionPoint = _terms->lower_bound(operant);
  if( insertionPoint == _terms->end() || insertionPoint->monomial() != operant.monomial() )
    _terms->insert(insertionPoint,operant.clone());
  else
  {
    //the item existed already
    Term * currentTerm = const_cast<Term*>(&*(insertionPoint));
    *(currentTerm) += operant;
    
    if( currentTerm->isZero() )
    {
      if( _terms->size() == 1 )
        setToZero();
      else
        _terms->erase(insertionPoint);
    }
  }
  
  return (*this);
}

polyjam::core::Poly &
polyjam::core::Poly::operator-=( const Poly & operant )
{
  //go through the terms and subtract each one
  for(
      terms_t::iterator iter = operant._terms->begin();
      iter != operant._terms->end();
      iter++ )
    (*this) -= (*iter);
  return (*this);
}

polyjam::core::Poly &
polyjam::core::Poly::operator-=( const Term & operant )
{
  //first check is the term has same characteristics
  if(!leadingTerm().isSimilar(operant))
  {
    cout << "Error: Attempt to add incompatible term" << endl;
    return (*this);
  }
  
  //also check if the zerm is zero, we won't have to add that :)
  if( operant.isZero() )
    return (*this);
  
  //check if this one is zero
  if( isZero() )
  {
    _terms->clear();
    _terms->insert(operant.negation());
    return (*this);
  }
  
  //now insert
  terms_t::iterator insertionPoint = _terms->lower_bound(operant);
  if( insertionPoint == _terms->end() || insertionPoint->monomial() != operant.monomial() )
    _terms->insert(insertionPoint,operant.negation());
  else
  {
    //the item existed already
    Term * currentTerm = const_cast<Term*>(&*(insertionPoint));
    *(currentTerm) -= operant;
    
    if( currentTerm->isZero() )
    {
      if( _terms->size() == 1 )
        setToZero();
      else
        _terms->erase(insertionPoint);
    }
  }
  
  return (*this);
}

polyjam::core::Poly &
polyjam::core::Poly::operator*=( const Poly & operant )
{
  Poly result(this->zero());
  
  for(
      terms_t::iterator iter = operant._terms->begin();
      iter != operant._terms->end();
      ++iter )
  {
    Poly temp(this->clone());
    temp *= (*iter);
    result += temp;
  }
  
  // now swap
  swap(result._terms,_terms);
  return (*this);
}

polyjam::core::Poly &
polyjam::core::Poly::operator*=( const Term & operant )
{
  //first check if the term has same characteristics
  if(!leadingTerm().isSimilar(operant))
  {
    cout << "Error: Attempt to multiply by incompatible term" << endl;
    return (*this);
  }
  
  //also check if it is a zero term, we won't have to do anything then :)  
  if( operant.isZero() )
  {
    setToZero();
    return (*this);
  }
  
  //now multiply each term
  for(
      terms_t::iterator iter = _terms->begin();
      iter != _terms->end();
      iter++ )
  {
    Term * currentTerm = const_cast<Term*>(&*iter);
    (*currentTerm) *= operant;
  }
  
  return (*this);
}

// comparisons

bool
polyjam::core::Poly::operator==( const Poly & operant ) const
{
  //We'll use the == function of the term, that will ensure that
  //monomial plus main field compatibility is checked (don't care about multiplicity)
  
  //first check that the order is similar
  if( leadingTerm().monomial().order() != operant.leadingTerm().monomial().order() )
    return false;
  
  //then check the size
  if( _terms->size() != operant._terms->size() )
    return false;
  
  //now loop through the terms
  terms_t::iterator iter1 = _terms->begin();
  terms_t::iterator iter2 = operant._terms->begin();
  
  while( iter1 != _terms->end() )
  {
    if( (*iter1) != (*iter2) )
      return false;
  	++iter1;
  	++iter2;
  }
  
  //if we reach this point, the polynomials are equal
  return true;
}

bool
polyjam::core::Poly::operator!=( const Poly & operant ) const
{
  return !( (*this) == operant );
}

bool
polyjam::core::Poly::isSimilar( const Poly & operant ) const
{
  //first check that the order is similar
  if(leadingTerm().monomial().order() != operant.leadingTerm().monomial().order())
    return false;
  
  //then check the size
  if( _terms->size() != operant._terms->size() )
    return false;
  
  //now loop through the terms, but check only monomials
  terms_t::iterator iter1 = _terms->begin();
  terms_t::iterator iter2 = operant._terms->begin();
  
  while( iter1 != _terms->end() )
  {
    //recall that this is only using the dominant term ...
    if( iter1->monomial() != iter2->monomial() )
      return false;
  	++iter1;
  	++iter2;
  }
  
  //if we reach this point, the polynomials are similar
  return true;
}

// handy operations

polyjam::core::Poly &
polyjam::core::Poly::setToOne()
{
  Term oneTerm(leadingTerm().one());
  _terms->clear();
  _terms->insert(oneTerm);
  return (*this);
}

polyjam::core::Poly &
polyjam::core::Poly::setToZero()
{
  Term zeroTerm(leadingTerm().zero());
  _terms->clear();
  _terms->insert(zeroTerm);
  return (*this);
}

// handy comparisons

bool
polyjam::core::Poly::isZero() const
{
  return _terms->size() == 1 && leadingTerm().isZero();
}

bool
polyjam::core::Poly::isOne() const
{
  return _terms->size() == 1 && leadingTerm().isOne();
}
