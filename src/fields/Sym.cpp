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

#include <polyjam/fields/Sym.hpp>
#include <sstream>
#include <iostream>

using namespace std;


polyjam::fields::Sym::Sym() : Field(Field::Sym)
{
  _combination = symCombinationPtr( new symCombination_t() );
  _combination->insert(SymProduct());
}

polyjam::fields::Sym::Sym( const std::string & name ) :
    Field(Field::Sym)
{
  _combination = symCombinationPtr( new symCombination_t() );
  _combination->insert(SymProduct(name));
}

polyjam::fields::Sym::Sym( const char * name ) :
    Field(Field::Sym)
{
  _combination = symCombinationPtr( new symCombination_t() );
  _combination->insert(SymProduct(std::string(name)));
}

polyjam::fields::Sym::Sym( int constant ) :
    Field(Field::Sym)
{
  _combination = symCombinationPtr( new symCombination_t() );
  _combination->insert(SymProduct(constant));
}

polyjam::fields::Sym::Sym( const Sym & copy ) :
    Field(Field::Sym),
    _combination( symCombinationPtr(new symCombination_t(*(copy._combination))) )
{}

polyjam::fields::Sym::Sym(
    const Field* copy ) : Field(Field::Sym)
{
  if(isWrong(copy))
  {
    cout << "Error: cannot create Sym-copy of non-Sym object" << endl;
    _combination = symCombinationPtr( new symCombination_t() );
    _combination->insert(SymProduct());
  }
  else
  {
    Sym * specCopy = (Sym *) copy;
    _combination =
        symCombinationPtr(new symCombination_t(*(specCopy->_combination)));
  }
}

polyjam::fields::Sym::~Sym()
{}

string
polyjam::fields::Sym::getString( bool c_version ) const
{
  stringstream temp;
  
  for(
      symCombination_t::iterator iter = _combination->begin();
      iter != _combination->end();
      ++iter )
  {
    if( iter != _combination->begin() && iter->_factor > 0 )
      temp << "+";
    
    bool firstAdded = false;
    
    if( abs(iter->_factor) != 1 || iter->_product.empty() )
    {
      temp << iter->_factor;
      firstAdded = true;
    }
    else
    {
      if(iter->_factor < 0)
        temp << "-";
    }
    
    for(
        symProduct_t::iterator iter2 = iter->_product.begin();
        iter2 != iter->_product.end();
        ++iter2 )
    {      
      if(firstAdded == true)
        temp << "*";
      
      if( iter2->_exponent == 1 )
        temp << iter2->_symbol;
      else
      {          
        if( c_version)
          temp << "pow(" << iter2->_symbol << "," << iter2->_exponent << ")";
        else
          temp << iter2->_symbol << "^" << iter2->_exponent;
      }
      
      firstAdded = true;
    }
  }
  
  return temp.str();
}

string
polyjam::fields::Sym::getStringSpecial( bool c_version ) const
{
  stringstream temp;
  
  for(
      symCombination_t::iterator iter = _combination->begin();
      iter != _combination->end();
      ++iter )
  {
    if( iter != _combination->begin() )
      temp << "+";
    
    bool firstAdded = false;
    
    if( iter->_product.empty() )
    {
      temp << "1";
      firstAdded = true;
    }
    
    for(
        symProduct_t::iterator iter2 = iter->_product.begin();
        iter2 != iter->_product.end();
        ++iter2 )
    {      
      if(firstAdded == true)
        temp << "*";
      
      if( iter2->_exponent == 1 )
        temp << iter2->_symbol;
      else
      {          
        if( c_version)
          temp << "pow(" << iter2->_symbol << "," << iter2->_exponent << ")";
        else
          temp << iter2->_symbol << "^" << iter2->_exponent;
      }
      
      firstAdded = true;
    }
  }
  
  return temp.str();
}

polyjam::fields::Field*
polyjam::fields::Sym::zero() const
{  
  Sym * s = new Sym();
  return s;
}

polyjam::fields::Field*
polyjam::fields::Sym::one() const
{
  Sym * s = new Sym(1);
  return s;
}

void
polyjam::fields::Sym::negation()
{  
  for(
      symCombination_t::iterator iter = _combination->begin();
      iter != _combination->end();
      ++iter )
    iter->_factor *= -1;
}

void
polyjam::fields::Sym::inversion()
{
  cout << "Error: inversion of symbolic coefficient is not yet supported!";
  cout << endl;
}

void
polyjam::fields::Sym::add( const Field* operant )
{
  if(isWrong(operant) )
    return;

  Sym * specOperant = (Sym *) operant;
  symCombinationPtr origin = specOperant->_combination;
  
  for(
      symCombination_t::iterator iter = origin->begin();
      iter != origin->end();
      ++iter )
  {
    if( iter->_factor != 0 )
    {
      std::pair<symCombination_t::iterator,bool> insertion =
          _combination->insert(*iter);
      if( !insertion.second )
      {
        insertion.first->_factor += iter->_factor;
        if( insertion.first->_factor == 0 )
          _combination->erase(insertion.first);
      }
    }
  }
  
  if( _combination->empty() )
    _combination->insert(SymProduct());
}

void
polyjam::fields::Sym::subtract( const Field* operant )
{
  if(isWrong(operant) )
    return;

  Sym * specOperant = (Sym *) operant;
  symCombinationPtr origin = specOperant->_combination;
  
  for(
      symCombination_t::iterator iter = origin->begin();
      iter != origin->end();
      ++iter )
  {
    if( iter->_factor != 0 )
    {
      std::pair<symCombination_t::iterator,bool> insertion =
          _combination->insert(*iter);
      if( !insertion.second )
      {
        insertion.first->_factor -= iter->_factor;
        if( insertion.first->_factor == 0 )
          _combination->erase(insertion.first);
      }
      else
        insertion.first->_factor *= -1;
    }
  }
  
  if( _combination->empty() )
    _combination->insert(SymProduct());
}

void
polyjam::fields::Sym::multiply( const Field* operant )
{
  if(isWrong(operant) )
    return;

  Sym * specOperant = (Sym *) operant;
  
  symCombinationPtr target = symCombinationPtr(new symCombination_t());
  symCombinationPtr list1 = _combination;
  symCombinationPtr list2 = specOperant->_combination;
  
  for(
      symCombination_t::iterator term1 = list1->begin();
      term1 != list1->end();
      ++term1 )
  {
    for(
        symCombination_t::iterator term2 = list2->begin();
        term2 != list2->end();
        ++term2 )
    {
      SymProduct copy = *term2;
      copy.multiply(*term1);
      
      if( copy._factor != 0 )
      {
        std::pair<symCombination_t::iterator,bool> insertion =
            target->insert(copy);
        if( !insertion.second )
        {
          insertion.first->_factor += copy._factor;
          if( insertion.first->_factor == 0 )
            target->erase(insertion.first);
        }
      }
    }
  }
  
  if( target->empty() )
    target->insert(SymProduct());
  
  swap(target,_combination);
}

void
polyjam::fields::Sym::divide( const Field* operant )
{
  cout << "Error: division of symbolic coefficients is not yet supported!";
  cout << endl;
}

bool
polyjam::fields::Sym::isEql( const Field* operant ) const
{
  if(isWrong(operant))
    return false;

  Sym * specOperant = (Sym *) operant;
  
  symCombinationPtr list1 = _combination;
  symCombinationPtr list2 = specOperant->_combination;
  
  //first compare the sizes of the lists
  if( list1->size() != list2->size() )
    return false;
  
  //then compare the elements of the list
  symCombination_t::iterator term2 = list2->begin();
  
  for(
      symCombination_t::iterator term1 = list1->begin();
      term1 != list1->end();
      ++term1 )
  {
    //Compare the "sign"
    if( term1->_factor != term2->_factor )
      return false;
    
    //Compare the compatibility of the product (equality)
    if( (*term1) != (*term2) )
      return false;
    
    ++term2;
  }
  
  return true;
}

int
polyjam::fields::Sym::compare( const Field* operant ) const
{
  cout << "Error: Comparison operation with Symbolic not allowed" << endl;
  return 0;
}
