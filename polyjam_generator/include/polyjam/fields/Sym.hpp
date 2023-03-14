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

/**
 * \file Sym.hpp
 * \brief Member of the algebraic field Symbolic.
 */

#ifndef POLYJAM_FIELDS_SYM_HPP_
#define POLYJAM_FIELDS_SYM_HPP_

#include <string>
#include <set>
#include <functional>
#include <utility>
#include <memory>

#include <polyjam/fields/Field.hpp>

/**
 * \brief The namespace of this library.
 */
namespace polyjam
{

/**
 * \brief The namespace for the fields.
 */
namespace fields
{

struct PoweredSym
{
  std::string _symbol;
  mutable unsigned int _exponent;
  
  PoweredSym( const std::string & symbol, unsigned int exponent = 1 ) :
      _symbol(symbol),
      _exponent(exponent) {};
  PoweredSym( const PoweredSym & copy ) :
      _symbol(copy._symbol), _exponent(copy._exponent) {};
  
  int compare( const PoweredSym & operant ) const
  {
    if( _symbol < operant._symbol )
      return -1;
    if( _symbol > operant._symbol )
      return  1;
    
    return 0;
  }
  
  bool operator< ( const PoweredSym & operant ) const
  { return compare(operant) <  0; };
  bool operator> ( const PoweredSym & operant ) const
  { return compare(operant) >  0; };
  bool operator<=( const PoweredSym & operant ) const
  { return compare(operant) <= 0; };
  bool operator>=( const PoweredSym & operant ) const
  { return compare(operant) >= 0; };
  bool operator==( const PoweredSym & operant ) const
  { return compare(operant) == 0; };
  bool operator!=( const PoweredSym & operant ) const
  { return compare(operant) != 0; };
};

typedef std::set< PoweredSym,std::greater<PoweredSym> > symProduct_t;

struct SymProduct
{
  symProduct_t _product;
  mutable int _factor;
  
  SymProduct() : _factor(0) {};
  SymProduct(
      const std::string & symbol, unsigned int exponent = 1, int factor = 1 ) :
      _factor(factor)
  {
    _product.insert(PoweredSym(symbol,exponent));
  }
  SymProduct( int factor ) : _factor(factor) {};
  SymProduct( const PoweredSym & poweredSym, int factor = 1 ) : _factor(factor)
  {
    _product.insert(poweredSym);
  }
  SymProduct( const SymProduct & copy ) :
      _product(copy._product), _factor(copy._factor) {};
  
  symProduct_t::iterator begin() const { return _product.begin(); };
  symProduct_t::iterator end()   const { return _product.end();   };
  size_t size() const { return _product.size(); };
  
  void multiply( const PoweredSym & operant )
  {
    std::pair< symProduct_t::iterator, bool> insertion =
        _product.insert(operant);
    
    if( !insertion.second )
      insertion.first->_exponent += operant._exponent;
  };
  
  void multiply( const SymProduct & operant )
  {
    _factor *= operant._factor;
    
    if( _factor == 0 )
      _product.clear();
    else
    {
      for(
          symProduct_t::iterator symIter = operant.begin();
          symIter != operant.end();
          ++symIter )
        multiply( *symIter );
    }
  };
  
  int compare( const SymProduct & operant ) const
  {    
    if( _product.size() < operant._product.size() )
      return -1;
    if( _product.size() > operant._product.size() )
      return  1;
    
    //the elements need to be checked left to right -> phone-book order
    symProduct_t::const_iterator symIter2 = operant._product.begin();
    for(
        symProduct_t::iterator symIter1 = _product.begin();
        symIter1 != _product.end();
        symIter1++ )
    {
      int subComparison = symIter1->compare(*symIter2);
      if( subComparison != 0)
        return subComparison;
      
      //compare the degree as well here!!
      subComparison = symIter1->_exponent - symIter2->_exponent;
      if( subComparison > 0 )
        return  1;
      if( subComparison < 0 )
        return -1;
      
      symIter2++;
    }
    
    return 0;
  };
  
  bool operator< ( const SymProduct & operant ) const
  { return compare(operant) <  0; };
  bool operator> ( const SymProduct & operant ) const
  { return compare(operant) >  0; };
  bool operator<=( const SymProduct & operant ) const
  { return compare(operant) <= 0; };
  bool operator>=( const SymProduct & operant ) const
  { return compare(operant) >= 0; };
  bool operator==( const SymProduct & operant ) const
  { return compare(operant) == 0; };
  bool operator!=( const SymProduct & operant ) const
  { return compare(operant) != 0; };
};

/**
 * The class Symbolic defines a member of the algebraic field Symbolic,
 * including addition, subtraction, multiplication, division.
 */
class Sym : public Field
{
public:

  typedef std::set<SymProduct,std::greater<SymProduct> > symCombination_t;
  typedef std::shared_ptr<symCombination_t> symCombinationPtr;

  /**
   * \brief Constructor.
   */
  Sym();
  /**
   * \brief Constructor.
   * \param[in] name The name.
   */
  Sym( const std::string & name );
  /**
   * \brief Constructor.
   * \param[in] name The name.
   */
  Sym( const char * name );
  /**
   * \brief Constructor for a symbolic constant.
   * \param[in] constant The constant.
   */
  Sym( int constant );
  /**
   * \brief Copy-Constructor (deep).
   * \param[in] copy Original.
   */
  Sym( const Sym & copy );
  /**
   * \brief Copy-Constructor (deep).
   * \param[in] copy Original.
   */
  Sym( const Field* copy );
  /**
   * Destructor.
   */
  virtual ~Sym();

  // output
  
  /** See base-class documentation */
  std::string getString( bool c_version = true ) const;

  /** See base-class documentation */
  std::string getStringSpecial( bool c_version = true ) const;

  // Get constants from this field
  
  /** See base-class documentation */
  Field* zero() const;
  /** See base-class documentation */
  Field* one() const;

  // standard self-operations

  /** See base-class documentation */
  void negation();
  /** See base-class documentation */
  void inversion();
  
  // standard operations
  
  /** See base-class documentation */
  void add( const Field* operant );
  /** See base-class documentation */
  void subtract( const Field* operant );
  /** See base-class documentation */
  void multiply( const Field* operant );
  /** See base-class documentation */
  void divide( const Field* operant );
  
  // comparison operations

  /** See base-class documentation */
  bool isEql( const Field* operant ) const;
  /** See base-class documentation */
  int compare( const Field* operant ) const;
  
private:
  /** The list of signed products of symbols. */
  symCombinationPtr _combination;
};

}
}

#endif /* POLYJAM_FIELDS_SYM_HPP_ */
