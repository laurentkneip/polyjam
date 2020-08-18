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

#include <polyjam/generator/CMatrix.hpp>

#include <set>
#include <algorithm>
#include <polyjam/math/GaussJordan.hpp>

#include <sstream>

namespace polyjam
{
namespace generator
{

class Comp
{
public:
  bool operator()(
      const polyjam::core::Monomial & m1,
      const polyjam::core::Monomial & m2 )
  {
    return m1 > m2;
  }
};

}
}

using namespace std;

//constructors
polyjam::generator::CMatrix::CMatrix( const polynomials_t & polynomials )
{
  fillMonomials( polynomials );
  fillMatrix( polynomials );
}

polyjam::generator::CMatrix::CMatrix( const polynomials_t & polynomials, const eqs_t & equations )
{
  //extract all polynomials
  std::vector<polynomials_t::const_iterator> polynomialsIters;
  polynomials_t::const_iterator it = polynomials.begin();
  while( it != polynomials.end() )
  {
    polynomialsIters.push_back(it);
    ++it;
  }
  
  polynomials_t allPolynomials;
  core::Coefficient one(polynomials.front()->leadingTerm().coefficient().one());
  for( size_t i = 0; i < equations.size(); i++ )
  {
    core::Term expander( one.clone(), equations[i].second );
    allPolynomials.push_back( new core::Poly( (**(polynomialsIters[equations[i].first]) ) * expander ) );
  }
  
  fillMonomials( allPolynomials );
  fillMatrix( allPolynomials );
  
  //delete all polynomials
  it = allPolynomials.begin();
  while( it != allPolynomials.end() )
  {
    delete (*it);
    ++it;
  }
}

polyjam::generator::CMatrix::CMatrix( const polynomials_t & polynomials, const monomials_t & order )
{
  _monomials = order;
  fillMatrix(polynomials,false);
}

polyjam::generator::CMatrix::CMatrix( const polynomials_t & polynomials, const monomials_t & order, const eqs_t & equations )
{
  //extract all polynomials
  std::vector<polynomials_t::const_iterator> polynomialsIters;
  polynomials_t::const_iterator it = polynomials.begin();
  while( it != polynomials.end() )
  {
    polynomialsIters.push_back(it);
    ++it;
  }
  
  polynomials_t allPolynomials;
  core::Coefficient one(polynomials.front()->leadingTerm().coefficient().one());
  for( size_t i = 0; i < equations.size(); i++ )
  {
    core::Term expander( one.clone(), equations[i].second );
    allPolynomials.push_back( new core::Poly( (**(polynomialsIters[equations[i].first]) ) * expander ) );
  }

  _monomials = order;
  fillMatrix(allPolynomials,false);
  
  //delete all polynomials
  it = allPolynomials.begin();
  while( it != allPolynomials.end() )
  {
    delete (*it);
    ++it;
  }
}

//destructor
polyjam::generator::CMatrix::~CMatrix()
{
  for( size_t row = 0; row < _matrix.size(); row++ )
    delete _matrix[row];
}

//modifiers
void
polyjam::generator::CMatrix::reduce()
{
  math::gaussReduction(_matrix);
}

//accessors
size_t
polyjam::generator::CMatrix::rows()
{
  return _matrix.size();
}

size_t
polyjam::generator::CMatrix::cols()
{
  return _monomials.size();
}

polyjam::core::Coefficient
polyjam::generator::CMatrix::operator()( size_t row, size_t col )
{
  return (*(_matrix[row]))[col].clone();
}

polyjam::generator::CMatrix
polyjam::generator::CMatrix::subMatrix( const std::list<int> & rows )
{
  cmatrix_t subMatrix;
  subMatrix.reserve(rows.size());
  
  for( std::list<int>::const_iterator i = rows.begin(); i != rows.end(); i++ )
  {
    crow_t* newRow = new crow_t();
    newRow->reserve(_matrix[*i]->size());
    for( size_t x = 0; x < _matrix[*i]->size(); x++ )
      newRow->push_back((*(_matrix[*i]))[x].clone());
    subMatrix.push_back(newRow);
  }
  
  return CMatrix(subMatrix,_monomials);
}

polyjam::generator::CMatrix::monomials_t
polyjam::generator::CMatrix::monomials()
{
  return _monomials;
}

polyjam::core::Poly
polyjam::generator::CMatrix::getPolynomial( int row )
{
  core::Poly result = core::Poly(
      core::Term(_matrix[row]->front().zero(),_monomials.front().one()));
  
  for( size_t col = 0; col < _monomials.size(); col++ )
  {
    if( !((*(_matrix[row]))[col].isZero()) )
      result += core::Term((*(_matrix[row]))[col],_monomials[col]);
  }
  
  return result;
}

polyjam::core::Poly
polyjam::generator::CMatrix::getSymbolicPolynomial( int row, const std::string & matrixName )
{
  core::Poly result = core::Poly::zeroS(_monomials.front().dimensions());
  
  for( size_t col = 0; col < _monomials.size(); col++ )
  {
    if( !((*(_matrix[row]))[col].isZero()) )
    {
      std::stringstream coeff;
      coeff << matrixName << "(" << row << "," << col << ")";
      result += core::Term(core::Coefficient(coeff.str()),_monomials[col]);
    }
  }
  
  return result;
}

polyjam::core::Poly
polyjam::generator::CMatrix::getSymbolicPolynomial2( int row )
{
  core::Poly result = core::Poly::zeroS(_monomials.front().dimensions());
  
  for( size_t col = 0; col < _monomials.size(); col++ )
  {
    if( !((*(_matrix[row]))[col].isZero()) )
    {
      result += core::Term(core::Coefficient((int) (col + 1),fields::Field::Sym),_monomials[col]);
    }
  }
  
  return result;
}

polyjam::generator::CMatrix::polynomials_t
polyjam::generator::CMatrix::getPolynomials()
{
  polynomials_t polynomials;
  
  for( size_t i = 0; i < _matrix.size(); i++ )
    polynomials.push_back(new core::Poly(getPolynomial(i)));
  
  return polynomials;
}

polyjam::generator::CMatrix::polynomials_t
polyjam::generator::CMatrix::getSymbolicPolynomials( const std::string & matrixName )
{
  polynomials_t polynomials;
  
  for( size_t i = 0; i < _matrix.size(); i++ )
    polynomials.push_back(new core::Poly(getSymbolicPolynomial(i,matrixName)));
  
  return polynomials;
}

polyjam::generator::CMatrix::polynomials_t
polyjam::generator::CMatrix::getSymbolicPolynomials2()
{
  polynomials_t polynomials;
  
  for( size_t i = 0; i < _matrix.size(); i++ )
    polynomials.push_back(new core::Poly(getSymbolicPolynomial2(i)));
  
  return polynomials;
}

//verifiers
bool
polyjam::generator::CMatrix::contains( const polynomials_t & polynomials )
{
  bool consolePrint = false;
  
  //get the polynomials of this matrix
  polynomials_t allPolynomials = getPolynomials();
  
  bool allFound = true;
  polynomials_t::const_iterator dontMissIter = polynomials.begin();
  int polyIndex = 1;
  
  while( dontMissIter != polynomials.end() )
  {
    if(consolePrint)
      std::cout << "Trying to find polynomial number " << polyIndex++ << "/" << polynomials.size() << " ... ";
    
    bool found = false;
    polynomials_t::iterator haveIter = allPolynomials.begin();
    while( haveIter != allPolynomials.end() )
    {
      if( (**dontMissIter) == (**haveIter) )
      {
        found = true;
        break;
      }
      
      haveIter++;
    }
    
    if(!found)
    {
      if(consolePrint)
        std::cout << "not found!" << std::endl;
      allFound = false;
      break;
    }
    else
    {
      if(consolePrint)
        std::cout << "found!" << std::endl;
    }
    
    dontMissIter++;
  }
  
  //delete the polynomials of this matrix again
  polynomials_t::iterator it = allPolynomials.begin();
  while( it != allPolynomials.end() )
  {
    delete (*it);
    ++it;
  }
  
  return allFound;
}

void
polyjam::generator::CMatrix::visualize( bool save )
{
  math::visualizeMatrix(_matrix,true,save);
}

//internal
polyjam::generator::CMatrix::CMatrix( cmatrix_t & matrix, monomials_t & monomials ) :
    _matrix(matrix),
    _monomials(monomials)
{};

void
polyjam::generator::CMatrix::fillMonomials( const polynomials_t & polynomials )
{
  //generate a tree of all monomials appearing in the problem
  //the good thing about a tree is that insertion of new elements is easy,
  //almost like in a list (exponential time). The other thing is that search
  //is also expontential in time.
  std::set<core::Monomial,Comp> monomialTree;
  
  polynomials_t::const_iterator polyIter = polynomials.begin();
  while( polyIter != polynomials.end() )
  {
    std::set<core::Term>::iterator termIter = (*polyIter)->begin();
    while( termIter != (*polyIter)->end() )
    {
      monomialTree.insert( termIter->monomial() );
      ++termIter;
    }
    ++polyIter;
  }
  
  //get the number of cols
  int cols = monomialTree.size();
  
  //efficiently convert the tree into a vector
  _monomials.reserve(cols);
  std::set<core::Monomial,Comp>::iterator treeIter = monomialTree.begin();
  while( treeIter != monomialTree.end() )
  {
    _monomials.push_back(*treeIter);
    ++treeIter;
  }
}

void
polyjam::generator::CMatrix::fillMatrix( const polynomials_t & polynomials, bool quickOrdering )
{
  int rows = polynomials.size();
  int cols = _monomials.size();
  
  //setup the matrix with zero coefficients
  core::Coefficient zero(
      polynomials.front()->leadingTerm().coefficient().zero() );
  
  _matrix.reserve(rows);
  for( int i = 0; i < rows; i++ )
  {
    crow_t* newRow = new crow_t();
    newRow->reserve(cols);
    for( int x = 0; x < cols; x++ )
      newRow->push_back(zero.clone());
    _matrix.push_back(newRow);
  }
  
  //now add all the terms by retrieving the index in the matrix via binary
  //search in the vector
  if( quickOrdering )
  {
    polynomials_t::const_iterator polyIter = polynomials.begin();
    int row = 0;
    Comp comp;
    while( polyIter != polynomials.end() )
    {
      std::set<core::Term>::iterator termIter = (*polyIter)->begin();
      monomials_t::iterator colIter = _monomials.begin();
      while( termIter != (*polyIter)->end() )
      {
        //find the element in the vector pe_monomials
        colIter =  std::lower_bound(
            colIter, _monomials.end(), termIter->monomial(), comp );
        int col = colIter - _monomials.begin();
        
        (*(_matrix[row]))[col] = termIter->coefficient().clone();
        ++termIter;
      }
      ++polyIter;
      ++row;
    }
  }
  else
  {
    polynomials_t::const_iterator polyIter = polynomials.begin();
    int row = 0;

    while( polyIter != polynomials.end() )
    {
      std::set<core::Term>::iterator termIter = (*polyIter)->begin();
      
      while( termIter != (*polyIter)->end() )
      {
        monomials_t::iterator colIter = _monomials.begin();

        while( *colIter != termIter->monomial() )
          colIter++;
        int col = colIter - _monomials.begin();
        
        (*(_matrix[row]))[col] = termIter->coefficient().clone();
        ++termIter;
      }
      ++polyIter;
      ++row;
    }
  }
}

