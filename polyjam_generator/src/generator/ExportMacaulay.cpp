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

#include <polyjam/generator/ExportMacaulay.hpp>
#include <polyjam/fields/Zp.hpp>

using namespace std;

polyjam::generator::ExportMacaulay::ExportMacaulay()
{
  _numberVariables = 0;
  _characteristic = 0;
  _numberEquations = 0;
}

polyjam::generator::ExportMacaulay::~ExportMacaulay()
{
}

void
polyjam::generator::ExportMacaulay::addPoly( const core::Poly & poly )
{  
  const core::Coefficient & coeff = poly.leadingTerm().coefficient();

  if( coeff.kind() != fields::Field::Zp )
  {
    cout << "Error: cannot create Macaulay file from non-Zp polynomials" << endl;
    cout << "Maybe you forgot to properly set the dominant index?" << endl;
    return;
  }

  if( _numberEquations == 0 )
  {
    _numberVariables = poly.leadingTerm().monomial().dimensions();
    _characteristic = coeff.characteristic();
  }
  else
  {    
    if(
        _numberVariables != poly.leadingTerm().monomial().dimensions() ||
        _characteristic != coeff.characteristic() )
    {
      cout << "Error: dimensions or characteristic don't match" << endl;
      return;
    }
  }
  
  _numberEquations++;
  
  if( _numberEquations > 1 )
    _set << " || ";
  _set << "f" << _numberEquations;
    
  _equations << "f" << _numberEquations << "=";
  _equations << poly.getString( false );
  _equations << ";" << endl;
}

void
polyjam::generator::ExportMacaulay::write( const std::string & fileName )
{
  ofstream file;
  file.open(fileName.c_str());
  
  file << "-- Macaulay2 code template for gbsolver generator"                    << endl;
  file << "-- by Laurent Kneip 2013"                                             << endl;
  file << ""                                                                     << endl;
  file << "KK = ZZ/" << _characteristic                                          << endl;
  file << "R = KK[x_1..x_" << _numberVariables << ", MonomialOrder=>GRevLex]"    << endl;
  file << ""                                                                     << endl;
  file << "-- equations"                                                         << endl;
  file << ""                                                                     << endl;
  file << _equations.rdbuf()                                                     << endl;
  file << "f = (" << _set.rdbuf() << ");"                                        << endl;
  file << ""                                                                     << endl;
  file << "-- computation of the basis"                                          << endl;
  file << ""                                                                     << endl;
  file << "gbTrace = 0;"                                                         << endl;
  file << "I1 = ideal(f);"                                                       << endl;
  file << "dm = dim I1;"                                                         << endl;
  file << "dg = degree I1;"                                                      << endl;
  file << ""                                                                     << endl;
  file << "--printing of the output"                                             << endl;
  //file << ""                                                                     << endl;
  //file << "print \"dim:\""                                                       << endl;
  //file << "print dm;"                                                            << endl;
  //file << "print \"deg:\""                                                       << endl;
  //file << "print dg;"                                                            << endl;
  file << ""                                                                     << endl;
  //file << "gens gb I1"                                                           << endl;
  file << "A = R/I1;"                                                            << endl;
  file << "Ab = basis A;"                                                        << endl;
  file << ""                                                                     << endl;
  //file << "print \"********\""                                                   << endl;
  file << "print Ab;"                                                            << endl;
  //file << "print \"********\""                                                   << endl;
  file << ""                                                                     << endl;
  file << "exit 0"                                                               << endl;
  
  file.close();
}
