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

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <sstream>

#include <polyjam/core/Poly.hpp>
#include <polyjam/core/PolyMatrix.hpp>
#include <polyjam/generator/methods.hpp>
#include <polyjam/generator/ExportMacaulay.hpp>

using namespace std;
using namespace polyjam;
using namespace polyjam::fields;
using namespace polyjam::core;
using namespace polyjam::generator;

namespace polyjam
{
  void initGenerator();
  void execGenerator( list<Poly*> & eqs, const string & solverName, const string & parameters, bool visualize = false );
  void execGenerator( list<Poly*> & eqs, list<Poly*> & eqs_sym, const string & solverName, const string & parameters, bool visualize = false );
  void execGenerator( list<Poly*> & eqs, const string & solverName, const string & suffix, const string & parameters, bool visualize = false );
  void execGenerator( list<Poly*> & eqs, list<Poly*> & eqs_sym, const string & solverName, const string & suffix, const string & parameters, bool visualize = false );

  //The following functions are for interaction with Macaulay and interpretation of the result
  string exec( const char * cmd );
  void extractMonomials( const string & input, vector<Monomial> & baseMonomials, size_t nu );
  vector<int> findOccurences( string input, char c);
  Monomial stringToMonomial( const string & input, size_t nu );
  vector<string> splitString( const string & input, const vector<int> & indices );
}
