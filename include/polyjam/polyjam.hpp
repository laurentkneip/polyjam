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

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <sstream>

#include <polyjam/core/Poly.hpp>
#include <polyjam/core/PolyMatrix.hpp>
#include <polyjam/elimination/methods.hpp>
#include <polyjam/exports/MacaulayTwo.hpp>

using namespace std;
using namespace polyjam;
using namespace polyjam::fields;
using namespace polyjam::core;
using namespace polyjam::elimination;
using namespace polyjam::exports;

namespace polyjam
{
  string exec( const char * cmd );
  void extractMonomials( const std::string & input, vector<Monomial> & baseMonomials, size_t nu );
  void initGenerator();
  void execGenerator( list<Poly*> & eqs, list<Poly*> & eqs_sym, const string & solverName, const string & parameters, bool visualize = false );
  void execGenerator( list<Poly*> & eqs, const string & solverName, const string & parameters, bool visualize = false );
  vector<int> findOccurences(string input, char c);
  Monomial stringToMonomial( const std::string & input, size_t nu );
  vector<string> splitString( const std::string & input, const vector<int> & indices );
}
