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

#include <polyjam/exports/MacaulayTwo.hpp>
#include <polyjam/fields/Zp.hpp>

using namespace std;

polyjam::exports::MacaulayTwo::MacaulayTwo()
{
  _numberVariables = 0;
  _characteristic = 0;
  _numberEquations = 0;
}

polyjam::exports::MacaulayTwo::~MacaulayTwo()
{
}

void
polyjam::exports::MacaulayTwo::addPoly( const core::Poly & poly )
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
polyjam::exports::MacaulayTwo::write( const std::string & fileName )
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
