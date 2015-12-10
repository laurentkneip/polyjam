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

#include <polyjam/polyjam.hpp>
#include <sys/time.h>
#include <cstdio>
#include <memory>

string
polyjam::exec( const char * cmd )
{
  std::shared_ptr<FILE> pipe( popen(cmd,"r"), pclose );
  if(!pipe)
    return "ERROR";
  char buffer[128];
  std::string result = "";
  while( !feof(pipe.get()) )
  {
    if( fgets( buffer,128,pipe.get() ) != NULL )
      result += buffer;
  }
  return result;
}

vector<int>
polyjam::findOccurences(string input, char c)
{
  vector<int> occurences;
  for( size_t i = 0; i < input.size(); i++ )
    if( input[i] == c )
      occurences.push_back(i);

  return occurences;
}

vector<string>
polyjam::splitString( const std::string & input, const vector<int> & indices )
{
  vector<string> subStrings;
  for( size_t i = 0; i < (indices.size()-1); i++ )
    subStrings.push_back( input.substr( indices[i], indices[i+1]-indices[i] ) );
  subStrings.push_back( input.substr(indices.back(),input.size()-indices.back()) );
  return subStrings;
}

Monomial
polyjam::stringToMonomial( const std::string & input, size_t nu )
{
  vector<unsigned int> exponents(nu,0);

  vector<int> x_occurences = findOccurences( input, 'x' );
  
  if( !x_occurences.empty() )
  {
    vector<string> subStrings = splitString( input, x_occurences );

    for( size_t i = 0; i < subStrings.size(); i++ )
    {
      //get the index for this part of monomial
      int dimension = 0;
      
      //get the exponent and the dimension
      int exponent = 1;
      int dimensionIndexEnd = subStrings[i].size();

      vector<int> hatIndices = polyjam::findOccurences(subStrings[i], '^');
      if( !hatIndices.empty() )
      {
        int exponentIndex = hatIndices[0]+1;
        string exponentString = subStrings[i].substr( exponentIndex, subStrings[i].size() - exponentIndex );
        istringstream iss(exponentString);
        iss >> exponent;
        dimensionIndexEnd = hatIndices[0];
      }
      
      vector<int> barIndices = polyjam::findOccurences(subStrings[i], '_');
      int dimensionIndex = barIndices[0]+1;
      string dimensionString = subStrings[i].substr( dimensionIndex, dimensionIndexEnd - dimensionIndex );
      istringstream iss(dimensionString);
      iss >> dimension;

      exponents[dimension-1] = exponent;
    }
  }

  return Monomial(exponents);
}

void
polyjam::extractMonomials( const std::string & input, vector<Monomial> & baseMonomials, size_t nu )
{
  //first split up into lines, and get rid of the ones that are not useful
  list<string> lines;
  string line;
  istringstream input_iss(input);
  while( getline(input_iss, line) )
  {
    if( line.substr(0,4) != "----" )
      lines.push_back(line);
  }

  //now remove the bars from the beginning and the end
  size_t firstLength = lines.front().size();
  string newFirst = lines.front().substr(2,firstLength-2);
  lines.front() = newFirst;

  size_t lastLength = lines.back().size();
  string newLast = lines.back().substr(0,lastLength-2);
  lines.back() = newLast;

  

  //now break up each row into all monomial strings
  list<string> monomialStrings;
  list<string>::iterator it = lines.begin();
  while( it != lines.end() )
  {
    istringstream iss(*it);
    do
    {
      string monomialString;
      iss >> monomialString;
      monomialStrings.push_back(monomialString);
    } while(iss);

    //for some reason the last element is nil
    monomialStrings.pop_back();
    it++;
  }

  //print the monomial strings for now
  list<string>::iterator it2 = monomialStrings.begin();
  while( it2 != monomialStrings.end() )
  {
    Monomial monomial = stringToMonomial( *it2, nu );
    baseMonomials.push_back(monomial);
    it2++;
  }
}

void
polyjam::initGenerator()
{
  struct timeval tic;
  gettimeofday( &tic, 0 );
  srand ( tic.tv_usec );
}

void
polyjam::execGenerator( list<Poly*> & eqs, const string & solverName, const string & parameters, bool visualize )
{
  //the goal of this function is to split up the polynomials into symbolic and finite field ones

  list<Poly*> eqs_zp;
  list<Poly*> eqs_sym;

  list<Poly*>::iterator it = eqs.begin();
  while( it != eqs.end() )
  {
    //the first coefficient in each term is the Zp one
    //the second coefficient in each term is the symbolic one

    (*it)->setDominant(0); //-> activate the Zp one
    eqs_zp.push_back(new Poly( (*it)->clone(false) )); //-> clone, but only dominant part
    (*it)->setDominant(1); //-> activate the Zp one
    eqs_sym.push_back(new Poly( (*it)->clone(false) )); //-> clone, but only dominant part

    it++;
  }

  //done, call the original execGenerator
  execGenerator( eqs_zp, eqs_sym, solverName, parameters, visualize );
}

void
polyjam::execGenerator( list<Poly*> & eqs, list<Poly*> & eqs_sym, const string & solverName, const string & parameters, bool visualize )
{
  //create a list of monomials for all the unknowns (those will expand the original system of equations)
  vector<Monomial> expanders;
  int nu = eqs.front()->leadingTerm().monomial().dimensions();
  vector<unsigned int> exponents(nu,0);
  for( int i = 0; i < nu; i++ )
  {
    exponents[i] = 1;
    expanders.push_back( Monomial(exponents) );
    exponents[i] = 0;
  }

  //export the Zp equations to a Macaulay2 script
  stringstream tempfile;
  tempfile << TEMPPATH << solverName << ".m2";
  MacaulayTwo exportMacaulay;
  
  std::list<Poly*>::const_iterator eqIter = eqs.begin();
  while( eqIter != eqs.end() )
  {
    exportMacaulay.addPoly(**eqIter);
    ++eqIter;
  }
  exportMacaulay.write(tempfile.str());

  stringstream cmd;
  cmd << MACAULAYCOMMAND << " --silent " << tempfile.str();
  string macaulayOutput = exec(cmd.str().c_str());

  std::cout << macaulayOutput << std::endl;

  //now extract the basis monomials from this string
  std::vector<Monomial> baseMonomials_temp;
  extractMonomials(macaulayOutput,baseMonomials_temp,nu);

  //now sort the base monomials
  std::vector<Monomial> baseMonomials;
  Poly orderingPolynomial = Poly::zeroZ(nu);
  Poly c1 = Poly::oneZ(nu);
  for( size_t i = 0; i < baseMonomials_temp.size(); i++ )
    orderingPolynomial += Term( c1.leadingTerm().coefficient().clone(), baseMonomials_temp[i] );
  Poly::terms_t::iterator it2 = orderingPolynomial.begin();
  while( it2 != orderingPolynomial.end() )
  {
    baseMonomials.push_back( it2->monomial() );
    it2++;
  }

  //create the action variable
  vector<unsigned int> action(nu,0);
  action.back() = 1;
  Monomial multiplier(action);

  int expanderDegree = methods::automaticDegreeFinder(
      eqs, expanders, baseMonomials, multiplier, visualize, false );

  //generate all expander based on the degree, and generate the solver
  methods::generateSuperlinearExpanders(expanders,expanderDegree);
  stringstream codeFile;
  codeFile << SOLVERPATH << solverName << ".cpp";
  methods::generate( eqs, eqs_sym, expanders, baseMonomials, multiplier, codeFile.str(), solverName, parameters, visualize );
}
