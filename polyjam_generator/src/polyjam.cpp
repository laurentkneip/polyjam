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

#include <polyjam/polyjam.hpp>
#include <sys/time.h>
#include <cstdio>
#include <memory>
#include <sys/types.h>
#include <sys/stat.h>


void
polyjam::initGenerator()
{
  std::cout << std::endl;
  std::cout << "polyjam" << std::endl;
  std::cout << "Copyright (C) 2015 Laurent Kneip, The Australian National University" << std::endl;
  std::cout << "This program comes with ABSOLUTELY NO WARRANTY; It is free software:" << std::endl;
  std::cout << "you can redistribute it and/or modify it under the terms of the GNU General Public License" << std::endl;
  std::cout << std::endl;

  struct timeval tic;
  gettimeofday( &tic, 0 );
  srand ( tic.tv_usec );
}

void
polyjam::execGenerator( list<Poly*> & eqs, const string & solverName, const string & parameters, bool visualize, const string & save_path )
{
  execGenerator(eqs, solverName, std::string(""), parameters, visualize, save_path);
}

void
polyjam::execGenerator( list<Poly*> & eqs, const string & solverName, const string & suffix, const string & parameters, bool visualize, const string & save_path )
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
  execGenerator( eqs_zp, eqs_sym, solverName, suffix, parameters, visualize, save_path );
}

void
polyjam::execGenerator( list<Poly*> & eqs, list<Poly*> & eqs_sym, const string & solverName, const string & parameters, bool visualize, const string & save_path )
{
  execGenerator(eqs, eqs_sym, solverName, std::string(""), parameters, visualize, save_path);
}

void
polyjam::execGenerator( list<Poly*> & eqs, list<Poly*> & eqs_sym, const string & solverName, const string & suffix, const string & parameters, bool visualize, const string & save_path )
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

  std::cout << "Analysing the Groebner basis in Macaulay2 ..." << std::endl;

  //export the Zp equations to a Macaulay2 script
  stringstream subdir;
  subdir << WORKSPACEPATH << solverName << "/M2script";

  struct stat info;
  if( stat( subdir.str().c_str(), &info ) != 0 )
  {
    stringstream dircmd;
    dircmd << "mkdir " << subdir.str();
    system(dircmd.str().c_str());
  }

  stringstream tempfile;
  if(suffix.empty()){
    tempfile << subdir.str() << "/" << solverName << ".m2";
  }
  else{
    tempfile << subdir.str() << "/" << solverName << "_" << suffix << ".m2";
  }
  ExportMacaulay exportMacaulay;
  
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

  std::cout << "The degree of the basis and the basis monomials (if existing) are:" << std::endl;
  std::cout << macaulayOutput << std::endl;

  //analyse the dimensionality
  string line;
  istringstream input_iss(macaulayOutput);
  getline(input_iss, line);
  int dim = atoi(line.c_str());
  if( dim != 0 )
  {
    if( dim > 0 )
      std::cout << "Stopping here. The dimensionality of the ideal is bigger than zero! This means that there are not enough equations, and the problem is underconstrained." << std::endl;
    
    if( dim < 0 )
      std::cout << "Stopping here. The dimensionality of the ideal is smaller than zero! This means that there are too many equations, and the problem is overconstrained." << std::endl;

    return;
  }

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

  std::cout << "Finding the degree of expansion." << std::endl;
  int expanderDegree = methods::automaticDegreeFinder(
      eqs, expanders, baseMonomials, multiplier, visualize, true );

  std::cout << "Generating the super-linear expanders." << std::endl;
  //generate all expander based on the degree, and generate the solver
  methods::generateSuperlinearExpanders(expanders,expanderDegree);

  //generate sub-directory if it does not exist
  stringstream subdir2;
  if(suffix.empty()){
    subdir2 << SOLVERPATH << solverName;
  }
  else{
    subdir2 << SOLVERPATH << solverName << "_" << suffix;
  }

  struct stat info2;
  if( stat( subdir2.str().c_str(), &info2 ) != 0 )
  {
    stringstream dircmd;
    dircmd << "mkdir " << subdir2.str();
    system(dircmd.str().c_str());
  }

  std::cout << "Starting the solver generation." << std::endl;
  stringstream codeFile;
  if(save_path != ""){
    codeFile << save_path << solverName << ".cpp";
  }
  else if(suffix.empty()){
    codeFile << SOLVERPATH << solverName << "/" << solverName << ".cpp";
  }
  else{
    codeFile << SOLVERPATH << solverName << "/" << solverName << "_" << suffix << ".cpp";
  }
  stringstream headerFile;
  if(save_path != ""){
    headerFile << save_path << solverName << ".hpp";
  }
  else if(suffix.empty()){
    headerFile << SOLVERPATH << solverName << "/" << solverName << ".hpp";
  }
  else{
    headerFile << SOLVERPATH << solverName << "/" << solverName << "_" << suffix << ".hpp";
  }

  if(suffix.empty()){
    methods::generate( eqs, eqs_sym, expanders, baseMonomials, multiplier, headerFile.str(), codeFile.str(), (solverName), parameters, visualize , save_path);
  }
  else{
    methods::generate( eqs, eqs_sym, expanders, baseMonomials, multiplier, headerFile.str(), codeFile.str(), (solverName + "_" + suffix), parameters, visualize , save_path );
  }
}

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

void
polyjam::extractMonomials( const std::string & input, vector<Monomial> & baseMonomials, size_t nu )
{
  //first split up into lines, and get rid of the ones that are not useful
  list<string> lines;
  string line;
  istringstream input_iss(input);
  //ATTENTION: Remove the first line which contains the dimensionality of the ideal
  getline(input_iss, line);
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
