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


//Todo (potential problems that I discovered:
//-The symbolic form of the polynomials may be different from the other one (even in finite prime field).
// This is for instance the case if some symbolic expression cancels out.
// In a sense the problem is worst if an entire column in the initial matrix M1 becomes zero.
// It may be not the case in the symbolic form, and the two matrices become different (different dimension)
//-The gauss elimination variant is now missing, and M2 is possibly too large

#include <polyjam/elimination/methods.hpp>
#include <sstream>
#include <fstream>

polyjam::elimination::CMatrix
polyjam::elimination::methods::experiment(
    const std::list<core::Poly*> & polynomials,
    const std::vector<core::Monomial> & expanders,
    bool visualization,
    bool consolePrint )
{
  CMatrix pe_matrix(polynomials);
  if(visualization)
    pe_matrix.visualize(true);

  pe_matrix.reduce();
  if(visualization)
    pe_matrix.visualize(true);
  
  if( consolePrint )
  {
    std::cout << "pre-elimination is done" << std::endl;
    std::cout << "the initial matrix size is: " << pe_matrix.rows() << "x" << pe_matrix.cols() << std::endl;
  }
  
  std::list<core::Poly*> newPolynomials = pe_matrix.getPolynomials();
  CMatrix::eqs_t equations = methods::transformExpanders(expanders,newPolynomials.size());
  CMatrix matrix(newPolynomials,equations);
  
  std::list<core::Poly*>::iterator polyIter = newPolynomials.begin();
  while( polyIter != newPolynomials.end() )
  {
    delete *polyIter;
    ++polyIter;
  }
  
  if( consolePrint )
  {
    std::cout << "starting elimination" << std::endl;
    std::cout << "the matrix size is: " << matrix.rows() << "x" << matrix.cols() << std::endl;
  }
  
  if(visualization)
    matrix.visualize();
  matrix.reduce();
  if(visualization)
    matrix.visualize();
  
  if( consolePrint )
  {
    std::cout << "elimination is done" << std::endl;
    std::cout << "the matrix size is: " << matrix.rows() << "x" << matrix.cols() << std::endl;
  }
  
  return matrix;
}

int
polyjam::elimination::methods::automaticDegreeFinder(
    const std::list<core::Poly*> & polynomials,
    const std::vector<core::Monomial> & expanders,
    const std::vector<core::Monomial> & baseMonomials,
    const core::Monomial & multiplier,
    bool visualization,
    bool consolePrint,
    bool evenOnly )
{
  std::vector<core::Monomial> leadingMonomials;

  for( size_t i = 0; i < baseMonomials.size(); i++ )
  {
    core::Monomial temp = baseMonomials[i] * multiplier;
    int index = -1;
    for( size_t j = 0; j < baseMonomials.size(); j++ )
    {
      if( baseMonomials[j] == temp )
      {
        index = j;
        break;
      }
    }
    if( index < 0 )
      leadingMonomials.push_back(temp);
  }
  
  if( consolePrint )
  {
    std::cout << "The leading Monomials that we are interested in are:" << std::endl;
    for( size_t i = 0; i < leadingMonomials.size(); i++ )
      std::cout << leadingMonomials[i].getString(false) << std::endl;
  }
  
  int expanderDegree = 1;
  if(evenOnly)
    expanderDegree = 0;
  bool allFound = false;
  while( !allFound )
  {
    allFound = true;
    if(evenOnly)
      expanderDegree += 2;
    else
      expanderDegree++;
    if(consolePrint)
      std::cout << "Trying out degree " << expanderDegree << std::endl;
    
    std::vector<core::Monomial> currentExpanders;
    if(evenOnly)
    {
      methods::generateEvendegreeExpanders(expanders,currentExpanders,expanderDegree);
    }
    else
    {
      currentExpanders = expanders;
      methods::generateSuperlinearExpanders(currentExpanders,expanderDegree);
    }
    CMatrix attempt = experiment(polynomials,currentExpanders,visualization,consolePrint);
    std::list<core::Poly*> polys = attempt.getPolynomials();
    
    //Now check if all polynomials have been found correctly
    for( size_t i = 0; i < leadingMonomials.size(); i++ )
    {
      core::Monomial & leadingMonomial = leadingMonomials[i];
      bool found = false;
      
      std::list<core::Poly*>::iterator polysIterator = polys.begin();
      
      while( polysIterator != polys.end() )
      {
        if( leadingMonomial == (*polysIterator)->leadingTerm().monomial() )
        {
          std::cout << "found a leading monomial" << std::endl;
          bool allOtherContained = true;
          core::Poly::terms_t::iterator monoIter = (*polysIterator)->begin();
          monoIter++;
          while( monoIter != (*polysIterator)->end() )
          {
            bool contained = false;
            
            //check this monomial
            for( size_t j = 0; j < baseMonomials.size(); j++ )
            {
              if( monoIter->monomial() == baseMonomials[j] )
              {
                contained = true;
                break;
              }
            }
            
            if( !contained )
            {
              allOtherContained = false;
              break;
            }
              
            monoIter++;
          }
          
          if( allOtherContained )
          {
            found = true;
            break;
          }
        }
        
        polysIterator++;
      }
      
      if( !found )
      {
        allFound = false;
        break;
      }
    }
    
    std::list<core::Poly*>::iterator polysIterator = polys.begin();
    while( polysIterator != polys.end() )
    {
      delete (*polysIterator);
      polysIterator++;
    }
  }
  
  return expanderDegree;
}

void
polyjam::elimination::methods::generate(
    const std::list<core::Poly*> & polynomials,
    const std::list<core::Poly*> & symPolynomials,
    const std::vector<core::Monomial> & expanders,
    const std::vector<core::Monomial> & baseMonomials,
    const core::Monomial & multiplier,
    const std::string & path,
    const std::string & solverName,
    const std::string & parameters,
    bool visualize )
{
  std::stringstream code;
  
  //setup the actual pre-elimination matrix
  CMatrix pe_matrix(polynomials);
   
  //setup the pre-elimination matrix in the code
  CMatrix pe_helper(symPolynomials);
  int M1rows = pe_helper.rows();
  int M1cols = pe_helper.cols();

  std::stringstream M1type;
  M1type << "std::vector<std::vector<double>*>";
  code << M1type.str() << " M1;" << std::endl;
  code << "for( int r = 0; r < " << M1rows << "; r++ )" << std::endl;
  code << "  M1.push_back( new std::vector<double>(" << M1cols << ",0.0) );" << std::endl;

  for( int r = 0; r < M1rows; r++ )
  {
    for( int c = 0; c < M1cols; c++ )
    {
      if( !pe_helper(r,c).isZero() )
        code << "(*(M1[" << r << "]))[" << c << "] = " << pe_helper(r,c).getString(true) << "; ";
    }
    code << std::endl;
  }
  code << std::endl;
  
  //Do the pre-elimination
  pe_matrix.reduce();
  std::list<core::Poly*> zp_polynomials  = pe_matrix.getPolynomials();
  //std::list<core::Poly*> sym_polynomials = pe_matrix.getSymbolicPolynomials( std::string("M1") );
  std::list<core::Poly*> sym_polynomials = pe_matrix.getSymbolicPolynomials2();
  code << "polyjam::elimination::gaussReduction(M1);" << std::endl << std::endl;

  //Now transform the vector of expanders and create the big matrix
  CMatrix::eqs_t equations = transformExpanders( expanders, zp_polynomials.size() );
  CMatrix big_matrix(zp_polynomials,equations);
  
  //extract the don't miss Polys automatically
  CMatrix attempt(zp_polynomials,equations);
  attempt.reduce();
  std::list<core::Poly*> goodPolynomials;
  for( size_t i = 0; i < baseMonomials.size(); i++ )
  {
    core::Monomial multipliedBase = baseMonomials[i] * multiplier;
    //check if we can find the multiplied base in the base
    bool inBase = false;
    for( size_t j = 0; j < baseMonomials.size(); j++ )
    {
      if( baseMonomials[j] == multipliedBase )
      {
        inBase = true;
        break;
      }
    }
    if( !inBase )
    {
      //ok, extract this polynomial
      for( size_t j = 0; j < attempt.rows(); j++ )
      {
        core::Poly* tempPoly = new core::Poly(attempt.getPolynomial(j));
        if( tempPoly->leadingTerm().monomial() == multipliedBase )
        {
          //we really need to find all of them here!
          goodPolynomials.push_back(tempPoly);
        }
        else
        {
          delete tempPoly;
        }
      }
    }
  }
  
  //ok, now we have the big matrix (plus the monomials), the polynomials that should remain (goodPolynomials),
  //plus a list of the origin of the equations
  //the goal is now to continuously remove polynomials such that all original equations remain
  std::list<int> usedEquations;
  for( size_t i = 0; i < equations.size(); i++ )
    usedEquations.push_back(i);
  
  bool removedSome = true;
  while(removedSome)
  {
  
  std::cout << "Starting with new round of computation of minimal equations." << std::endl;
  int originalNumber = usedEquations.size();
  std::list<int>::iterator ueIt = usedEquations.begin();
  int ueInd = 0;
  int toRemove = 1;
  
  while(ueIt != usedEquations.end())
  {
    std::cout << "Current size of expanders is " << usedEquations.size() << ". Original size was ";
    std::cout << originalNumber << ". Trying to remove " << toRemove << " expanders at index " << ueInd << "." << std::endl;
  
    //Remove a couple of Monomials
    std::vector<int> removed;
    for( int i = 0; i < toRemove; i++ )
    {
      removed.push_back(*ueIt);
      ueIt = usedEquations.erase(ueIt);
      if( ueIt == usedEquations.end() )
        break;
    }
    
    std::cout << "Removed " << removed.size() << " expanders. Now computing the polynomials." << std::endl;
    
    //now check if all required polynomials are still around
    //copy the correct rows of _origMatrix, and perform gaussReduction
    CMatrix subMatrix = big_matrix.subMatrix(usedEquations);
    subMatrix.reduce();
    
    if( !subMatrix.contains(goodPolynomials) )
    {
      std::cout << "I did not find all polynomials. This trial was unsuccessful." << std::endl;
      std::cout << "Readding the removed expanders." << std::endl;
      
      for( int i = removed.size()-1; i >= 0; i-- )
        ueIt = usedEquations.insert(ueIt,removed[i]);
      
      if( removed.size() > 1 )
        toRemove /= 2;
      else
      {
        ueIt++; ueInd++;
      }
    }
    else
    {
      std::cout << "I found all polynomials. I am increasing the size of polynomials to remove." << std::endl;
      toRemove *= 2;
      while( toRemove > (int) usedEquations.size() )
      {
        std::cout << "not possible, need to decrease less!" << std::endl;
        toRemove /= 2;
      }
    }
  }
  
  std::cout << "I am done with this round. Original size of expanders was " << originalNumber << ". Now it is " << usedEquations.size() << "." << std::endl;
  
  if( true )//usedEquations.size() >= originalNumber )
    removedSome = false;
  
  }
  
  //verify that the final matrix does not change in size anymore!
  //in any case, this can be enforced (vanishing equations are simply redundant)
  std::cout << "I am done with the whole reduction. The original number of equations is " << usedEquations.size() << "." << std::endl;
  CMatrix subMatrix = big_matrix.subMatrix(usedEquations);
  if(visualize)
    subMatrix.visualize();
  subMatrix.reduce();
  if(visualize)
    subMatrix.visualize();  
  std::cout << "The number of equations after the reduction is: " << subMatrix.rows() << std::endl;
  
  //ok, now we have to redo the computation such that _origMonomials contains only the present monomials
  std::vector<std::pair<int,core::Monomial> > finalEquations;
  for( std::list<int>::iterator i = usedEquations.begin(); i != usedEquations.end(); i++ )
    finalEquations.push_back(equations[*i]);
  
  CMatrix final_matrix(zp_polynomials,finalEquations);
  if(visualize)
    final_matrix.visualize();
  final_matrix.reduce();
  if(visualize)
    final_matrix.visualize();
  
  //once we are done with that, we should identify the leading monomials, and reorder the monomials
  std::list<core::Poly*> finalPolynomials = final_matrix.getPolynomials();
  
  std::vector<core::Monomial> finalMonomials;
  std::list<core::Poly*>::iterator p = finalPolynomials.begin();
  while( p != finalPolynomials.end() )
  {
    finalMonomials.push_back((**p).leadingTerm().monomial());
    ++p;
  }
  std::vector<core::Monomial> intMonomials = final_matrix.monomials();
  for( size_t i = 0; i < intMonomials.size(); i++ )
  {
    core::Monomial temp = intMonomials[i];
    bool found = false;
    for( size_t j = 0; j < finalMonomials.size(); j++ )
    {
      if( temp == finalMonomials[j] )
      {
        found = true;
        break;
      }
    }
    if(!found)
      finalMonomials.push_back(temp);
  }

  //verify that the reordered matrix gives a good result
  CMatrix test_matrix( zp_polynomials, finalMonomials, finalEquations );
  if(visualize)
    test_matrix.visualize();
  test_matrix.reduce();
  if(visualize)
    test_matrix.visualize();

  /*
  std::cout << "starting print of monomials" << std::endl;
  for( int i = 0; i < finalMonomials.size(); i++ )
  {
    finalMonomials[i].print();
    std::cout << std::endl;
  }
  */
  
  //everything alright
  //summary:
  //we have the size of the elimination matrix from final_matrix
  //we have the monomial order on top (finalMonomials)
  //we have a list of the finalEquations, containing index of initial equation, plus monomial expander
  //we have a list of the initial polynomials in symbolic form (sym_polynomials)
  
  //now extract the code
  
  int M2rows = final_matrix.rows();
  int M2cols = final_matrix.cols();
  int M3cols = M2cols - M2rows;

  std::stringstream M2type;
  M2type << "Eigen::Matrix<double," << M2rows << "," << M2cols << ">";
  code << M2type.str() << " M2 = " << M2type.str() <<  "::Zero();" << std::endl;
  CMatrix helper(sym_polynomials,finalMonomials,finalEquations);

  for( int r = 0; r < M2rows; r++ )
  {
    code << "static const int ind_2_" << r << " [] = {";
    bool insertComma = false;
    for( int c = 0; c < M2cols; c++ )
    {
      if( !helper(r,c).isZero() )
      {
        if( insertComma )
          code << ",";
        code << c;
        insertComma = true;
      }
    }
    code << "};" << std::endl;

    size_t numberCoefficients = 0;

    code << "static const int ind_1_" << r << " [] = {";
    insertComma = false;
    for( int c = 0; c < M2cols; c++ )
    {
      if( !helper(r,c).isZero() )
      {
        if( insertComma )
          code << ",";
        core::Coefficient temp = helper(r,c) - helper(r,c).one();
        code << temp.getString();
        insertComma = true;
        numberCoefficients++;
      }
    }
    code << "};" << std::endl;
    code << "initRow( M2, M1, " << r << ", " << finalEquations[r].first << ", ind_2_" << r << ", ind_1_" << r << ", " << numberCoefficients << "  );" << std::endl;
  }
  code << std::endl;
  
  //add the matrix inversion and multiplication
  code << "Eigen::PartialPivLU< Eigen::Matrix<double," << M2rows << "," << M2rows << "> > lu(M2.block<" << M2rows << "," << M2rows << ">(0,0));" << std::endl;
  code << "Eigen::Matrix<double," << M2rows << "," << M3cols << "> M3 = ";
  code << "lu.solve(M2.block<" << M2rows << "," << M3cols << ">(0," << M2rows << "));" << std::endl;
  
  //now, add the action matrix extraction
  int solNbr = baseMonomials.size();
  
  std::stringstream Actiontype;
  Actiontype << "Eigen::Matrix<double," << solNbr << "," << solNbr << ">";
  code << "Action = " << Actiontype.str() << "::Zero();" << std::endl;
  
  for( int i = 0; i < solNbr; i++ )
  {
    core::Monomial temp = baseMonomials[i] * multiplier;
    int index = -1;
    for( int j = 0; j < solNbr; j++ )
    {
      if( baseMonomials[j] == temp )
      {
        index = j;
        break;
      }
    }
    if( index >= 0 )
      code << "Action(" << i << "," << index << ") = 1.0;" << std::endl;
    else
    {
      //get the values from the correct equation in M3
      for( size_t j = 0; j < finalMonomials.size(); j++ )
      {
        if( temp == finalMonomials[j] )
        {
          index = j;
          break;
        }
      }
      
      code << "Action.row(" << i << ") -= M3.block<1," << solNbr << ">(" << index << "," << (M3cols - solNbr) << ");" << std::endl;
    }
  }
  
  //and finally add a comment block that describes the meaning of the individual columns of Action
  code << "//columns of Action mean:" << std::endl << "//";
  for( int i = 0; i < solNbr; i++ )
    code << " " << baseMonomials[i].getString(false);
  
  //finally print the code to path
  std::ofstream file;
  file.open(path.c_str());
  file << "/******************************************************************************" << std::endl;
  file << " * Author:   Laurent Kneip                                                    *" << std::endl;
  file << " * Contact:  kneip.laurent@gmail.com                                          *" << std::endl;
  file << " * License:  Copyright (c) 2013 Laurent Kneip, ANU. All rights reserved.      *" << std::endl;
  file << " *                                                                            *" << std::endl;
  file << " * Redistribution and use in source and binary forms, with or without         *" << std::endl;
  file << " * modification, are permitted provided that the following conditions         *" << std::endl;
  file << " * are met:                                                                   *" << std::endl;
  file << " * * Redistributions of source code must retain the above copyright           *" << std::endl;
  file << " *   notice, this list of conditions and the following disclaimer.            *" << std::endl;
  file << " * * Redistributions in binary form must reproduce the above copyright        *" << std::endl;
  file << " *   notice, this list of conditions and the following disclaimer in the      *" << std::endl;
  file << " *   documentation and/or other materials provided with the distribution.     *" << std::endl;
  file << " * * Neither the name of ANU nor the names of its contributors may be         *" << std::endl;
  file << " *   used to endorse or promote products derived from this software without   *" << std::endl;
  file << " *   specific prior written permission.                                       *" << std::endl;
  file << " *                                                                            *" << std::endl;
  file << " * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS \"AS IS\"*" << std::endl;
  file << " * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE  *" << std::endl;
  file << " * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE *" << std::endl;
  file << " * ARE DISCLAIMED. IN NO EVENT SHALL ANU OR THE CONTRIBUTORS BE LIABLE        *" << std::endl;
  file << " * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL *" << std::endl;
  file << " * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR *" << std::endl;
  file << " * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER *" << std::endl;
  file << " * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT         *" << std::endl;
  file << " * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY  *" << std::endl;
  file << " * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF     *" << std::endl;
  file << " * SUCH DAMAGE.                                                               *" << std::endl;
  file << " ******************************************************************************/" << std::endl;
  file << std::endl;
  file << "//This code is automatically generated by polyjam for solving " << solverName << std::endl;
  file << std::endl;
  file << "#include <stdlib.h>" << std::endl;
  file << "#include <Eigen/Eigen>" << std::endl;
  file << "#include <polyjam/elimination/GaussJordan.hpp>" << std::endl;
  file << std::endl;
  file << std::endl;
  file << "void" << std::endl;
  file << "polyjam::" << solverName << "::initRow(" << std::endl;
  file << "    " << M2type.str() << " & M2," << std::endl;
  file << "    const " << M1type.str() << " & M1," << std::endl;
  file << "    int row2," << std::endl;
  file << "    int row1," << std::endl;
  file << "    const int * cols2," << std::endl;
  file << "    const int * cols1," << std::endl;
  file << "    size_t numberCols )" << std::endl;
  file << "{" << std::endl;
  file << "  for( int i = 0; i < numberCols; i++ )" << std::endl;
  file << "    M2(row2,cols2[i]) = (*(M1[row1]))[cols1[i]];" << std::endl;
  file << "}" << std::endl;
  file << std::endl;
  file << "void" << std::endl;
  file << "polyjam::" << solverName << "::solve( " << parameters << ", " << Actiontype.str() << " & Action )" << std::endl;
  file << "{" << std::endl;
  file << code.str() << std::endl;
  file << "}";
  file.close();
}

polyjam::elimination::CMatrix::eqs_t
polyjam::elimination::methods::transformExpanders(
    const std::vector<core::Monomial> & expanders, size_t polynomials )
{
  CMatrix::eqs_t equations;
  
  for( int row = 0; row < (int) polynomials; row++ )
    equations.push_back(std::pair<int,core::Monomial>(row,expanders.front().dimensions()));
  
  std::vector<core::Monomial>::const_iterator expIter = expanders.begin();
  while( expIter != expanders.end() )
  {    
    CMatrix::eqs_t::iterator firstEq = equations.begin();
    for( size_t polyInd = 0; polyInd < polynomials; polyInd++ )
    {
      firstEq = equations.insert( firstEq, std::pair<int,core::Monomial>(polyInd,(*expIter) ) );
      ++firstEq;
    }
    ++expIter;
  }
  
  return equations;
}

void
polyjam::elimination::methods::generateSuperlinearExpanders( std::vector<core::Monomial> & expanders, int maxDegree )
{
  int numberLinearExpanders = expanders.size();
  
  //quadratic expanders
  if( maxDegree >= 2 )
  {
    for( int i = 0; i < numberLinearExpanders; i++ )
    {
      for( int j = i; j < numberLinearExpanders; j++ )
        expanders.push_back( expanders[i]*expanders[j] );
    }
  }
  
  //third-order expanders
  if( maxDegree >= 3 )
  {
    for( int i = 0; i < numberLinearExpanders; i++ )
    {
      for( int j = i; j < numberLinearExpanders; j++ )
      {
        for( int k = j; k < numberLinearExpanders; k++ )
          expanders.push_back(expanders[i]*expanders[j]*expanders[k]);
      }
    }
  }
  
  //fourth-order expanders
  if( maxDegree >= 4 )
  {
    for( int i = 0; i < numberLinearExpanders; i++ )
    {
      for( int j = i; j < numberLinearExpanders; j++ )
      {
        for( int k = j; k < numberLinearExpanders; k++ )
        {
          for( int l = k; l < numberLinearExpanders; l++ )
            expanders.push_back(expanders[i]*expanders[j]*expanders[k]*expanders[l]);
        }
      }
    }
  }
    
  //fifth-order expanders
  if( maxDegree >= 5 )
  {
    for( int i = 0; i < numberLinearExpanders; i++ )
    {
      for( int j = i; j < numberLinearExpanders; j++ )
      {
        for( int k = j; k < numberLinearExpanders; k++ )
        {
          for( int l = k; l < numberLinearExpanders; l++ )
          {
            for( int m = l; m < numberLinearExpanders; m++ )
              expanders.push_back(expanders[i]*expanders[j]*expanders[k]*expanders[l]*expanders[m]);
          }
        }
      }
    }
  }
    
  //sixth-order expanders
  if( maxDegree >= 6 )
  {
    for( int i = 0; i < numberLinearExpanders; i++ )
    {
      for( int j = i; j < numberLinearExpanders; j++ )
      {
        for( int k = j; k < numberLinearExpanders; k++ )
        {
          for( int l = k; l < numberLinearExpanders; l++ )
          {
            for( int m = l; m < numberLinearExpanders; m++ )
            {
              for( int n = m; n < numberLinearExpanders; n++ )
                expanders.push_back(expanders[i]*expanders[j]*expanders[k]*expanders[l]*expanders[m]*expanders[n]);
            }
          }
        }
      }
    }
  }
    
  //seventh-order expanders
  if( maxDegree >= 7 )
  {
    for( int i = 0; i < numberLinearExpanders; i++ )
    {
      for( int j = i; j < numberLinearExpanders; j++ )
      {
        for( int k = j; k < numberLinearExpanders; k++ )
        {
          for( int l = k; l < numberLinearExpanders; l++ )
          {
            for( int m = l; m < numberLinearExpanders; m++ )
            {
              for( int n = m; n < numberLinearExpanders; n++ )
              {
                for( int o = n; o < numberLinearExpanders; o++ )
                  expanders.push_back(expanders[i]*expanders[j]*expanders[k]*expanders[l]*expanders[m]*expanders[n]*expanders[o]);
              }
            }
          }
        }
      }
    }
  }
    
  //eigth-order expanders
  if( maxDegree >= 8 )
  {
    for( int i = 0; i < numberLinearExpanders; i++ )
    {
      for( int j = i; j < numberLinearExpanders; j++ )
      {
        for( int k = j; k < numberLinearExpanders; k++ )
        {
          for( int l = k; l < numberLinearExpanders; l++ )
          {
            for( int m = l; m < numberLinearExpanders; m++ )
            {
              for( int n = m; n < numberLinearExpanders; n++ )
              {
                for( int o = n; o < numberLinearExpanders; o++ )
                {
                  for( int p = o; p < numberLinearExpanders; p++ )
                    expanders.push_back(expanders[i]*expanders[j]*expanders[k]*expanders[l]*expanders[m]*expanders[n]*expanders[o]*expanders[p]);
                }
              }
            }
          }
        }
      }
    }
  }
    
  //ninth-order expanders
  if( maxDegree >= 9 )
  {
    for( int i = 0; i < numberLinearExpanders; i++ )
    {
      for( int j = i; j < numberLinearExpanders; j++ )
      {
        for( int k = j; k < numberLinearExpanders; k++ )
        {
          for( int l = k; l < numberLinearExpanders; l++ )
          {
            for( int m = l; m < numberLinearExpanders; m++ )
            {
              for( int n = m; n < numberLinearExpanders; n++ )
              {
                for( int o = n; o < numberLinearExpanders; o++ )
                {
                  for( int p = o; p < numberLinearExpanders; p++ )
                  {
                    for( int q = p; q < numberLinearExpanders; q++ )
                      expanders.push_back(expanders[i]*expanders[j]*expanders[k]*expanders[l]*expanders[m]*expanders[n]*expanders[o]*expanders[p]*expanders[q]);
                  }
                }
              }
            }
          }
        }
      }
    }
  }
    
  //tenth-order expanders
  if( maxDegree >= 10 )
  {
    for( int i = 0; i < numberLinearExpanders; i++ )
    {
      for( int j = i; j < numberLinearExpanders; j++ )
      {
        for( int k = j; k < numberLinearExpanders; k++ )
        {
          for( int l = k; l < numberLinearExpanders; l++ )
          {
            for( int m = l; m < numberLinearExpanders; m++ )
            {
              for( int n = m; n < numberLinearExpanders; n++ )
              {
                for( int o = n; o < numberLinearExpanders; o++ )
                {
                  for( int p = o; p < numberLinearExpanders; p++ )
                  {
                    for( int q = p; q < numberLinearExpanders; q++ )
                    {
                      for( int r = q; r < numberLinearExpanders; r++ )
                        expanders.push_back(expanders[i]*expanders[j]*expanders[k]*expanders[l]*expanders[m]*expanders[n]*expanders[o]*expanders[p]*expanders[q]*expanders[r]);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  //eleventh-order expanders
  if( maxDegree >= 11 )
  {
    for( int i = 0; i < numberLinearExpanders; i++ )
    {
      for( int j = i; j < numberLinearExpanders; j++ )
      {
        for( int k = j; k < numberLinearExpanders; k++ )
        {
          for( int l = k; l < numberLinearExpanders; l++ )
          {
            for( int m = l; m < numberLinearExpanders; m++ )
            {
              for( int n = m; n < numberLinearExpanders; n++ )
              {
                for( int o = n; o < numberLinearExpanders; o++ )
                {
                  for( int p = o; p < numberLinearExpanders; p++ )
                  {
                    for( int q = p; q < numberLinearExpanders; q++ )
                    {
                      for( int r = q; r < numberLinearExpanders; r++ )
                      {
                        for( int s = r; s < numberLinearExpanders; s++ )
                          expanders.push_back(expanders[i]*expanders[j]*expanders[k]*expanders[l]*expanders[m]*expanders[n]*expanders[o]*expanders[p]*expanders[q]*expanders[r]*expanders[s]);
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  //twelveth-order expanders
  if( maxDegree >= 12 )
  {
    for( int i = 0; i < numberLinearExpanders; i++ )
    {
      for( int j = i; j < numberLinearExpanders; j++ )
      {
        for( int k = j; k < numberLinearExpanders; k++ )
        {
          for( int l = k; l < numberLinearExpanders; l++ )
          {
            for( int m = l; m < numberLinearExpanders; m++ )
            {
              for( int n = m; n < numberLinearExpanders; n++ )
              {
                for( int o = n; o < numberLinearExpanders; o++ )
                {
                  for( int p = o; p < numberLinearExpanders; p++ )
                  {
                    for( int q = p; q < numberLinearExpanders; q++ )
                    {
                      for( int r = q; r < numberLinearExpanders; r++ )
                      {
                        for( int s = r; s < numberLinearExpanders; s++ )
                        {
                          for( int t = s; t < numberLinearExpanders; t++ )
                            expanders.push_back(expanders[i]*expanders[j]*expanders[k]*expanders[l]*expanders[m]*expanders[n]*expanders[o]*expanders[p]*expanders[q]*expanders[r]*expanders[s]*expanders[t]);
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

void
polyjam::elimination::methods::generateEvendegreeExpanders(
    const std::vector<core::Monomial> & originalMonomials,
    std::vector<core::Monomial> & expanders,
    int maxDegree )
{
  int numberLinearExpanders = originalMonomials.size();
  
  //quadratic expanders
  if( maxDegree >= 2 )
  {
    for( int i = 0; i < numberLinearExpanders; i++ )
    {
      for( int j = i; j < numberLinearExpanders; j++ )
        expanders.push_back( originalMonomials[i]*originalMonomials[j] );
    }
  }
    
  //fourth-order expanders
  if( maxDegree >= 4 )
  {
    for( int i = 0; i < numberLinearExpanders; i++ )
    {
      for( int j = i; j < numberLinearExpanders; j++ )
      {
        for( int k = j; k < numberLinearExpanders; k++ )
        {
          for( int l = k; l < numberLinearExpanders; l++ )
            expanders.push_back(originalMonomials[i]*originalMonomials[j]*originalMonomials[k]*originalMonomials[l]);
        }
      }
    }
  }
    
  //sixth-order expanders
  if( maxDegree >= 6 )
  {
    for( int i = 0; i < numberLinearExpanders; i++ )
    {
      for( int j = i; j < numberLinearExpanders; j++ )
      {
        for( int k = j; k < numberLinearExpanders; k++ )
        {
          for( int l = k; l < numberLinearExpanders; l++ )
          {
            for( int m = l; m < numberLinearExpanders; m++ )
            {
              for( int n = m; n < numberLinearExpanders; n++ )
                expanders.push_back(originalMonomials[i]*originalMonomials[j]*originalMonomials[k]*originalMonomials[l]*originalMonomials[m]*originalMonomials[n]);
            }
          }
        }
      }
    }
  }
    
  //eigth-order expanders
  if( maxDegree >= 8 )
  {
    for( int i = 0; i < numberLinearExpanders; i++ )
    {
      for( int j = i; j < numberLinearExpanders; j++ )
      {
        for( int k = j; k < numberLinearExpanders; k++ )
        {
          for( int l = k; l < numberLinearExpanders; l++ )
          {
            for( int m = l; m < numberLinearExpanders; m++ )
            {
              for( int n = m; n < numberLinearExpanders; n++ )
              {
                for( int o = n; o < numberLinearExpanders; o++ )
                {
                  for( int p = o; p < numberLinearExpanders; p++ )
                    expanders.push_back(originalMonomials[i]*originalMonomials[j]*originalMonomials[k]*originalMonomials[l]*originalMonomials[m]*originalMonomials[n]*originalMonomials[o]*originalMonomials[p]);
                }
              }
            }
          }
        }
      }
    }
  }
    
  //tenth-order expanders
  if( maxDegree >= 10 )
  {
    for( int i = 0; i < numberLinearExpanders; i++ )
    {
      for( int j = i; j < numberLinearExpanders; j++ )
      {
        for( int k = j; k < numberLinearExpanders; k++ )
        {
          for( int l = k; l < numberLinearExpanders; l++ )
          {
            for( int m = l; m < numberLinearExpanders; m++ )
            {
              for( int n = m; n < numberLinearExpanders; n++ )
              {
                for( int o = n; o < numberLinearExpanders; o++ )
                {
                  for( int p = o; p < numberLinearExpanders; p++ )
                  {
                    for( int q = p; q < numberLinearExpanders; q++ )
                    {
                      for( int r = q; r < numberLinearExpanders; r++ )
                        expanders.push_back(originalMonomials[i]*originalMonomials[j]*originalMonomials[k]*originalMonomials[l]*originalMonomials[m]*originalMonomials[n]*originalMonomials[o]*originalMonomials[p]*originalMonomials[q]*originalMonomials[r]);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}
