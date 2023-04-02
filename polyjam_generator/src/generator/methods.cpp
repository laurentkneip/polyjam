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


#include <polyjam/generator/methods.hpp>
#include <sstream>
#include <fstream>

polyjam::generator::CMatrix
polyjam::generator::methods::experiment(
    const std::list<core::Poly*> & polynomials,
    const std::vector<core::Monomial> & expanders,
    bool visualization,
    bool consolePrint )
{
  bool preVisualization = false;
  CMatrix pe_matrix(polynomials);
  if(visualization && preVisualization)
    pe_matrix.visualize();

  pe_matrix.reduce();
  if(visualization && preVisualization)
    pe_matrix.visualize();
  
  if( consolePrint ) {
    std::cout << "pre-elimination is done" << std::endl;
    std::cout << "the initial matrix size is: " << pe_matrix.rows() << "x" << pe_matrix.cols() << std::endl;
  }
  
  std::list<core::Poly*> newPolynomials = pe_matrix.getPolynomials();
  CMatrix::eqs_t equations = methods::transformExpanders(expanders,newPolynomials.size());
  CMatrix matrix(newPolynomials,equations);
  
  std::list<core::Poly*>::iterator polyIter = newPolynomials.begin();
  while( polyIter != newPolynomials.end() ) {
    delete *polyIter;
    ++polyIter;
  }
  
  if( consolePrint ) {
    std::cout << "starting elimination" << std::endl;
    std::cout << "the matrix size is: " << matrix.rows() << "x" << matrix.cols() << std::endl;
  }
  
  if(visualization)
    matrix.visualize();
  matrix.reduce();
  if(visualization)
    matrix.visualize();
  
  if( consolePrint ) {
    std::cout << "elimination is done" << std::endl;
    std::cout << "the matrix size is: " << matrix.rows() << "x" << matrix.cols() << std::endl;
  }
  
  return matrix;
}

int
polyjam::generator::methods::automaticDegreeFinder(
    const std::list<core::Poly*> & polynomials,
    const std::vector<core::Monomial> & expanders,
    const std::vector<core::Monomial> & baseMonomials,
    const core::Monomial & multiplier,
    bool visualization,
    bool consolePrint,
    bool evenOnly ) {
  std::vector<core::Monomial> leadingMonomials;

  for( size_t i = 0; i < baseMonomials.size(); i++ ) {
    core::Monomial temp = baseMonomials[i] * multiplier;
    int index = -1;
    for( size_t j = 0; j < baseMonomials.size(); j++ ) {
      if( baseMonomials[j] == temp ) {
        index = j;
        break;
      }
    }
    if( index < 0 )
      leadingMonomials.push_back(temp);
  }
  
  if( consolePrint ) {
    std::cout << "The leading Monomials that we are interested in are:" << std::endl;
    for( size_t i = 0; i < leadingMonomials.size(); i++ )
      std::cout << leadingMonomials[i].getString(false) << std::endl;
  }
  
  int expanderDegree = 1;
  if(evenOnly)
    expanderDegree = 0;
  bool allFound = false;
  while( !allFound ) {
    allFound = true;
    if(evenOnly)
      expanderDegree += 2;
    else
      expanderDegree++;
    if(consolePrint)
      std::cout << "Trying out degree " << expanderDegree << std::endl;
    
    std::vector<core::Monomial> currentExpanders;
    if(evenOnly) {
      methods::generateEvendegreeExpanders(expanders,currentExpanders,expanderDegree);
    } else {
      currentExpanders = expanders;
      methods::generateSuperlinearExpanders(currentExpanders,expanderDegree);
    }
    CMatrix attempt = experiment(polynomials,currentExpanders,visualization,false);
    if(consolePrint)
      std::cout << "Template size: " << attempt.rows() << "x" << attempt.cols() << std::endl;
    CMatrix::polynomials_t polys = attempt.getPolynomials();
    
    //Now check if all polynomials have been found correctly
    for( size_t i = 0; i < leadingMonomials.size(); i++ ) {
      core::Monomial & leadingMonomial = leadingMonomials[i];
      bool found = false;
      
      std::list<core::Poly*>::iterator polysIterator = polys.begin();
      
      while( polysIterator != polys.end() ) {
        if( leadingMonomial == (*polysIterator)->leadingTerm().monomial() ) {
          bool allOtherContained = true;
          core::Poly::terms_t::iterator monoIter = (*polysIterator)->begin();
          monoIter++;
          while( monoIter != (*polysIterator)->end() )
          {
            bool contained = false;
            
            //check this monomial
            for( size_t j = 0; j < baseMonomials.size(); j++ ) {
              if( monoIter->monomial() == baseMonomials[j] ) {
                contained = true;
                break;
              }
            }
            
            if( !contained ) {
              allOtherContained = false;
              break;
            }
              
            monoIter++;
          }
          
          if( allOtherContained ) {
            found = true;
            break;
          }
        }
        
        polysIterator++;
      }
      
      if( !found ) {
        std::cout << "Did not find all monomials." << std::endl;
        allFound = false;
        break;
      }
    }
    
    auto polysIterator = polys.begin();
    while( polysIterator != polys.end() ) {
      delete (*polysIterator);
      polysIterator++;
    }
  }
  
  std::cout << "Found all monomials." << std::endl;
  return expanderDegree;
}

void
polyjam::generator::methods::generate(
    const std::list<core::Poly*> & polynomials,
    const std::list<core::Poly*> & symPolynomials,
    const std::vector<core::Monomial> & expanders,
    const std::vector<core::Monomial> & baseMonomials,
    const core::Monomial & multiplier,
    const std::string & headerFile,
    const std::string & codeFile,
    const std::string & solverName,
    const std::string & parameters,
    const std::string & save_path,
    bool visualize )
{
  /////////////////////////////
  //general configuration for additional saved data
  /////////////////////////////
  bool saveBeforePreElimination = false; // this is useful for a hand-crafted pre-elimination
  bool saveAfterPreElimination = false;  // this is useful for a hand-crafted pre-elimination
  bool saveBeforeElimination = false;    // this is useful for a hand-crafted pre-elimination
  bool saveAfterElimination = false;     // this is useful if the some manual code is required because Action matrix won't work
  bool savePolynomials = false;          // this is useful if the some manual code is required because Action matrix won't work
  bool saveMonomials = false;            // this is useful if the some manual code is required because Action matrix won't work

  bool useGaussJordan = false;
  /////////////////////////////
  /////////////////////////////
  /////////////////////////////

  std::stringstream code;
  
  //setup the actual pre-elimination matrix
  CMatrix pe_matrix(polynomials);
   
  //setup the pre-elimination matrix in the code
  CMatrix pe_helper(symPolynomials);
  int M1rows = pe_helper.rows();
  int M1cols = pe_helper.cols();

  std::stringstream M1type;
  M1type << "Eigen::MatrixXd";
  code << M1type.str() << " M1(" << M1rows << "," << M1cols << ");" << std::endl;
  code << "M1.fill(0.0);" << std::endl;

  for( int r = 0; r < M1rows; r++ ) {
    for( int c = 0; c < M1cols; c++ ) {
      if( !pe_helper(r,c).isZero() )
        code << "M1(" << r << "," << c << ") = " << pe_helper(r,c).getString(true) << "; ";
    }
    code << std::endl;
  }
  code << std::endl;

  // save the pre-elimination as pictures
  if( visualize )
    pe_matrix.visualize();
  if( saveBeforePreElimination )
    pe_matrix.save( std::string("M1before"), save_path );

  pe_matrix.reduce();
  
  if( visualize )
    pe_matrix.visualize();
  if( saveAfterPreElimination )
    pe_matrix.save( std::string("M1after"), save_path );

  CMatrix::polynomials_t zp_polynomials  = pe_matrix.getPolynomials();
  //CMatrix::polynomials_t sym_polynomials = pe_matrix.getSymbolicPolynomials( std::string("M1") );
  CMatrix::polynomials_t sym_polynomials = pe_matrix.getSymbolicPolynomials2();
  
  if(useGaussJordan) {
    code << "polyjam::gaussReduction(M1);" << std::endl << std::endl;
  } else {
    //Gauss reduction may be not stable, due to the numerics it may do something else than what we expect

    //start by extracting the reordering indices
    CMatrix::monomials_t monomials = pe_matrix.monomials();
    std::vector<int> shufflingIndices;

    //start by the leading ones
    auto it = zp_polynomials.begin();
    while( it != zp_polynomials.end() ) {
      core::Monomial leadingMonomial = (*it)->leadingMonomial().monomial();

      auto it2 = monomials.begin();
      int counter = 0;
      while( it2 != monomials.end() ) {
        if( (*it2) == leadingMonomial ) {
          shufflingIndices.push_back(counter);
          break;
        }
        counter++;
        it2++;
      }

      it++;
    }

    //Now follow up with the remaining ones
    for( int i = 0; i < M1cols; i++ ) {
      bool found = false;
      for( int j = 0; j < shufflingIndices.size(); j++ ) {
        if( i == shufflingIndices[j] ) {
          found = true;
          break;
        }
      }
      if(!found)
        shufflingIndices.push_back(i);
    }

    //Add the code for the reordering of the matrix
    code << "//swap the columns to have the leading monomials in the front" << std::endl;
    code << M1type.str() << " M1temp(" << M1rows << "," << M1cols << ");" << std::endl;
    for( int i = 0; i < shufflingIndices.size(); i++ )
      code << "M1temp.col(" << i << ") = M1.col(" << shufflingIndices[i] << ");" << std::endl;
    code << std::endl;

    //Add the actual elimination
    if( zp_polynomials.size() < M1rows ) {
      //this is the complicated case of over determination

      int M1rows2 = zp_polynomials.size();
      code << "//Use the pseudo-inverse to replace the Gauss-Jordan elimination (this here is the over determined case)" << std::endl;
      code << "Eigen::Matrix<double," << M1rows2 << "," << M1cols << "> P = M1temp.block<" << M1rows << "," << M1rows2 << ">(0,0).transpose() * M1temp;" << std::endl;
      code << "Eigen::PartialPivLU< Eigen::Matrix<double," << M1rows2 << "," << M1rows2 << "> > lupre(P.block<" << M1rows2 << "," << M1rows2 << ">(0,0));" << std::endl;
      code << "M1temp.block<" << M1rows2 << "," << (M1cols - M1rows2) << ">(0," << M1rows2 << ") = lupre.solve(P.block<" << M1rows2 << "," << (M1cols - M1rows2) << ">(0," << M1rows2 << "));" << std::endl;
      code << "M1temp.block<" << M1rows2 << "," << M1rows2 << ">(0,0) = Eigen::MatrixXd::Identity(" << M1rows2 << "," << M1rows2 << ");" << std::endl;
      code << "M1temp.block<" << (M1rows - M1rows2) << "," << M1cols << ">(" << M1rows2 << ",0) = Eigen::MatrixXd::Zero(" << (M1rows - M1rows2) << "," << M1cols << ");" << std::endl;
      code << std::endl;

    } else {

      //This here is the old version, which only works if the system is not over-determined
      code << "Eigen::Matrix<double," << M1rows << "," << M1rows << "> temp = M1temp.topLeftCorner(" << M1rows << "," << M1rows << ").inverse();\n";
      code << "Eigen::Matrix<double," << M1rows << "," << M1cols << "> temp2 = temp * M1temp;\n";
      code << "M1temp = temp2;\n\n";

    }

    //Now add the code for swapping back
    code << "//Swap the columns back to the original form" << std::endl;
    for( int i = 0; i < shufflingIndices.size(); i++ )
      code << "M1.col(" << shufflingIndices[i] << ") = M1temp.col(" << i << ");" << std::endl;
    code << std::endl;
  }

  std::cout << "Pre-elimination is done." << std::endl;

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

  std::cout << "Extracted the polynomials that are needed for composing the Action matrix." << std::endl;
  
  //ok, now we have the big matrix (plus the monomials), the polynomials that should remain (goodPolynomials),
  //plus a list of the origin of the equations
  //the goal is now to continuously remove polynomials such that all original equations remain
  std::list<int> usedEquations;
  for( size_t i = 0; i < equations.size(); i++ )
    usedEquations.push_back(i);
  
  bool removedSome = true;
  while(removedSome)
  {
  
  std::cout << "Trying to remove equations that are unnecessary." << std::endl;
  int originalNumber = usedEquations.size();
  std::list<int>::iterator ueIt = usedEquations.begin();
  int ueInd = 0;
  int toRemove = 1;
  
  while(ueIt != usedEquations.end())
  {
    std::cout << usedEquations.size() << " .. " << std::flush;
    //std::cout << "Current size of expanders is " << usedEquations.size() << ". Original size was ";
    //std::cout << originalNumber << ". Trying to remove " << toRemove << " expanders at index " << ueInd << "." << std::endl;
  
    //Remove a couple of Monomials
    std::vector<int> removed;
    for( int i = 0; i < toRemove; i++ )
    {
      removed.push_back(*ueIt);
      ueIt = usedEquations.erase(ueIt);
      if( ueIt == usedEquations.end() )
        break;
    }
    
    //std::cout << "Removed " << removed.size() << " expanders. Now computing the polynomials." << std::endl;
    
    //now check if all required polynomials are still around
    //copy the corresponding rows, and perform gaussReduction
    CMatrix subMatrix = big_matrix.subMatrix(usedEquations);
    subMatrix.reduce();
    
    if( !subMatrix.contains(goodPolynomials) )
    {
      //std::cout << "I did not find all polynomials. This trial was unsuccessful." << std::endl;
      //std::cout << "Readding the removed expanders." << std::endl;
      
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
      //std::cout << "I found all polynomials. I am increasing the size of polynomials to remove." << std::endl;
      toRemove *= 2;
      while( toRemove > (int) usedEquations.size() )
      {
        //std::cout << "not possible, need to decrease less!" << std::endl;
        toRemove /= 2;
      }
    }
  }
  
  std::cout << std::endl;
  std::cout << "I am done with this round. Original height of template was " << originalNumber << ". Now it is " << usedEquations.size() << "." << std::endl;
  
  if( true )//usedEquations.size() >= originalNumber )
    removedSome = false;
  
  }
  
  //verify that the final matrix does not change in size anymore!
  //in any case, this can be enforced (vanishing equations are simply redundant)
  CMatrix subMatrix = big_matrix.subMatrix(usedEquations);
  if(visualize)
    subMatrix.visualize();
  subMatrix.reduce();
  if(visualize)
    subMatrix.visualize();
  
  std::cout << "Removing unused monomials." << std::endl;

  //ok, now we have to redo the computation such that we use only the final equations and monomials
  std::vector<std::pair<int,core::Monomial> > finalEquations;
  for( std::list<int>::iterator i = usedEquations.begin(); i != usedEquations.end(); i++ )
    finalEquations.push_back(equations[*i]);
  
  CMatrix final_matrix(zp_polynomials,finalEquations);
  if(visualize)
    final_matrix.visualize();
  final_matrix.reduce();
  if(visualize)
    final_matrix.visualize();

  std::cout << "Final template size: " << final_matrix.rows() << "x" << final_matrix.cols() << std::endl;
  
  //once we are done with that, we should identify the leading monomials, and reorder the monomials
  std::cout << "Reordering the monomials." << std::endl;
  CMatrix::polynomials_t finalPolynomials = final_matrix.getPolynomials();
  
  CMatrix::monomials_t finalMonomials;
  auto p = finalPolynomials.begin();
  while( p != finalPolynomials.end() ) {
    finalMonomials.push_back((**p).leadingTerm().monomial());
    ++p;
  }
  CMatrix::monomials_t intMonomials = final_matrix.monomials();
  for( size_t i = 0; i < intMonomials.size(); i++ ) {
    core::Monomial temp = intMonomials[i];
    bool found = false;
    for( size_t j = 0; j < finalMonomials.size(); j++ ) {
      if( temp == finalMonomials[j] ) {
        found = true;
        break;
      }
    }
    if(!found)
      finalMonomials.push_back(temp);
  }

  //New block: reorder the equations such that we are closest possible to row-echelon form
  //finalEquations is of form: std::vector< std::pair<int,core::Monomial> >
  //use test_matrix_temp to reorder the stuff
  CMatrix test_matrix_temp( zp_polynomials, finalMonomials, finalEquations );
  std::vector<int> preIndices; preIndices.reserve(test_matrix_temp.rows());
  std::vector<int> postIndices; postIndices.reserve(test_matrix_temp.rows());
  for( int i = 0; i < test_matrix_temp.rows(); i++ )
    preIndices.push_back(i);
  int currentCol = 0;
  while( !preIndices.empty() ) {
    std::vector<int>::iterator nonzeroSearcher = preIndices.begin();
    while( nonzeroSearcher != preIndices.end() )
    {
      if( !test_matrix_temp(*nonzeroSearcher,currentCol).isZero() )
      {
        postIndices.push_back(*nonzeroSearcher);
        nonzeroSearcher = preIndices.erase(nonzeroSearcher);
      }
      else
        nonzeroSearcher++;
    }
    currentCol++;
  }

  std::vector<std::pair<int,core::Monomial> > finalReorderedEquations; finalReorderedEquations.reserve(finalEquations.size());
  for( int i = 0; i < postIndices.size(); i++ )
    finalReorderedEquations.push_back( finalEquations[postIndices[i]] );

  //verify that the reordered matrix gives a good result
  CMatrix test_matrix( zp_polynomials, finalMonomials, finalReorderedEquations );
  if(visualize)
    test_matrix.visualize();
  if(saveBeforeElimination)
    test_matrix.save( std::string("M2before"), save_path );
  test_matrix.reduce();
  if(visualize)
    test_matrix.visualize();
  if(saveAfterElimination)
    test_matrix.save( std::string("M2after"), save_path );

  if(savePolynomials) {
    //use this block here in order to print the polynomials into a file
    std::ofstream polynomialsFile;
    polynomialsFile.open( (save_path + std::string("polynomials.txt")).c_str() );
    polynomialsFile << "Printing all polynomials into the console" << std::endl;
    CMatrix::polynomials_t printPolys = test_matrix.getPolynomials();
    auto ititit = printPolys.begin();
    while( ititit != printPolys.end() ) {
      polynomialsFile << (*ititit)->getString(false);
      ititit++;
    }
    polynomialsFile << std::endl;
    polynomialsFile.close();
  }

  if(saveMonomials) {
    //use this block here to print all the monomials into a file
    std::ofstream monomialsFile;
    monomialsFile.open( (save_path + std::string("monomials.txt")).c_str() );
    monomialsFile << "Printing all monomials into the console" << std::endl;
    for( int i = 0; i < finalMonomials.size(); i++ )
    {
      monomialsFile << finalMonomials[i].getString(false);
      monomialsFile << std::endl;
    }
    monomialsFile << std::endl;
    monomialsFile.close();
  }
  
  std::cout << "Extracting the code" << std::endl;

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
  M2type << "Eigen::MatrixXd ";
  code << M2type.str() << " M2(" << M2rows << "," << M2cols << ");" << std::endl;
  code << "M2.fill(0.0);" << std::endl;
  CMatrix helper(sym_polynomials,finalMonomials,finalReorderedEquations);

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
    code << "initRow( M2, M1, " << r << ", " << finalReorderedEquations[r].first << ", ind_2_" << r << ", ind_1_" << r << ", " << numberCoefficients << "  );" << std::endl;
  }
  code << std::endl;
  
  //add the matrix inversion and multiplication
  code << "Eigen::PartialPivLU<Eigen::MatrixXd> lu(M2.block(0,0," << M2rows << "," << M2rows << "));" << std::endl;
  code << "Eigen::MatrixXd M3 = lu.solve(M2.block(0," << M2rows << "," << M2rows << "," << M3cols << "));" << std::endl;
  
  //now, add the action matrix extraction
  int solNbr = baseMonomials.size();
  
  int unknownNbr = baseMonomials[0].dimensions();
  std::stringstream solutionsType;
  solutionsType << "std::vector< Eigen::Matrix<double," << unknownNbr << ",1> >";
  std::stringstream Actiontype;
  Actiontype << "Eigen::Matrix<double," << solNbr << "," << solNbr << ">";
  code << Actiontype.str() << " Action = " << Actiontype.str() << "::Zero();" << std::endl;
  
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
      
      code << "Action.row(" << i << ") -= M3.block(" << index << "," << (M3cols - solNbr) << ",1," << solNbr << ");" << std::endl;
    }
  }
  
  //and finally add a comment block that describes the meaning of the individual columns of Action
  code << "//columns of Action mean:" << std::endl << "//";
  for( int i = 0; i < solNbr; i++ )
    code << " " << baseMonomials[i].getString(false);
  
  //finally print the code to path
  std::ofstream file;
  file.open(codeFile.c_str());
  std::ofstream header;
  header.open(headerFile.c_str());

  //ATTENTION: DO NOT CHANGE THE FOLLOWING 50 LINES. ANY DERIVATION OF THIS POLYJAM DISTRIBUTION (INCLUDING GENERATED CODE)
  //HAS TO AGREE WITH THE TERMS OF THE GNU GPL LICENSE, AND THE DISCLAIMER THEREFORE NEEDS TO BE COPIED TO THE GENERATED SOLVER
  file << "/*************************************************************************" << std::endl;
  file << " *                                                                       *" << std::endl;
  file << " * polyjam, a polynomial solver generator for C++                        *" << std::endl;
  file << " * Copyright (C) 2015 Laurent Kneip, The Australian National University  *" << std::endl;
  file << " *                                                                       *" << std::endl;
  file << " * This program is free software: you can redistribute it and/or modify  *" << std::endl;
  file << " * it under the terms of the GNU General Public License as published by  *" << std::endl;
  file << " * the Free Software Foundation, either version 3 of the License, or     *" << std::endl;
  file << " * (at your option) any later version.                                   *" << std::endl;
  file << " *                                                                       *" << std::endl;
  file << " * This program is distributed in the hope that it will be useful,       *" << std::endl;
  file << " * but WITHOUT ANY WARRANTY; without even the implied warranty of        *" << std::endl;
  file << " * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *" << std::endl;
  file << " * GNU General Public License for more details.                          *" << std::endl;
  file << " *                                                                       *" << std::endl;
  file << " * You should have received a copy of the GNU General Public License     *" << std::endl;
  file << " * along with this program.  If not, see <http://www.gnu.org/licenses/>. *" << std::endl;
  file << " *                                                                       *" << std::endl;
  file << " *************************************************************************/" << std::endl;
  file << std::endl;
  file << "//This code is automatically generated by polyjam for solving " << solverName << "." << std::endl;
  file << "//It is licensed under the GNU GPL terms." << std::endl;
  file << "//Please contact the author of polyjam for proprietary use." << std::endl;

  header << "/*************************************************************************" << std::endl;
  header << " *                                                                       *" << std::endl;
  header << " * polyjam, a polynomial solver generator for C++                        *" << std::endl;
  header << " * Copyright (C) 2015 Laurent Kneip, The Australian National University  *" << std::endl;
  header << " *                                                                       *" << std::endl;
  header << " * This program is free software: you can redistribute it and/or modify  *" << std::endl;
  header << " * it under the terms of the GNU General Public License as published by  *" << std::endl;
  header << " * the Free Software Foundation, either version 3 of the License, or     *" << std::endl;
  header << " * (at your option) any later version.                                   *" << std::endl;
  header << " *                                                                       *" << std::endl;
  header << " * This program is distributed in the hope that it will be useful,       *" << std::endl;
  header << " * but WITHOUT ANY WARRANTY; without even the implied warranty of        *" << std::endl;
  header << " * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *" << std::endl;
  header << " * GNU General Public License for more details.                          *" << std::endl;
  header << " *                                                                       *" << std::endl;
  header << " * You should have received a copy of the GNU General Public License     *" << std::endl;
  header << " * along with this program.  If not, see <http://www.gnu.org/licenses/>. *" << std::endl;
  header << " *                                                                       *" << std::endl;
  header << " *************************************************************************/" << std::endl;
  header << std::endl;
  header << "//This code is automatically generated by polyjam for solving " << solverName << "." << std::endl;
  header << "//It is licensed under the GNU GPL terms." << std::endl;
  header << "//Please contact the author of polyjam for proprietary use." << std::endl;

  file << std::endl;
  file << "#include \"" << solverName << ".hpp\"" << std::endl;
  file << std::endl;
  file << std::endl;
  if( useGaussJordan ) {
    file << "void" << std::endl;
    file << "polyjam::gaussReduction( Eigen::MatrixXd & matrix ) {" << std::endl;
    file << std::endl;
    file << "  //create some constant for zero-checking etc." << std::endl;
    file << "  double precision = 0.0000000001;" << std::endl;
    file << "  " << std::endl;
    file << "  int rows = matrix.rows();" << std::endl;
    file << "  int cols = matrix.cols();" << std::endl;
    file << "  int maxIndentation = cols-1;" << std::endl;
    file << "  " << std::endl;
    file << "  //first compute the number of steps to do" << std::endl;
    file << "  int steps = rows;" << std::endl;
    file << "  if( rows > cols )" << std::endl;
    file << "    steps = cols;" << std::endl;
    file << "  " << std::endl;
    file << "  //working variables" << std::endl;
    file << "  int currentIndentation = 0;" << std::endl;
    file << "  int frontRow = 0;" << std::endl;
    file << "  " << std::endl;
    file << "  //first step down" << std::endl;
    file << "  for( int dummy = 0; dummy < steps; dummy++ )" << std::endl;
    file << "  {" << std::endl;
    file << "    // first iterate through the rows and find a row that has a" << std::endl;
    file << "    // non-zero coefficient in the relevant column" << std::endl;
    file << "    " << std::endl;
    file << "    int row;" << std::endl;
    file << "    " << std::endl;
    file << "    while(currentIndentation <= maxIndentation)" << std::endl;
    file << "    {" << std::endl;
    file << "      double maxValue = -1.0;" << std::endl;
    file << "      " << std::endl;
    file << "      for( int tempRow = frontRow; tempRow < rows; tempRow++ )" << std::endl;
    file << "      {" << std::endl;
    file << "        double value = fabs(matrix(tempRow,currentIndentation));" << std::endl;
    file << "        if( value > maxValue )" << std::endl;
    file << "        {" << std::endl;
    file << "          row = tempRow;" << std::endl;
    file << "          maxValue = value;" << std::endl;
    file << "        }" << std::endl;
    file << "      }" << std::endl;
    file << "      " << std::endl;
    file << "      //if, after checking this column, we found a non-zero element, we are done" << std::endl;
    file << "      if(maxValue > precision)" << std::endl;
    file << "        break;" << std::endl;
    file << "      " << std::endl;
    file << "      //if not, we have to move on to the next column" << std::endl;
    file << "      currentIndentation++;" << std::endl;
    file << "    }" << std::endl;
    file << "    " << std::endl;
    file << "    //if we are beyond the maxIndentation, the whole rest is zero. Break!" << std::endl;
    file << "    if(currentIndentation > maxIndentation)" << std::endl;
    file << "      break;" << std::endl;
    file << "    " << std::endl;
    file << "    //row is now the row that should go in the place of frontRow->swap" << std::endl;
    file << "    Eigen::MatrixXd rowCopy = matrix.row(row);" << std::endl;
    file << "    matrix.row(row) = matrix.row(frontRow);" << std::endl;
    file << "    matrix.row(frontRow) = rowCopy;" << std::endl;
    file << "    " << std::endl;
    file << "    //ok, now use frontRow!" << std::endl;
    file << "    int col = currentIndentation;" << std::endl;
    file << "    double leadingCoefficient = matrix(frontRow,col);" << std::endl;
    file << "    matrix.row(frontRow) /= leadingCoefficient;" << std::endl;
    file << "    " << std::endl;
    file << "    //iterate through all remaining rows, and subtract correct multiple of " << std::endl;
    file << "    //first row (if leading coefficient is non-zero!)" << std::endl;
    file << "    row = frontRow;" << std::endl;
    file << "    ++row;" << std::endl;
    file << "    while( row < rows )" << std::endl;
    file << "    {" << std::endl;
    file << "      double leadingCoefficient = matrix(row,currentIndentation);" << std::endl;
    file << "      if( fabs(leadingCoefficient) > precision )" << std::endl;
    file << "        matrix.row(row) -= leadingCoefficient * matrix.row(frontRow);" << std::endl;
    file << "      ++row;" << std::endl;
    file << "    }" << std::endl;
    file << "    " << std::endl;
    file << "    //increment row and currentIndentation" << std::endl;
    file << "    ++frontRow;" << std::endl;
    file << "    currentIndentation++;" << std::endl;
    file << "    if( currentIndentation > maxIndentation )" << std::endl;
    file << "      break;" << std::endl;
    file << "  }" << std::endl;
    file << "  " << std::endl;
    file << "  //if the following, this is the zero matrix -> nothing to do!" << std::endl;
    file << "  if( frontRow == 0 )" << std::endl;
    file << "    return;" << std::endl;
    file << "  " << std::endl;
    file << "  //actually, everything starting from frontRow should be zero -> delete" << std::endl;
    file << "  Eigen::MatrixXd croppedMatrix = matrix.block(0,0,frontRow,cols);" << std::endl;
    file << "  matrix = croppedMatrix;" << std::endl;
    file << "  " << std::endl;
    file << "  //set index to the last non-zero row" << std::endl;
    file << "  --frontRow;" << std::endl;
    file << "  " << std::endl;
    file << "  //Now step up" << std::endl;
    file << "  while( frontRow != 0 )" << std::endl;
    file << "  {" << std::endl;  
    file << "    //indent until we find a non-zero element in frontRow" << std::endl;
    file << "    int indentations = 0;" << std::endl;
    file << "    while( fabs(matrix(frontRow,indentations)) < precision && indentations < cols )" << std::endl;
    file << "      indentations++;" << std::endl;
    file << "    " << std::endl;
    file << "    //if the following, there is a problem, maybe zero-matrix!" << std::endl;
    file << "    if( indentations == cols )" << std::endl;
    file << "      break;" << std::endl;
    file << "    " << std::endl;
    file << "    //get the working row" << std::endl;
    file << "    int row = frontRow;" << std::endl;
    file << "    " << std::endl;
    file << "    do" << std::endl;
    file << "    {" << std::endl;
    file << "      //decrement working row" << std::endl;
    file << "      --row;" << std::endl;
    file << "      " << std::endl;
    file << "      //working column" << std::endl;
    file << "      " << std::endl;
    file << "      //now get the leading coefficient" << std::endl;
    file << "      double leadingCoefficient = matrix(row,indentations);" << std::endl;
    file << "      " << std::endl;
    file << "      //Now iterator until the end, and subtract each time the multiplied" << std::endl;
    file << "      //front-row" << std::endl;
    file << "      if( fabs(leadingCoefficient) > precision )" << std::endl;
    file << "        matrix.row(row) -= leadingCoefficient * matrix.row(frontRow);" << std::endl;
    file << "    }" << std::endl;
    file << "    while( row != 0 );" << std::endl;
    file << "    " << std::endl;
    file << "    --frontRow;" << std::endl;
    file << "  }" << std::endl;
    file << "}" << std::endl;
    file << std::endl;
    file << std::endl;
  }
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
  file << "    M2(row2,cols2[i]) = M1(row1,cols1[i]);" << std::endl;
  file << "}" << std::endl;
  file << std::endl;
  file << "void" << std::endl;
  file << "polyjam::" << solverName << "::solve( " << parameters << ", " << solutionsType.str() << " & solutions )" << std::endl;
  file << "{" << std::endl;
  file << code.str() << std::endl;
  file << std::endl;
  file << "  Eigen::EigenSolver< Eigen::Matrix<double," << solNbr << "," << solNbr << "> > Eig(Action,true);" << std::endl;
  file << "  Eigen::Matrix<std::complex<double>," << solNbr << ",1> D = Eig.eigenvalues();" << std::endl;
  file << "  Eigen::Matrix<std::complex<double>," << solNbr << "," << solNbr << "> V = Eig.eigenvectors();" << std::endl;
  file << std::endl;
  file << "  for( int c = 0; c < " << solNbr << "; c++ )" << std::endl;
  file << "  {" << std::endl;
  file << "    std::complex<double> eigValue = D[c];" << std::endl;
  file << std::endl;
  file << "    if( fabs(eigValue.imag()) < 0.0001 )" << std::endl;
  file << "    {" << std::endl;
  file << "      Eigen::Matrix<double," << unknownNbr << ",1> sol;" << std::endl;
  file << std::endl;
  file << "      std::complex<double> temp;" << std::endl;

  for( int d = 0; d < unknownNbr; d++ )
  {
    //search for the index of a certain unknown in the base-monomials
    std::vector<unsigned int> exponents( unknownNbr, 0 );
    exponents[d] = 1;
    core::Monomial unknown( exponents );

    int b = -1;
    for( b = 0; b < baseMonomials.size(); b++ )
    {
      if( unknown == baseMonomials[b] )
        break;
    }

  file << "      temp = V(" << b << ",c) / V(" << solNbr - 1 << ",c);" << std::endl;
  file << "      sol(" << d << "," << 0 << ") = temp.real();" << std::endl;

  }

  file << "      solutions.push_back(sol);" << std::endl;
  file << "    }" << std::endl;
  file << "  }" << std::endl;
  file << "}";
  file.close();

  header << std::endl;
  header << "#ifndef POLYJAM_" << solverName << "_HPP_" << std::endl;
  header << "#define POLYJAM_" << solverName << "_HPP_" << std::endl;
  header << std::endl;
  header << "#include <stdlib.h>" << std::endl;
  header << "#include <Eigen/Eigen>" << std::endl;
  header << "#include <vector>" << std::endl;
  header << "#include <list>" << std::endl;
  header << std::endl;
  header << std::endl;
  header << "namespace polyjam" << std::endl;
  header << "{" << std::endl;

  if( useGaussJordan ) {
    header << std::endl;
    header << "void gaussReduction( Eigen::MatrixXd & matrix );" << std::endl;
    header << std::endl;
  }

  header << "namespace " << solverName << std::endl;
  header << "{" << std::endl;
  header << std::endl;
  header << "  void initRow(" << std::endl;
  header << "      " << M2type.str() << " & M2," << std::endl;
  header << "      const " << M1type.str() << " & M1," << std::endl;
  header << "      int row2," << std::endl;
  header << "      int row1," << std::endl;
  header << "      const int * cols2," << std::endl;
  header << "      const int * cols1," << std::endl;
  header << "      size_t numberCols );" << std::endl;
  header << std::endl;
  header << "  void solve( " << parameters << ", " << solutionsType.str() << " & solutions );" << std::endl;
  header << std::endl;
  header << "}" << std::endl;
  header << "}" << std::endl;
  header << std::endl;
  header << "#endif /* POLYJAM_" << solverName << "_HPP_ */";
  header.close();
}

polyjam::generator::CMatrix::eqs_t
polyjam::generator::methods::transformExpanders(
    const std::vector<core::Monomial> & expanders, size_t polynomials )
{
  CMatrix::eqs_t equations;
  
  for( int row = 0; row < (int) polynomials; row++ )
    equations.push_back(std::pair<int,core::Monomial>( row,expanders.front().dimensions() ));
  
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
polyjam::generator::methods::generateSuperlinearExpanders( std::vector<core::Monomial> & expanders, int maxDegree )
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
polyjam::generator::methods::generateEvendegreeExpanders(
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
