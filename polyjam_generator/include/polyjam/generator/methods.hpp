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
 
#ifndef POLYJAM_GENERATOR_METHODS_HPP_
#define POLYJAM_GENERATOR_METHODS_HPP_

#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <list>

#include <polyjam/core/Poly.hpp>
#include <polyjam/generator/CMatrix.hpp>


/**
 * \brief The namespace of this library.
 */
namespace polyjam
{

/**
 * \brief The namespace for the elimination templates.
 */
namespace generator
{

/**
 * \brief The namespace of all the methods related to elimination and code
 *        extraction.
 */
namespace methods
{

CMatrix experiment(
    const std::list<core::Poly*> & polynomials,
    const std::vector<core::Monomial> & expanders,
    bool visualization = false,
    bool consolePrint = false );

int automaticDegreeFinder(
    const std::list<core::Poly*> & polynomials,
    const std::vector<core::Monomial> & expanders,
    const std::vector<core::Monomial> & baseMonomials,
    const core::Monomial & multiplier,
    bool visualization = false,
    bool consolePrint = false,
    bool evenOnly = false );
    
void generate(
    const std::list<core::Poly*> & polynomials,
    const std::list<core::Poly*> & symPolynomials,
    const std::vector<core::Monomial> & expanders,
    const std::vector<core::Monomial> & baseMonomials,
    const core::Monomial & multiplier,
    const std::string & headerFile,
    const std::string & codeFile,
    const std::string & solverName,
    const std::string & parameters,
    bool visualize = false );

void generateSuperlinearExpanders( std::vector<core::Monomial> & expanders, int maxDegree );

void generateEvendegreeExpanders( const std::vector<core::Monomial> & originalMonomials, std::vector<core::Monomial> & expanders, int maxDegree );

CMatrix::eqs_t transformExpanders(
    const std::vector<core::Monomial> & expanders, size_t polynomials );
    
}
}
}

#endif /* POLYJAM_GENERATOR_METHODS_HPP_ */
