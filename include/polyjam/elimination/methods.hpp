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
 
#ifndef POLYJAM_ELIMINATION_METHODS_HPP_
#define POLYJAM_ELIMINATION_METHODS_HPP_

#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <list>

#include <polyjam/core/Poly.hpp>
#include <polyjam/elimination/CMatrix.hpp>


/**
 * \brief The namespace of this library.
 */
namespace polyjam
{

/**
 * \brief The namespace for the elimination templates.
 */
namespace elimination
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
    const std::string & path,
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

#endif /* POLYJAM_ELIMINATION_METHODS_HPP_ */
