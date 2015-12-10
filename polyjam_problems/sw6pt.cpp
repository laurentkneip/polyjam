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

int main( int argc, char** argv )
{  
  //initialize the random generator
  initGenerator();
  size_t nu = 3; //the number of unknowns in the problem

  //****** Part 1: get all equations with both random Zp and symbolic measurements *******//

  //initialize the input with random coefficients
  PolyMatrix F1(Poly::zeroSZ(nu),3,3);
  PolyMatrix F2(Poly::zeroSZ(nu),3,3);
  PolyMatrix F3(Poly::zeroSZ(nu),3,3);

  for( int r = 0; r < 3; r++ )
  {
    for( int c = 0; c < 3; c++ )
    {
      //we chose smart symbolic names that actually represent the location of the measurement in the final code
      //(e.g. inside a vector or-in this case-a matrix)
      stringstream name1; name1 << "F1(" << r << "," << c << ")";
      stringstream name2; name2 << "F2(" << r << "," << c << ")";
      stringstream name3; name3 << "F3(" << r << "," << c << ")";

      F1(r,c) = Poly::SrandZ(name1.str(),nu);
      F2(r,c) = Poly::SrandZ(name2.str(),nu);
      F3(r,c) = Poly::SrandZ(name3.str(),nu);
    }
  }

  //initialize unknowns
  Poly x = Poly::uSZ(1,nu);
  Poly y = Poly::uSZ(2,nu);
  Poly w = Poly::uSZ(3,nu);

  //prepare some intermediate variables
  PolyMatrix F = F1 * x + F2 * y + F3;
  PolyMatrix Ft = F.transpose();
  PolyMatrix Q(Poly::oneSZ(nu),3,3,true);
  Q(2,2) = w.clone();
  PolyMatrix FQFtQ = F*Q*Ft*Q;
  PolyMatrix te = FQFtQ * F * Poly::constSZ(2,nu) - F * FQFtQ.trace();

  //store all equations
  list<Poly*> eqs;
  eqs.push_back(new Poly(F.determinant()));
  for( int r = 0; r < 3; r++ )
  {
    for( int c = 0; c < 3; c++ )
      eqs.push_back(new Poly(te(r,c)));
  }

  //****** Part 2: Generate the solver ***************
  execGenerator( eqs, string("sw6pt"), string("Eigen::Matrix3d & F1, Eigen::Matrix3d & F2, Eigen::Matrix3d & F3") );
}
