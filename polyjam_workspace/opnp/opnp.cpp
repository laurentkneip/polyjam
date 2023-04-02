#include <polyjam/polyjam.hpp>

int main( int argc, char** argv )
{  
  //initialize the random generator
  initGenerator();
  size_t nu = 4; //the number of unknowns in the problem
  size_t numberBearingVectors = 3;

  //****** Part 1: get independent random variables *******//

  //random rotation
  std::cout << "generating random rotation" << std::endl;
  PolyMatrix q_gt(Poly::zeroZ(nu),4,1);
  for( int j = 0; j < 4; j++ )
    q_gt[j] = Poly::randZ(nu);
  
  Poly scale = (q_gt[0]*q_gt[0]) + (q_gt[1]*q_gt[1]) + (q_gt[2]*q_gt[2]) + (q_gt[3]*q_gt[3]);
  //We need to make sure that the scale is the inverse of the average depth
  //Let us calculate the total depth, it will be used later to make sure things stay consistent
  Poly averageDepth = Poly::oneZ(nu).leadingTerm() / scale.leadingTerm();
  Poly totalDepth = Poly::constZ(numberBearingVectors,nu) * averageDepth;
  
  //Derive the scaled rotation matrix
  PolyMatrix R_gt(Poly::zeroZ(nu),3,3);
  R_gt(0,0) = (q_gt[0]*q_gt[0]) + (q_gt[1]*q_gt[1]) - (q_gt[2]*q_gt[2]) - (q_gt[3]*q_gt[3]);
  R_gt(1,1) = (q_gt[0]*q_gt[0]) - (q_gt[1]*q_gt[1]) + (q_gt[2]*q_gt[2]) - (q_gt[3]*q_gt[3]);
  R_gt(2,2) = (q_gt[0]*q_gt[0]) - (q_gt[1]*q_gt[1]) - (q_gt[2]*q_gt[2]) + (q_gt[3]*q_gt[3]);
  R_gt(0,1) = Poly::constZ(2,nu) * (q_gt[1]*q_gt[2]-q_gt[3]*q_gt[0]);
  R_gt(0,2) = Poly::constZ(2,nu) * (q_gt[1]*q_gt[3]+q_gt[2]*q_gt[0]);
  R_gt(1,2) = Poly::constZ(2,nu) * (q_gt[2]*q_gt[3]-q_gt[1]*q_gt[0]);
  R_gt(1,0) = Poly::constZ(2,nu) * (q_gt[1]*q_gt[2]+q_gt[3]*q_gt[0]);
  R_gt(2,0) = Poly::constZ(2,nu) * (q_gt[1]*q_gt[3]-q_gt[2]*q_gt[0]);
  R_gt(2,1) = Poly::constZ(2,nu) * (q_gt[2]*q_gt[3]+q_gt[1]*q_gt[0]);
  
  //random translation (this is not the true translation, but a scale one again)
  std::cout << "generating random translation" << std::endl;
  PolyMatrix t_gt(Poly::zeroZ(nu),3,1);
  for( int i = 0; i < 3; i++ )
    t_gt[i] = Poly::randZ(nu);
  
  //random image points and random depths (for opnp, we need to start from the image points)
  std::cout << "generating random image points" << std::endl;
  std::vector<PolyMatrix> fs;
  std::vector<Poly> depths;
  for( int i = 0; i < (int) numberBearingVectors; i++ ) {
    fs.push_back(PolyMatrix(Poly::zeroZ(nu),3,1));
    fs.back()[0] = Poly::randZ(nu);
    fs.back()[1] = Poly::randZ(nu);
    fs.back()[2] = Poly::oneZ(nu);
    depths.push_back( Poly::randZ(nu) );
  }

  //Important trick: reset the third depth such that we obtain a consistent total depth
  depths[2] = totalDepth - depths[0] - depths[1];

  //****** Part 2: extract consistent world points wihin Zp *****//
  std::cout << "extracting the consistent world points" << std::endl;
  
  //Calculate the inverse of R_gt (again, R_gt is not an orthonormal matrix)
  Poly scaleSquared = scale * scale;
  Poly scaleSquaredInverse = Poly::oneZ(nu).leadingTerm() / scaleSquared.leadingTerm();
  PolyMatrix R_gt_inv = R_gt.transpose() * scaleSquaredInverse;

  std::vector<PolyMatrix> wps;
  for( int i = 0; i < (int) numberBearingVectors; i++ ) {
    Poly scaledDepth = depths[i].leadingTerm() / averageDepth.leadingTerm();
    wps.push_back( R_gt_inv * ( fs[i] * scaledDepth - t_gt) );
  }

  /////////////////////////////////////////////////////////////
  //////////////////generate double wps////////////////////////
  /////////////////////////////////////////////////////////////
  std::vector<PolyMatrix> wps_sz;
  for( int i = 0; i < (int) numberBearingVectors; i++ ) {
    wps_sz.push_back(PolyMatrix(Poly::zeroSZ(nu),3,1));

    for( int j = 0; j < 3; j++ ) {      
      std::stringstream name; name << "wps[" << i << "][" << j << "]";
      wps_sz.back()[j] = Poly(Term(
        Coefficient(name.str()),
        wps[i][j].leadingTerm().coefficient().clone(),
        Monomial(nu) ));
    }
  }

  /////////////////////////////////////////////////////////////
  //////////////////generate double fs/////////////////////////
  /////////////////////////////////////////////////////////////
  std::vector<PolyMatrix> fs_sz;
  for( int i = 0; i < (int) numberBearingVectors; i++ ) {
    fs_sz.push_back(PolyMatrix(Poly::zeroSZ(nu),3,1));

    for( int j = 0; j < 3; j++ ) {
      std::stringstream name; name << "fs[" << i << "][" << j << "]";
      fs_sz.back()[j] = Poly(Term(
        Coefficient(name.str()),
        fs[i][j].leadingTerm().coefficient().clone(),
        Monomial(nu) ));
    }
  }

  //******** Part 3: start the calculation of the polynomial coefficients ******//

  //1) find center of image points
  Poly n = Poly::constZ(numberBearingVectors,nu);

  PolyMatrix f_center(Poly::zeroZ(nu),3,1);
  for( int i = 0; i < numberBearingVectors; i++ )
    f_center += fs[i];
  for( int j = 0; j < 3; j++ )
    f_center[j] = f_center[j].leadingTerm()/n.leadingTerm();

  /////////////////////////////////////////////////////////////
  //////////////////generate double f_center///////////////////
  /////////////////////////////////////////////////////////////
  PolyMatrix f_center_sz(Poly::zeroSZ(nu),3,1);

  for( int j = 0; j < 3; j++ ) {
    std::stringstream name; name << "f_center[" << j << "]";
    f_center_sz[j] = Poly(Term(
        Coefficient(name.str()),
        f_center[j].leadingTerm().coefficient().clone(),
        Monomial(nu) ));
  }
  /////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////

  //2) find the center of the world points
  PolyMatrix wp_center(Poly::zeroZ(nu),3,1);
  for( int i = 0; i < numberBearingVectors; i++ )
    wp_center += wps[i];
  for( int j = 0; j < 3; j++ )
    wp_center[j] = wp_center[j].leadingTerm()/n.leadingTerm();

  /////////////////////////////////////////////////////////////
  //////////////////generate double wp_center//////////////////
  /////////////////////////////////////////////////////////////
  PolyMatrix wp_center_sz(Poly::zeroSZ(nu),3,1);

  for( int j = 0; j < 3; j++ ) {
    std::stringstream name; name << "wp_center[" << j << "]";
    wp_center_sz[j] = Poly(Term(
        Coefficient(name.str()),
        wp_center[j].leadingTerm().coefficient().clone(),
        Monomial(nu) ));
  }
  /////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////

  //2) create inner summation terms
  vector<PolyMatrix> uwps, vwps;
  PolyMatrix uwp_center(Poly::zeroZ(nu),3,1);
  PolyMatrix vwp_center(Poly::zeroZ(nu),3,1);
  for( int i = 0; i < numberBearingVectors; i++ ) {
    uwps.push_back( (wps[i]-wp_center) * fs[i][0] ); uwp_center += uwps.back();
    vwps.push_back( (wps[i]-wp_center) * fs[i][1] ); vwp_center += vwps.back();
  }
  for( int j = 0; j < 3; j++ ) {
    uwp_center[j] = uwp_center[j].leadingTerm()/n.leadingTerm();
    vwp_center[j] = vwp_center[j].leadingTerm()/n.leadingTerm();
  }

  /////////////////////////////////////////////////////////////
  //////////////////generate double uwps///////////////////////
  /////////////////////////////////////////////////////////////
  vector<PolyMatrix> uwps_sz;
  for( int i = 0; i < (int) numberBearingVectors; i++ ) {
    uwps_sz.push_back(PolyMatrix(Poly::zeroSZ(nu),3,1));

    for( int j = 0; j < 3; j++ ) {      
      std::stringstream name; name << "uwps[" << i << "][" << j << "]";
      uwps_sz.back()[j] = Poly(Term(
        Coefficient(name.str()),
        uwps[i][j].leadingTerm().coefficient().clone(),
        Monomial(nu) ));
    }
  }

  /////////////////////////////////////////////////////////////
  //////////////////generate double vwps///////////////////////
  /////////////////////////////////////////////////////////////
  vector<PolyMatrix> vwps_sz;
  for( int i = 0; i < (int) numberBearingVectors; i++ ) {
    vwps_sz.push_back(PolyMatrix(Poly::zeroSZ(nu),3,1));

    for( int j = 0; j < 3; j++ ) {      
      std::stringstream name; name << "vwps[" << i << "][" << j << "]";
      vwps_sz.back()[j] = Poly(Term(
        Coefficient(name.str()),
        vwps[i][j].leadingTerm().coefficient().clone(),
        Monomial(nu) ));
    }
  }

  /////////////////////////////////////////////////////////////
  //////////////////generate double uwp_center/////////////////
  /////////////////////////////////////////////////////////////
  PolyMatrix uwp_center_sz(Poly::zeroSZ(nu),3,1);

  for( int j = 0; j < 3; j++ ) {
    std::stringstream name; name << "uwp_center[" << j << "]";
    uwp_center_sz[j] = Poly(Term(
        Coefficient(name.str()),
        uwp_center[j].leadingTerm().coefficient().clone(),
        Monomial(nu) ));
  }
  /////////////////////////////////////////////////////////////
  /////////////////generate double vwp_center//////////////////
  /////////////////////////////////////////////////////////////
  PolyMatrix vwp_center_sz(Poly::zeroSZ(nu),3,1);

  for( int j = 0; j < 3; j++ ) {
    std::stringstream name; name << "vwp_center[" << j << "]";
    vwp_center_sz[j] = Poly(Term(
        Coefficient(name.str()),
        vwp_center[j].leadingTerm().coefficient().clone(),
        Monomial(nu) ));
  }
  /////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////


  //3) construct unknowns
  Poly a = Poly::uSZ(1,nu);
  Poly b = Poly::uSZ(2,nu);
  Poly c = Poly::uSZ(3,nu);
  Poly d = Poly::uSZ(4,nu);

  PolyMatrix r1t(Poly::zeroSZ(nu),1,3);
  PolyMatrix r2t(Poly::zeroSZ(nu),1,3);
  PolyMatrix r3t(Poly::zeroSZ(nu),1,3);
  r1t(0,0) = a*a + b*b - c*c - d*d;
  r1t(0,1) = Poly::constSZ(2,nu)*b*c - Poly::constSZ(2,nu)*a*d;
  r1t(0,2) = Poly::constSZ(2,nu)*b*d + Poly::constSZ(2,nu)*a*c;
  r2t(0,0) = Poly::constSZ(2,nu)*b*c + Poly::constSZ(2,nu)*a*d;
  r2t(0,1) = a*a - b*b + c*c - d*d;
  r2t(0,2) = Poly::constSZ(2,nu)*c*d - Poly::constSZ(2,nu)*a*b;
  r3t(0,0) = Poly::constSZ(2,nu)*b*d - Poly::constSZ(2,nu)*a*c;
  r3t(0,1) = Poly::constSZ(2,nu)*c*d + Poly::constSZ(2,nu)*a*b;
  r3t(0,2) = a*a - b*b - c*c + d*d;

  //4) construct the main equation
  Poly t1 = f_center_sz[0] + (r3t * uwp_center_sz)(0,0) - (r1t * wp_center_sz)(0,0);
  Poly t2 = f_center_sz[1] + (r3t * vwp_center_sz)(0,0) - (r2t * wp_center_sz)(0,0);
  Poly energy = Poly::zeroSZ(nu);

  for( int i = 0; i < numberBearingVectors; i++ ) {
    Poly part1 = fs_sz[i][0] + (r3t*uwps_sz[i])(0,0) - (r1t * wps_sz[i])(0,0) - t1;
    Poly part2 = fs_sz[i][1] + (r3t*vwps_sz[i])(0,0) - (r2t * wps_sz[i])(0,0) - t2;
    energy += part1 * part1 + part2 * part2;
  }

  //5) take the derivatives
  list<Poly*> eqs; vector<Poly*> eqs_vec;
  for( int p = 0; p < nu; p++ ) {
    Poly * newPolynomial = new Poly(Poly::zeroSZ(nu));
    eqs.push_back(newPolynomial);
    eqs_vec.push_back(newPolynomial);
  }

  auto it = energy.begin();
  while( it != energy.end() ) {
    Term t = it->clone();

    for( int p = 0; p < nu; p++ ) {
      unsigned int constant = t.monomial().exponents()[p];
      if( constant > 0 ) {
        std::vector<unsigned int> newExponents = t.monomial().exponents();
        newExponents[p]--;
        Coefficient coeff_zp = t.coefficient().clone();
        t.setDominant(1);
        Coefficient coeff_sym = t.coefficient().clone();
        t.setDominant(0);
        (*eqs_vec[p]) += Poly::constSZ(constant,nu) * Poly(Term(coeff_sym,coeff_zp,Monomial(newExponents)));
      }
    }

    it++;
  }

  //******** Part 5: Generate the solver ******

  //find solver that respects symmetry (the following code is a bit ad-hoc as it is used to find the baseMonomials and the degree of expansion)

  /*
  //define the expanders (with a hypothetical expansion degree)
  std::vector<Monomial> expanders;
  expanders.push_back(a.leadingTerm().monomial());
  expanders.push_back(b.leadingTerm().monomial());
  expanders.push_back(c.leadingTerm().monomial());
  expanders.push_back(d.leadingTerm().monomial());

  int expanderDegree = 6;
  std::vector<Monomial> expanders2;
  methods::generateEvendegreeExpanders(expanders,expanders2,expanderDegree);
  expanders = expanders2;

  //attempt solution and print all the monomials
  CMatrix attempt = methods::experiment(eqs,expanders,true);
  std::vector<core::Monomial> monomials2 = attempt.monomials();
  for( int i = 0; i < monomials2.size(); i++ )
  {
    monomials2[i].print(); std::cout << std::endl;
  }
  return 0;
  */
  

  //Once the above code has been used to find the degree of expansion and the base monomials, proceed

  //expansion degree
  int expanderDegree = 6;

  //monomials (found by visual inspection)
  std::vector<Monomial> baseMonomials;

  unsigned int exp01[4] = {1,0,0,6}; baseMonomials.push_back(Monomial(nu,exp01)); // x_1*x_4^6
  unsigned int exp02[4] = {0,1,0,6}; baseMonomials.push_back(Monomial(nu,exp02)); // x_2*x_4^6
  unsigned int exp03[4] = {0,0,1,6}; baseMonomials.push_back(Monomial(nu,exp03)); // x_3*x_4^6
  unsigned int exp04[4] = {0,0,0,7}; baseMonomials.push_back(Monomial(nu,exp04)); // x_4^7
  unsigned int exp05[4] = {2,0,1,2}; baseMonomials.push_back(Monomial(nu,exp05)); // x_1^2*x_3*x_4^2
  unsigned int exp06[4] = {1,1,1,2}; baseMonomials.push_back(Monomial(nu,exp06)); // x_1*x_2*x_3*x_4^2
  unsigned int exp07[4] = {0,2,1,2}; baseMonomials.push_back(Monomial(nu,exp07)); // x_2^2*x_3*x_4^2
  unsigned int exp08[4] = {1,0,2,2}; baseMonomials.push_back(Monomial(nu,exp08)); // x_1*x_3^2*x_4^2
  unsigned int exp09[4] = {0,1,2,2}; baseMonomials.push_back(Monomial(nu,exp09)); // x_2*x_3^2*x_4^2
  unsigned int exp10[4] = {0,0,3,2}; baseMonomials.push_back(Monomial(nu,exp10)); // x_3^3*x_4^2
  unsigned int exp11[4] = {2,0,0,3}; baseMonomials.push_back(Monomial(nu,exp11)); // x_1^2*x_4^3
  unsigned int exp12[4] = {1,1,0,3}; baseMonomials.push_back(Monomial(nu,exp12)); // x_1*x_2*x_4^3
  unsigned int exp13[4] = {0,2,0,3}; baseMonomials.push_back(Monomial(nu,exp13)); // x_2^2*x_4^3
  unsigned int exp14[4] = {1,0,1,3}; baseMonomials.push_back(Monomial(nu,exp14)); // x_1*x_3*x_4^3
  unsigned int exp15[4] = {0,1,1,3}; baseMonomials.push_back(Monomial(nu,exp15)); // x_2*x_3*x_4^3
  unsigned int exp16[4] = {0,0,2,3}; baseMonomials.push_back(Monomial(nu,exp16)); // x_3^2*x_4^3
  unsigned int exp17[4] = {1,0,0,4}; baseMonomials.push_back(Monomial(nu,exp17)); // x_1*x_4^4
  unsigned int exp18[4] = {0,1,0,4}; baseMonomials.push_back(Monomial(nu,exp18)); // x_2*x_4^4
  unsigned int exp19[4] = {0,0,1,4}; baseMonomials.push_back(Monomial(nu,exp19)); // x_3*x_4^4
  unsigned int exp20[4] = {0,0,0,5}; baseMonomials.push_back(Monomial(nu,exp20)); // x_4^5
  unsigned int exp21[4] = {2,0,1,0}; baseMonomials.push_back(Monomial(nu,exp21)); // x_1^2*x_3
  unsigned int exp22[4] = {1,1,1,0}; baseMonomials.push_back(Monomial(nu,exp22)); // x_1*x_2*x_3
  unsigned int exp23[4] = {0,2,1,0}; baseMonomials.push_back(Monomial(nu,exp23)); // x_2^2*x_3
  unsigned int exp24[4] = {1,0,2,0}; baseMonomials.push_back(Monomial(nu,exp24)); // x_1*x_3^2
  unsigned int exp25[4] = {0,1,2,0}; baseMonomials.push_back(Monomial(nu,exp25)); // x_2*x_3^2
  unsigned int exp26[4] = {0,0,3,0}; baseMonomials.push_back(Monomial(nu,exp26)); // x_3^3
  unsigned int exp27[4] = {2,0,0,1}; baseMonomials.push_back(Monomial(nu,exp27)); // x_1^2*x_4
  unsigned int exp28[4] = {1,1,0,1}; baseMonomials.push_back(Monomial(nu,exp28)); // x_1*x_2*x_4
  unsigned int exp29[4] = {0,2,0,1}; baseMonomials.push_back(Monomial(nu,exp29)); // x_2^2*x_4
  unsigned int exp30[4] = {1,0,1,1}; baseMonomials.push_back(Monomial(nu,exp30)); // x_1*x_3*x_4
  unsigned int exp31[4] = {0,1,1,1}; baseMonomials.push_back(Monomial(nu,exp31)); // x_2*x_3*x_4
  unsigned int exp32[4] = {0,0,2,1}; baseMonomials.push_back(Monomial(nu,exp32)); // x_3^2*x_4
  unsigned int exp33[4] = {1,0,0,2}; baseMonomials.push_back(Monomial(nu,exp33)); // x_1*x_4^2
  unsigned int exp34[4] = {0,1,0,2}; baseMonomials.push_back(Monomial(nu,exp34)); // x_2*x_4^2
  unsigned int exp35[4] = {0,0,1,2}; baseMonomials.push_back(Monomial(nu,exp35)); // x_3*x_4^2
  unsigned int exp36[4] = {0,0,0,3}; baseMonomials.push_back(Monomial(nu,exp36)); // x_4^3
  unsigned int exp37[4] = {1,0,0,0}; baseMonomials.push_back(Monomial(nu,exp37)); // x_1
  unsigned int exp38[4] = {0,1,0,0}; baseMonomials.push_back(Monomial(nu,exp38)); // x_2
  unsigned int exp39[4] = {0,0,1,0}; baseMonomials.push_back(Monomial(nu,exp39)); // x_3
  unsigned int exp40[4] = {0,0,0,1}; baseMonomials.push_back(Monomial(nu,exp40)); // x_4

  //generate the symmetric solver
  std::string solverName("opnp");
  string parameters("std::vector<Eigen::Vector3d> & fs, std::vector<Eigen::Vector3d> & wps, Eigen::Vector3d & f_center, Eigen::Vector3d & wp_center, std::vector<Eigen::Vector3d> & uwps, std::vector<Eigen::Vector3d> & vwps, Eigen::Vector3d & uwp_center, Eigen::Vector3d & vwp_center");
  execGeneratorSym( eqs, expanderDegree, baseMonomials, solverName, parameters, true );
}
  