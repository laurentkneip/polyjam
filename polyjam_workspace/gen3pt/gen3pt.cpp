#include <polyjam/polyjam.hpp>

int main( int argc, char** argv )
{  
  //initialize the random generator
  initGenerator();
  size_t nu = 6; //the number of unknowns in the problem
  size_t numberBearingVectors = 3;

  //****** Part 1: get independent random variables *******//

  //random rotation
  std::cout << "generating random rotation" << std::endl;
  PolyMatrix cay_gt(Poly::zeroZ(nu),3,1);
  for( int i = 0; i < 3; i++ )
    cay_gt[i] = Poly::randZ(nu);
  
  Poly scale12 = Poly::oneZ(nu) + (cay_gt[0]*cay_gt[0]) + (cay_gt[1]*cay_gt[1]) + (cay_gt[2]*cay_gt[2]);
  Poly scale12_inv = Poly::oneZ(nu).leadingTerm() / scale12.leadingTerm();
  
  PolyMatrix R_gt(Poly::zeroZ(nu),3,3);
  R_gt(0,0) = scale12_inv * ( Poly::oneZ(nu) + (cay_gt[0]*cay_gt[0]) - (cay_gt[1]*cay_gt[1]) - (cay_gt[2]*cay_gt[2]) );
  R_gt(1,1) = scale12_inv * ( Poly::oneZ(nu) - (cay_gt[0]*cay_gt[0]) + (cay_gt[1]*cay_gt[1]) - (cay_gt[2]*cay_gt[2]) );
  R_gt(2,2) = scale12_inv * ( Poly::oneZ(nu) - (cay_gt[0]*cay_gt[0]) - (cay_gt[1]*cay_gt[1]) + (cay_gt[2]*cay_gt[2]) );
  R_gt(0,1) = scale12_inv * ( Poly::constZ(2,nu) * (cay_gt[0]*cay_gt[1]-cay_gt[2]) );
  R_gt(0,2) = scale12_inv * ( Poly::constZ(2,nu) * (cay_gt[0]*cay_gt[2]+cay_gt[1]) );
  R_gt(1,2) = scale12_inv * ( Poly::constZ(2,nu) * (cay_gt[1]*cay_gt[2]-cay_gt[0]) );
  R_gt(1,0) = scale12_inv * ( Poly::constZ(2,nu) * (cay_gt[0]*cay_gt[1]+cay_gt[2]) );
  R_gt(2,0) = scale12_inv * ( Poly::constZ(2,nu) * (cay_gt[0]*cay_gt[2]-cay_gt[1]) );
  R_gt(2,1) = scale12_inv * ( Poly::constZ(2,nu) * (cay_gt[1]*cay_gt[2]+cay_gt[0]) );
  
  //random translation
  std::cout << "generating random translation" << std::endl;
  PolyMatrix t_gt(Poly::zeroZ(nu),3,1);
  for( int i = 0; i < 3; i++ )
    t_gt[i] = Poly::randZ(nu);
  
  //random world points and offsets
  std::cout << "generating random world points and offsets" << std::endl;
  std::vector<PolyMatrix> wps;
  std::vector<PolyMatrix> vs;
  for( int i = 0; i < (int) numberBearingVectors; i++ )
  {
    wps.push_back(PolyMatrix(Poly::zeroZ(nu),3,1));
    vs.push_back(PolyMatrix(Poly::zeroZ(nu),3,1));

    for( int j = 0; j < 3; j++ )
    {
      wps.back()[j] = Poly::randZ(nu);
      vs.back()[j] = Poly::randZ(nu);
    }
  }

  //****** Part 2: extract the measurements with both Zp and symbolic representation *****//

  std::cout << "extracting the image measurements" << std::endl;
  std::vector<PolyMatrix> fs;
  
  for( int i = 0; i < (int) numberBearingVectors; i++ )
  {
    PolyMatrix f = R_gt.transpose() * ( wps[i] - t_gt ) - vs[i];

    //do some normalization here (just divide by the last coordinate!)
    Poly f_scale = f[2].clone();
    for( int j = 0; j < 3; j++ )
      f[j] = f[j].leadingTerm()/f_scale.leadingTerm();
    
    fs.push_back(PolyMatrix(Poly::zeroSZ(nu),3,1));
    
    for( int j = 0; j < 3; j++ )
    {      
      std::stringstream name; name << "fs[" << i << "][" << j << "]";
      fs.back()[j] = Poly(Term(
        Coefficient(name.str()),
        f[j].leadingTerm().coefficient().clone(),
        Monomial(nu) ));
    }
  }

  //****** Part 3: also lift the 3D points to the double representation ******//

  std::cout << "lifting the 3D points" << std::endl;
  std::vector<PolyMatrix> wps_double;
  std::vector<PolyMatrix> vs_double;
  for( int i = 0; i < (int) numberBearingVectors; i++ )
  {
    wps_double.push_back(PolyMatrix(Poly::zeroSZ(nu),3,1));
    vs_double.push_back(PolyMatrix(Poly::zeroSZ(nu),3,1));

    for( int j = 0; j < 3; j++ )
    {
      std::stringstream name1; name1 << "wps[" << i << "][" << j << "]";
      wps_double.back()[j] = Poly(Term(
        Coefficient(name1.str()),
        wps[i][j].leadingTerm().coefficient().clone(),
        Monomial(nu) ));
      
      std::stringstream name2; name2 << "vs[" << i << "][" << j << "]";
      vs_double.back()[j] = Poly(Term(
        Coefficient(name2.str()),
        vs[i][j].leadingTerm().coefficient().clone(),
        Monomial(nu) ));
    }
  }

  //****** Part 4: get all the equations ****************//

  std::cout << "Extracting the equations" << std::endl;

  std::cout << "defining unknowns" << std::endl;
  Poly x = Poly::uSZ(4,nu);
  Poly y = Poly::uSZ(5,nu);
  Poly z = Poly::uSZ(6,nu);
  
  Poly n1 = Poly::uSZ(1,nu);
  Poly n2 = Poly::uSZ(2,nu);
  Poly n3 = Poly::uSZ(3,nu);

  std::cout << "defining rotation" << std::endl;
  PolyMatrix R(Poly::zeroSZ(nu),3,3);
  R(0,0) = Poly::oneSZ(nu) + (x*x) - (y*y) - (z*z);
  R(1,1) = Poly::oneSZ(nu) - (x*x) + (y*y) - (z*z);
  R(2,2) = Poly::oneSZ(nu) - (x*x) - (y*y) + (z*z);
  R(0,1) = Poly::constSZ(2,nu) * (x*y-z);
  R(0,2) = Poly::constSZ(2,nu) * (x*z+y);
  R(1,2) = Poly::constSZ(2,nu) * (y*z-x);
  R(1,0) = Poly::constSZ(2,nu) * (x*y+z);
  R(2,0) = Poly::constSZ(2,nu) * (x*z-y);
  R(2,1) = Poly::constSZ(2,nu) * (y*z+x);
  
  std::cout << "defining scale" << std::endl;
  Poly scale = Poly::oneSZ(nu) + (x*x) + (y*y) + (z*z);

  std::cout << "defining the M matrices" << std::endl;
  PolyMatrix M1 = ( fs[0]*n1 - fs[1]*n2 + vs_double[0] - vs_double[1] ) * scale -
      R.transpose()*wps_double[0] + R.transpose()*wps_double[1];
  PolyMatrix M2 = ( fs[1]*n2 - fs[2]*n3 + vs_double[1] - vs_double[2] ) * scale -
      R.transpose()*wps_double[1] + R.transpose()*wps_double[2];
  PolyMatrix M3 = ( fs[2]*n3 - fs[0]*n1 + vs_double[2] - vs_double[0] ) * scale -
      R.transpose()*wps_double[2] + R.transpose()*wps_double[0];
  
  std::cout << "extracting the equations" << std::endl;
  std::list<Poly*> eqs;
  for( size_t i = 0; i < 3; i++ )
    eqs.push_back(new Poly(M1[i]));
  for( size_t i = 0; i < 3; i++ )
    eqs.push_back(new Poly(M2[i]));
  for( size_t i = 0; i < 3; i++ )
    eqs.push_back(new Poly(M3[i]));

  //****** Part 5: generate the solver ********************//

  execGenerator( eqs, string("gen3pt"), string("std::vector<Eigen::Vector3d> & fs, std::vector<Eigen::Vector3d> & vs, std::vector<Eigen::Vector3d> & wps"), true, "." );
}
