#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include "sw6pt.hpp"
#include <sstream>
#include <fstream>

#include "random_generators.hpp"


using namespace std;
using namespace Eigen;
using namespace opengv;


int main( int argc, char** argv )
{
  //initialize random seed
  initializeRandomSeed();
  
  //Initialize three random matrices F1, F2, and F3
  Eigen::Matrix3d F1, F2, F3;

  for( int r = 0; r < 3; r++ ) {
    for( int c = 0; c < 3; c++ ) {
      
      F1(r,c) = (((double) rand())/ ((double) RAND_MAX)-0.5)*2.0;
      F2(r,c) = (((double) rand())/ ((double) RAND_MAX)-0.5)*2.0;
      F3(r,c) = (((double) rand())/ ((double) RAND_MAX)-0.5)*2.0;

    }
  }

  //call the solver
  std::vector< Eigen::Matrix<double,3,1> > solutions;
  polyjam::sw6pt::solve( F1, F2, F3, solutions );

  //verify the solution
  for( int i = 0; i < solutions.size(); i++ ) {

    std::cout << "Verifying the solution " << i <<  ":\n";

    double x,y,w;
    x = solutions[i](0,0);
    y = solutions[i](1,0);
    w = solutions[i](2,0);

    Eigen::Matrix3d F = F1 * x + F2 * y + F3;
    Eigen::Matrix3d Ft = F.transpose();
    Eigen::Matrix3d Q = Eigen::Matrix3d::Identity();
    Q(2,2) = w;
    Eigen::Matrix3d FQFtQ = F*Q*Ft*Q;
    Eigen::Matrix3d te = FQFtQ * F * 2.0 - F * (FQFtQ(0,0) + FQFtQ(1,1) + FQFtQ(2,2));

    std::cout << "Determinant constraint: " << F.determinant() << std::endl;

    for( int r = 0; r < 3; r++ ) {
      for( int c = 0; c < 3; c++ )
        std::cout << "Next value: " << te(r,c) << std::endl;
    }
    std::cout << std::endl;

  }

  return 0;
}
