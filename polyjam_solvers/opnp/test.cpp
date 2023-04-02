#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <sys/time.h>

#include "opnp.hpp"


using namespace std;
using namespace Eigen;


int main( int argc, char** argv )
{
  //initialize random seed
  struct timeval tic;
  gettimeofday( &tic, 0 );
  srand ( tic.tv_usec );

  //settings
  size_t numberBearingVectors = 3;

  //first create a normalized quaternion
  Eigen::Matrix<double,4,1> q_gt;
  for( int i = 0; i < 4; i++ )
    q_gt[i] = (((double) rand())/ ((double) RAND_MAX)-0.5)*2.0;
  q_gt /= q_gt.norm();

  //Then sample random depths and get average depth
  std::vector<double> depths;
  double averageDepth = 0.0;
  for( int i = 0; i < numberBearingVectors; i++ ) {
    depths.push_back( ((double) rand())/((double) RAND_MAX)*1.0 );
    averageDepth += depths.back();
  }
  averageDepth /= (double) numberBearingVectors;

  //now resample the quaternion to make sure it is what we want
  q_gt *= sqrt(averageDepth);

  //derive the scaled rotation matrix
  Eigen::Matrix3d R_gt;
  R_gt(0,0) = (q_gt[0]*q_gt[0]) + (q_gt[1]*q_gt[1]) - (q_gt[2]*q_gt[2]) - (q_gt[3]*q_gt[3]);
  R_gt(1,1) = (q_gt[0]*q_gt[0]) - (q_gt[1]*q_gt[1]) + (q_gt[2]*q_gt[2]) - (q_gt[3]*q_gt[3]);
  R_gt(2,2) = (q_gt[0]*q_gt[0]) - (q_gt[1]*q_gt[1]) - (q_gt[2]*q_gt[2]) + (q_gt[3]*q_gt[3]);
  R_gt(0,1) = 2.0 * (q_gt[1]*q_gt[2]-q_gt[3]*q_gt[0]);
  R_gt(0,2) = 2.0 * (q_gt[1]*q_gt[3]+q_gt[2]*q_gt[0]);
  R_gt(1,2) = 2.0 * (q_gt[2]*q_gt[3]-q_gt[1]*q_gt[0]);
  R_gt(1,0) = 2.0 * (q_gt[1]*q_gt[2]+q_gt[3]*q_gt[0]);
  R_gt(2,0) = 2.0 * (q_gt[1]*q_gt[3]-q_gt[2]*q_gt[0]);
  R_gt(2,1) = 2.0 * (q_gt[2]*q_gt[3]+q_gt[1]*q_gt[0]);

  //random translation (this is not the true translation, but a scale one again)
  Eigen::Vector3d t_gt;
  for( int i = 0; i < 3; i++ )
    t_gt[i] = (((double) rand())/ ((double) RAND_MAX)-0.5)*2.0 / averageDepth;

  //random image measurements
  vector<Eigen::Vector3d> fs(numberBearingVectors);
  for( int i = 0; i < numberBearingVectors; i++ ) {
    fs[i][0] = (((double) rand())/ ((double) RAND_MAX)-0.5)*2.0;
    fs[i][1] = (((double) rand())/ ((double) RAND_MAX)-0.5)*2.0;
    fs[i][2] = 1.0;
  }

  //derive the world points
  vector<Eigen::Vector3d> wps(numberBearingVectors);
  for( int i = 0; i < numberBearingVectors; i++ ) {
    double scaledDepth = depths[i] / averageDepth;
    wps[i] = R_gt.inverse() * ( fs[i] * scaledDepth - t_gt);
  }

  //now derive all the other matrices that are needed

  //center of image points
  Vector3d f_center = Vector3d::Zero();
  for( int i = 0; i < numberBearingVectors; i++ )
    f_center += fs[i];
  f_center /= (double) numberBearingVectors;

  //center of the world points
  Vector3d wp_center = Vector3d::Zero();
  for( int i = 0; i < numberBearingVectors; i++ )
    wp_center += wps[i];
  wp_center /= (double) numberBearingVectors;

  //create inner summation terms
  vector<Vector3d> uwps, vwps;
  Vector3d uwp_center = Vector3d::Zero();
  Vector3d vwp_center = Vector3d::Zero();
  for( int i = 0; i < numberBearingVectors; i++ ) {
    uwps.push_back( (wps[i]-wp_center) * fs[i][0] ); uwp_center += uwps.back();
    vwps.push_back( (wps[i]-wp_center) * fs[i][1] ); vwp_center += vwps.back();
  }
  uwp_center /= (double) numberBearingVectors;
  vwp_center /= (double) numberBearingVectors;

  //call the solver
  std::vector< Eigen::Matrix<double,4,1> > solutions;
  polyjam::opnp::solve( fs, wps, f_center, wp_center, uwps, vwps, uwp_center, vwp_center, solutions );

  //print the results
  cout << "groundtruth result: " << q_gt.transpose() / q_gt.norm() << endl << endl;

  for( int i = 0; i < solutions.size(); i++ )
    cout << "solution " << i << ": " << solutions[i].transpose() / solutions[i].norm() << endl;

  return 0;
}
