#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <sys/time.h>

#include "gen3pt.hpp"


using namespace std;
using namespace Eigen;


Matrix3d
cayley2rot( const Vector3d & cayley)
{
  Matrix3d R;
  double scale = 1+pow(cayley[0],2)+pow(cayley[1],2)+pow(cayley[2],2);

  R(0,0) = 1+pow(cayley[0],2)-pow(cayley[1],2)-pow(cayley[2],2);
  R(0,1) = 2*(cayley[0]*cayley[1]-cayley[2]);
  R(0,2) = 2*(cayley[0]*cayley[2]+cayley[1]);
  R(1,0) = 2*(cayley[0]*cayley[1]+cayley[2]);
  R(1,1) = 1-pow(cayley[0],2)+pow(cayley[1],2)-pow(cayley[2],2);
  R(1,2) = 2*(cayley[1]*cayley[2]-cayley[0]);
  R(2,0) = 2*(cayley[0]*cayley[2]-cayley[1]);
  R(2,1) = 2*(cayley[1]*cayley[2]+cayley[0]);
  R(2,2) = 1-pow(cayley[0],2)-pow(cayley[1],2)+pow(cayley[2],2);

  R = (1/scale) * R;
  return R;
}

int main( int argc, char** argv )
{
  //initialize random seed
  struct timeval tic;
  gettimeofday( &tic, 0 );
  srand ( tic.tv_usec );
  
  //set experiment parameters
  size_t numberPoints = 10;
  int numberCameras = 4;

  //generate random translation
  Vector3d translation;
  for( int i = 0; i < 3; i++ )
    translation[i] = (((double) rand())/ ((double) RAND_MAX)-0.5)*4.0;

  //generate random rotation
  Vector3d cay;
  for( int i = 0; i < 3; i++ )
    cay[i] = ((double) rand())/ ((double) RAND_MAX)-0.5;
  Matrix3d R = cayley2rot(cay);
  
  //create a random camera-system
  vector<Vector3d> camOffsets;
  vector<Matrix3d> camRotations;
  
  for( int i = 0; i < numberCameras; i++ ) {
    Vector3d c1, c2;
    for( int i = 0; i < 3; i++ ) {
      c1[i] = (((double) rand())/ ((double) RAND_MAX)-0.5)*2.0;
      c2[i] = (((double) rand())/ ((double) RAND_MAX)-0.5)*2.0;
    }
    camOffsets.push_back(cayley2rot(c1).col(0)*0.5);
    camRotations.push_back(cayley2rot(c2));
  }
  
  //generate random 3D points
  MatrixXd gt(3,numberPoints);
  for( size_t i = 0; i < (size_t) gt.cols(); i++ ) {
    Vector3d cleanPoint;
    for( int j = 0; j < 3; j++ )
      cleanPoint[j] = (((double) rand())/ ((double) RAND_MAX)-0.5)*2.0;
    Vector3d direction = cleanPoint / cleanPoint.norm();
    gt.col(i) = (8.0-4.0) * cleanPoint + 4.0 * direction;
  }

  //derive image points (and camera correspondence vector)
  vector<Vector3d> bearingVectors;
  vector<Vector3d> points;
  vector<int> camCorrespondences;

  size_t camCorrespondence = 0;
  
  for( size_t i = 0; i < (size_t) gt.cols(); i++ ) {
    
    //get the camera transformation
    Eigen::Vector3d camOffset = camOffsets[camCorrespondence];
    Eigen::Matrix3d camRotation = camRotations[camCorrespondence];

    //store the point
    points.push_back(gt.col(i));
    
    //project the point into the viewpoint frame
    Eigen::Vector3d bodyPoint = R.transpose()*(gt.col(i) - translation);
    
    //project the point into the camera frame
    bearingVectors.push_back(camRotation.transpose()*(bodyPoint - camOffset));

    //normalize the bearing-vector to 1
    bearingVectors[i] = bearingVectors[i] / bearingVectors[i].norm();
    
    //push back the camera correspondence
    camCorrespondences.push_back(camCorrespondence++);
    if(camCorrespondence > (numberCameras-1) )
      camCorrespondence = 0;
  }

  //print ground truth
  cout << translation.transpose() << endl << endl;
  cout << R << endl << endl;

  //Change the format for our derivation (and limit to minimum case of 3 points)
  vector<Vector3d> realBearingVectors;
  vector<Vector3d> realCamOffsets;
  vector<Vector3d> realPoints;
  for( int i = 0; i < 3; i++ ) {
    realBearingVectors.push_back(camRotations[camCorrespondences[i]] * bearingVectors[i]);
    realCamOffsets.push_back(camOffsets[camCorrespondences[i]]);
    realPoints.push_back(points[i]);
  }

  //run the experiments
  cout << "running the solver\n";
  vector< Matrix<double,6,1> > solutions;
  polyjam::gen3pt::solve( realBearingVectors, realCamOffsets, realPoints, solutions );

  //transform from Action matrix to real solutions
  vector< Matrix<double,3,4> > gp3p_transformations;
  
  for( int c = 0; c < solutions.size(); c++ ) {
    Vector3d cayley;
    Vector3d n;

    for(size_t i = 0; i < 3; i++) {
      complex<double> cay = solutions[c](i+3,0);
      cayley[i] = cay.real();
      complex<double> depth = solutions[c](i,0);
      n[i] = depth.real();
    }

    Matrix3d rotation = cayley2rot(cayley);

    Vector3d center_cam = Vector3d::Zero();
    Vector3d center_world = Vector3d::Zero();
    for( size_t i = 0; i < (size_t) realBearingVectors.size(); i++ ) {
      Vector3d temp = rotation*(n[i]*realBearingVectors[i]+realCamOffsets[i]);
      center_cam = center_cam + temp;
      center_world = center_world + realPoints[i];
    }

    center_cam = center_cam/realBearingVectors.size();
    center_world = center_world/realBearingVectors.size();
    Vector3d translation = center_world - center_cam;

    Matrix<double,3,4> transformation;
    transformation.block<3,3>(0,0) = rotation;
    transformation.col(3) = translation;
    gp3p_transformations.push_back(transformation);
  }

  //print the results
  cout << "results from gp3p algorithm:" << endl;
  for(size_t i = 0; i < gp3p_transformations.size(); i++)
    cout << gp3p_transformations[i] << endl << endl;
}
