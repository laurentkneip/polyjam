#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include "gen3pt.hpp"
#include <sstream>
#include <fstream>

#include "random_generators.hpp"
#include "experiment_helpers.hpp"
#include "time_measurement.hpp"


using namespace std;
using namespace Eigen;
using namespace opengv;


Eigen::Matrix3d
cayley2rot( const Eigen::Vector3d & cayley)
{
  Eigen::Matrix3d R;
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
  initializeRandomSeed();
  
  //set experiment parameters
  double noise = 0.0;
  double outlierFraction = 0.0;
  size_t numberPoints = 100;
  int numberCameras = 4;

  //create a random viewpoint pose
  Eigen::Vector3d position = generateRandomTranslation(2.0);
  Eigen::Matrix3d rotation = generateRandomRotation(0.5);
  
  //create a random camera-system
  std::vector<Eigen::Vector3d> camOffsets;
  std::vector<Eigen::Matrix3d> camRotations;
  generateRandomCameraSystem( numberCameras, camOffsets, camRotations );
  
  //derive correspondences based on random point-cloud
  std::vector<Eigen::Vector3d> bearingVectors;
  std::vector<Eigen::Vector3d> points;
  std::vector<int> camCorrespondences;
  Eigen::MatrixXd gt(3,numberPoints);
  generateRandom2D3DCorrespondences(
      position, rotation, camOffsets, camRotations, numberPoints, noise, outlierFraction,
      bearingVectors, points, camCorrespondences, gt );

  std::vector<Eigen::Vector3d> realBearingVectors;
  std::vector<Eigen::Vector3d> realCamOffsets;
  std::vector<Eigen::Vector3d> realPoints;
  for( int i = 0; i < 3; i++ )
  {
    realBearingVectors.push_back(camRotations[camCorrespondences[i]] * bearingVectors[i]);
    realCamOffsets.push_back(camOffsets[camCorrespondences[i]]);
    realPoints.push_back(points[i]);
  }

  //print the experiment characteristics
  printExperimentCharacteristics(
      position, rotation, noise, outlierFraction );

  //timer
  struct timeval tic;
  struct timeval toc;
  size_t iterations = 50;

  //run the experiments
  std::cout << "running Kneip's GP3P (using first three correspondences/";
  std::cout << std::endl;
  std::vector< Eigen::Matrix<double,6,1> > solutions;
  gettimeofday( &tic, 0 );
  for(size_t i = 0; i < iterations; i++)
  {
    solutions.clear();
    polyjam::gen3pt::solve( realBearingVectors, realCamOffsets, realPoints, solutions );
  }
  gettimeofday( &toc, 0 );
  double gp3p_time = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;

  //todo: transform from Action matrix to real solutions
  std::vector< Eigen::Matrix<double,3,4> > gp3p_transformations;
  
  for( int c = 0; c < solutions.size(); c++ )
  {
      Eigen::Vector3d cayley;
      Eigen::Vector3d n;

      for(size_t i = 0; i < 3; i++)
      {
        std::complex<double> cay = solutions[c](i+3,0);
        cayley[i] = cay.real();
        std::complex<double> depth = solutions[c](i,0);
        n[i] = depth.real();
      }

      Eigen::Matrix3d rotation = cayley2rot(cayley);

      Eigen::Vector3d center_cam = Eigen::Vector3d::Zero();
      Eigen::Vector3d center_world = Eigen::Vector3d::Zero();
      for( size_t i = 0; i < (size_t) realBearingVectors.size(); i++ )
      {
        Eigen::Vector3d temp = rotation*(n[i]*realBearingVectors[i]+realCamOffsets[i]);
        center_cam = center_cam + temp;
        center_world = center_world + realPoints[i];
      }

      center_cam = center_cam/realBearingVectors.size();
      center_world = center_world/realBearingVectors.size();
      Eigen::Vector3d translation = center_world - center_cam;

      Eigen::Matrix<double,3,4> transformation;
      transformation.block<3,3>(0,0) = rotation;
      transformation.col(3) = translation;
      gp3p_transformations.push_back(transformation);
  }

  //print the results
  std::cout << "results from gp3p algorithm:" << std::endl;
  for(size_t i = 0; i < gp3p_transformations.size(); i++)
    std::cout << gp3p_transformations[i] << std::endl << std::endl;

  std::cout << "timings from gp3p algorithm: ";
  std::cout << gp3p_time << std::endl;
}
