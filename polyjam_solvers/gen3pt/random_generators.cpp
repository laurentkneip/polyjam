#include "random_generators.hpp"
#include "time_measurement.hpp"
#include <math.h>
#include <iostream>


using namespace Eigen;

void
opengv::initializeRandomSeed()
{
  struct timeval tic;
  gettimeofday( &tic, 0 );
  srand ( tic.tv_usec );
}

Eigen::Vector3d
opengv::generateRandomPoint( double maximumDepth, double minimumDepth )
{
  Eigen::Vector3d cleanPoint;
  cleanPoint[0] = (((double) rand())/ ((double) RAND_MAX)-0.5)*2.0;
  cleanPoint[1] = (((double) rand())/ ((double) RAND_MAX)-0.5)*2.0;
  cleanPoint[2] = (((double) rand())/ ((double) RAND_MAX)-0.5)*2.0;
  Eigen::Vector3d direction = cleanPoint / cleanPoint.norm();
  cleanPoint =
      (maximumDepth-minimumDepth) * cleanPoint + minimumDepth * direction;
  return cleanPoint;
}

Eigen::Vector3d
opengv::generateRandomPointPlane()
{
  Eigen::Vector3d cleanPoint;
  cleanPoint[0] = (((double) rand())/ ((double) RAND_MAX)-0.5)*2.0;
  cleanPoint[1] = (((double) rand())/ ((double) RAND_MAX)-0.5)*2.0;
  cleanPoint[2] = (((double) rand())/ ((double) RAND_MAX)-0.5)*2.0;

  cleanPoint[0] = 6*cleanPoint[0];
  cleanPoint[1] = 6*cleanPoint[1];
  cleanPoint[2] = 2*cleanPoint[2]-6.0;

  return cleanPoint;
}

Eigen::Vector3d
opengv::addNoise( double noiseLevel, Eigen::Vector3d cleanPoint )
{
  //compute a vector in the normal plane (based on good conditioning)
  Eigen::Vector3d normalVector1;
  if(
      (fabs(cleanPoint[0]) > fabs(cleanPoint[1])) &&
      (fabs(cleanPoint[0]) > fabs(cleanPoint[2])) )
  {
    normalVector1[1] = 1.0;
    normalVector1[2] = 0.0;
    normalVector1[0] = -cleanPoint[1]/cleanPoint[0];
  }
  else
  {
    if(
        (fabs(cleanPoint[1]) > fabs(cleanPoint[0])) &&
        (fabs(cleanPoint[1]) > fabs(cleanPoint[2])) )
    {
      normalVector1[2] = 1.0;
      normalVector1[0] = 0.0;
      normalVector1[1] = -cleanPoint[2]/cleanPoint[1];
    }
    else
    {
      normalVector1[0] = 1.0;
      normalVector1[1] = 0.0;
      normalVector1[2] = -cleanPoint[0]/cleanPoint[2];
    }
  }

  normalVector1 = normalVector1 / normalVector1.norm();
  Eigen::Vector3d normalVector2 = cleanPoint.cross(normalVector1);
  double noiseX =
      noiseLevel * (((double) rand())/ ((double) RAND_MAX)-0.5)*2.0 / 1.4142;
  double noiseY =
      noiseLevel * (((double) rand())/ ((double) RAND_MAX)-0.5)*2.0 / 1.4142;

  Eigen::Vector3d noisyPoint =
      800 * cleanPoint + noiseX *normalVector1 + noiseY * normalVector2;
  noisyPoint = noisyPoint / noisyPoint.norm();
  return noisyPoint;

}

Eigen::Vector3d
opengv::generateRandomTranslation( double maximumParallax )
{
  Eigen::Vector3d translation;
  translation[0] = (((double) rand())/ ((double) RAND_MAX)-0.5)*2.0;
  translation[1] = (((double) rand())/ ((double) RAND_MAX)-0.5)*2.0;
  translation[2] = (((double) rand())/ ((double) RAND_MAX)-0.5)*2.0;
  return maximumParallax * translation;
}

Eigen::Vector3d
opengv::generateRandomDirectionTranslation( double parallax )
{
  Matrix3d rotation = generateRandomRotation();
  Vector3d translation;
  translation << 1.0, 0.0, 0.0;
  translation = rotation * translation;
  translation = parallax * translation;
  return translation;
}

Eigen::Matrix3d
opengv::generateRandomRotation( double maxAngle )
{
  Eigen::Vector3d rpy;
  rpy[0] = ((double) rand())/ ((double) RAND_MAX);
  rpy[1] = ((double) rand())/ ((double) RAND_MAX);
  rpy[2] = ((double) rand())/ ((double) RAND_MAX);

  rpy[0] = maxAngle*2.0*(rpy[0]-0.5);
  rpy[1] = maxAngle*2.0*(rpy[1]-0.5);
  rpy[2] = maxAngle*2.0*(rpy[2]-0.5);

  Matrix3d R1;
  R1(0,0) = 1.0;
  R1(0,1) = 0.0;
  R1(0,2) = 0.0;
  R1(1,0) = 0.0;
  R1(1,1) = cos(rpy[0]);
  R1(1,2) = -sin(rpy[0]);
  R1(2,0) = 0.0;
  R1(2,1) = -R1(1,2);
  R1(2,2) = R1(1,1);

  Matrix3d R2;
  R2(0,0) = cos(rpy[1]);
  R2(0,1) = 0.0;
  R2(0,2) = sin(rpy[1]);
  R2(1,0) = 0.0;
  R2(1,1) = 1.0;
  R2(1,2) = 0.0;
  R2(2,0) = -R2(0,2);
  R2(2,1) = 0.0;
  R2(2,2) = R2(0,0);

  Matrix3d R3;
  R3(0,0) = cos(rpy[2]);
  R3(0,1) = -sin(rpy[2]);
  R3(0,2) = 0.0;
  R3(1,0) =-R3(0,1);
  R3(1,1) = R3(0,0);
  R3(1,2) = 0.0;
  R3(2,0) = 0.0;
  R3(2,1) = 0.0;
  R3(2,2) = 1.0;

  Matrix3d rotation = R3 * R2 * R1;

  rotation.col(0) = rotation.col(0) / rotation.col(0).norm();
  rotation.col(2) = rotation.col(0).cross(rotation.col(1));
  rotation.col(2) = rotation.col(2) / rotation.col(2).norm();
  rotation.col(1) = rotation.col(2).cross(rotation.col(0));
  rotation.col(1) = rotation.col(1) / rotation.col(1).norm();

  return rotation;
}

Eigen::Matrix3d
opengv::generateRandomRotation()
{
  Eigen::Vector3d rpy;
  rpy[0] = ((double) rand())/ ((double) RAND_MAX);
  rpy[1] = ((double) rand())/ ((double) RAND_MAX);
  rpy[2] = ((double) rand())/ ((double) RAND_MAX);

  rpy[0] = 2*M_PI*(rpy[0]-0.5);
  rpy[1] = M_PI*(rpy[1]-0.5);
  rpy[2] = 2*M_PI*(rpy[2]-0.5);

  Matrix3d R1;
  R1(0,0) = 1.0;
  R1(0,1) = 0.0;
  R1(0,2) = 0.0;
  R1(1,0) = 0.0;
  R1(1,1) = cos(rpy[0]);
  R1(1,2) = -sin(rpy[0]);
  R1(2,0) = 0.0;
  R1(2,1) = -R1(1,2);
  R1(2,2) = R1(1,1);

  Matrix3d R2;
  R2(0,0) = cos(rpy[1]);
  R2(0,1) = 0.0;
  R2(0,2) = sin(rpy[1]);
  R2(1,0) = 0.0;
  R2(1,1) = 1.0;
  R2(1,2) = 0.0;
  R2(2,0) = -R2(0,2);
  R2(2,1) = 0.0;
  R2(2,2) = R2(0,0);

  Matrix3d R3;
  R3(0,0) = cos(rpy[2]);
  R3(0,1) = -sin(rpy[2]);
  R3(0,2) = 0.0;
  R3(1,0) =-R3(0,1);
  R3(1,1) = R3(0,0);
  R3(1,2) = 0.0;
  R3(2,0) = 0.0;
  R3(2,1) = 0.0;
  R3(2,2) = 1.0;

  Matrix3d rotation = R3 * R2 * R1;

  rotation.col(0) = rotation.col(0) / rotation.col(0).norm();
  rotation.col(2) = rotation.col(0).cross(rotation.col(1));
  rotation.col(2) = rotation.col(2) / rotation.col(2).norm();
  rotation.col(1) = rotation.col(2).cross(rotation.col(0));
  rotation.col(1) = rotation.col(1) / rotation.col(1).norm();

  return rotation;
}
