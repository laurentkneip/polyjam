#include "experiment_helpers.hpp"
#include "random_generators.hpp"
#include <iostream>
#include <iomanip>


void
opengv::generateRandomCameraSystem(
    int numberCameras,
    std::vector<Eigen::Vector3d> & camOffsets,
    std::vector<Eigen::Matrix3d> & camRotations )
{
  double offset = 0.5; //this is the distance from the viewpoint origin
  
  for( int i = 0; i < numberCameras; i++ )
  {
    Eigen::Vector3d camOffset = generateRandomDirectionTranslation(offset);
    Eigen::Matrix3d camRotation = generateRandomRotation();
    camOffsets.push_back(camOffset);
    camRotations.push_back(camRotation);
  }
}

void
opengv::printExperimentCharacteristics(
    const Eigen::Vector3d & position,
    const Eigen::Matrix3d & rotation,
    double noise,
    double outlierFraction )
{
  std::cout << "the random position is:" << std::endl;
  std::cout << position << std::endl << std::endl;
  std::cout << "the random rotation is:" << std::endl;
  std::cout << rotation << std::endl << std::endl;
  std::cout << "the noise in the data is:" << std::endl;
  std::cout << noise << std::endl;
  std::cout << "the outlier fraction is:" << std::endl;
  std::cout << outlierFraction << std::endl;
}

void
opengv::generateRandom2D3DCorrespondences(
    const Eigen::Vector3d & position,
    const Eigen::Matrix3d & rotation,
    const std::vector<Eigen::Vector3d> & camOffsets,
    const std::vector<Eigen::Matrix3d> & camRotations,
    size_t numberPoints,
    double noise,
    double outlierFraction,
    std::vector<Eigen::Vector3d> & bearingVectors,
    std::vector<Eigen::Vector3d> & points,
    std::vector<int> & camCorrespondences,
    Eigen::MatrixXd & gt )
{
  //initialize point-cloud
  double minDepth = 4;
  double maxDepth = 8;
  
  for( size_t i = 0; i < (size_t) gt.cols(); i++ )
    gt.col(i) = generateRandomPoint( maxDepth, minDepth );
    
  //create the 2D3D-correspondences by looping through the cameras
  size_t numberCams = camOffsets.size();
  size_t camCorrespondence = 0;
  
  for( size_t i = 0; i < (size_t) gt.cols(); i++ )
  {
    //get the camera transformation
    Eigen::Vector3d camOffset = camOffsets[camCorrespondence];
    Eigen::Matrix3d camRotation = camRotations[camCorrespondence];

    //store the point
    points.push_back(gt.col(i));
    
    //project the point into the viewpoint frame
    Eigen::Vector3d bodyPoint = rotation.transpose()*(gt.col(i) - position);
    
    //project the point into the camera frame
    bearingVectors.push_back(camRotation.transpose()*(bodyPoint - camOffset));

    //normalize the bearing-vector to 1
    bearingVectors[i] = bearingVectors[i] / bearingVectors[i].norm();

    //add noise
    if( noise > 0.0 )
      bearingVectors[i] = addNoise(noise,bearingVectors[i]);
    
    //push back the camera correspondence
    camCorrespondences.push_back(camCorrespondence++);
    if(camCorrespondence > (numberCams-1) )
      camCorrespondence = 0;  
  }
  
  //add outliers
  //compute the number of outliers based on fraction
  size_t numberOutliers = (size_t) floor(outlierFraction*numberPoints);
  //make the first numberOutliers points be outliers
  for(size_t i = 0; i < numberOutliers; i++)
  {
    //extract the camera transformation
    Eigen::Vector3d camOffset = camOffsets[camCorrespondences[i]];
    Eigen::Matrix3d camRotation = camRotations[camCorrespondences[i]];

    //create random point
    Eigen::Vector3d p = generateRandomPoint(8,4);
    
    //project into viewpoint frame
    Eigen::Vector3d bodyPoint = rotation.transpose()*(p - position);
    
    //project into camera-frame and use as outlier measurement
    bearingVectors[i] = camRotation.transpose()*(bodyPoint - camOffset);
    
    //normalize the bearing vector
    bearingVectors[i] = bearingVectors[i] / bearingVectors[i].norm();
  }
}