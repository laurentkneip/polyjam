#ifndef OPENGV_RANDOM_GENERATORS_HPP_
#define OPENGV_RANDOM_GENERATORS_HPP_

#include <stdlib.h>
#include <vector>
#include <Eigen/Eigen>

namespace opengv
{

void initializeRandomSeed();
Eigen::Vector3d generateRandomPoint( double maximumDepth, double minimumDepth );
Eigen::Vector3d generateRandomPointPlane();
Eigen::Vector3d addNoise( double noiseLevel, Eigen::Vector3d cleanPoint );
Eigen::Vector3d generateRandomTranslation( double maximumParallax );
Eigen::Vector3d generateRandomDirectionTranslation( double parallax );
Eigen::Matrix3d generateRandomRotation( double maxAngle );
Eigen::Matrix3d generateRandomRotation();

}

#endif /* OPENGV_RANDOM_GENERATORS_HPP_ */
