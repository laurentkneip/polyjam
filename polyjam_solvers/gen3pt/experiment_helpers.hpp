#ifndef OPENGV_EXPERIMENT_HELPERS_HPP_
#define OPENGV_EXPERIMENT_HELPERS_HPP_

#include <stdlib.h>
#include <vector>
#include <Eigen/Eigen>
#include <boost/shared_ptr.hpp>

namespace opengv
{

void generateRandomCameraSystem(
    int numberCameras,
    std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> > & camOffsets,
    std::vector<Eigen::Matrix3d,Eigen::aligned_allocator<Eigen::Matrix3d> > & camRotations );

void printExperimentCharacteristics(
    const Eigen::Vector3d & position,
    const Eigen::Matrix3d & rotation,
    double noise,
    double outlierFraction );

void generateRandom2D3DCorrespondences(
    const Eigen::Vector3d & position,
    const Eigen::Matrix3d & rotation,
    const std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> > & camOffsets,
    const std::vector<Eigen::Matrix3d, Eigen::aligned_allocator<Eigen::Matrix3d> > & camRotations,
    size_t numberPoints,
    double noise,
    double outlierFraction,
    std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> > & bearingVectors,
    std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> > & points,
    std::vector<int> & camCorrespondences,
    Eigen::MatrixXd & gt );

}

#endif /* OPENGV_EXPERIMENT_HELPERS_HPP_ */
