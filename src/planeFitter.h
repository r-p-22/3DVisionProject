/*
 * PlaneFitter.h
 *
 *  Created on: Mar 10, 2016
 *      Author: lianos91
 */

#ifndef PLANEFITTER_H_
#define PLANEFITTER_H_

#include <Eigen/Dense>
#include <vector>


/**
 * \class planeFitter
 *
 *
 * This class executes RANSAC for the calculation of a plane to fit a number of points.
 *
 * 
 */

class PlaneFitter {

private:

	const float RANSAC_THRESH = 0.006;

	Eigen::Vector4d* fittedplane;

	Eigen::Matrix<double,4,Eigen::Dynamic>* Dinliers;


public:

	int numberOfInliers;
	PlaneFitter();

	virtual ~PlaneFitter();

	/*! 
   	 * getter method to retrieve the fitted plane ( in a Vector4d).
	 *
	*/
	Eigen::Vector4d getFittedPlane() const;

	/*! 
   	 * getter method to retrieve the inlier points after fitting a lattice.
	 *
	*/
	Eigen::Matrix<double,4,Eigen::Dynamic> getInlierPoints() const;

	/*! 
   	 * method to execute the ransac and fit a plane.
	 * @param[in] points3d the 3d points to fit a plane
	 * @param[in] inputIndices the indices of the input points (indices pointing the initial dataset, inputManager.pointModel)
	 *
	 * @return a vector of integers, containing the indices of the inliers (indices pointing the initial dataset, inputManager.pointModel)
	*/
	std::vector<int> ransacFit(std::vector<Eigen::Vector3d> points3d, std::vector<int> inputIndices);

	void fitPlane(Eigen::Matrix<double,Eigen::Dynamic,4> Points, Eigen::Vector4d *plane);



	/*! 
   	 * getter method to get the inliers 3D points, projected into the fitted plane.
	 *
	*/
	std::vector<Eigen::Vector3d> getProjectedInliers();

	static void vecToEigenMat(std::vector<Eigen::Vector3d> in,
			Eigen::Matrix<double,4,Eigen::Dynamic> &out, std::vector<int> idvec = std::vector<int>()){

		if (idvec.size()==0){
			for (size_t i = 0;i < in.size();i++){
				out.block<3,1>(0,i) = in.at(i);
				out(3,i) = 1.;
			}
		}else{
			for (size_t i = 0;i < idvec.size();i++){
				out.block<3,1>(0,i) = in.at(idvec[i]);
				out(3,i) = 1.;
			}
		}
		return;
	}


};

#endif /* PLANEFITTER_H_ */
