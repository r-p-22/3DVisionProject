//============================================================================
// Name        : planeFit.cpp
// Author      : Konstantinos Nektarios Lianos
// Version     :
// Copyright   : 
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <Eigen/Dense>
#include <vector>

#include "main.h"

using namespace std;

int planeFit() {

	int n = 4;

	//Eigen::Matrix<float,Eigen::Dynamic,4> D(n,4);
	Eigen::Matrix<float,3,4> D;
	D << 1,1,1,1,
	0,0,0,1,
	1,0,0,1;


	/*
	Insert points?..
	int i=0;
	for (int i=0;i<n;i++){
	Eigen::Vector4f p = new Eigen::Vector4f(3);
	}
	 */

	Eigen::JacobiSVD<Eigen::MatrixXf> svd(D, Eigen::ComputeThinU | Eigen::ComputeFullV);
	cout << "The eigenvalues of A are:\n" << svd.singularValues() << endl;
	cout << "Its left singular vectors are the columns of the thin U matrix:" << endl << svd.matrixU() << endl;
	cout << "Its right singular vectors are the columns of the full V matrix:" << endl << svd.matrixV() << endl;

	Eigen::Vector4f Plane = svd.matrixV().block<4,1>(0,3);//<sizeRows,sizeCols>(beginRow,beginCol)

	cout << "Its right singular vectors are the columns of the full V matrix:" << endl << Plane << endl;

	// u is not in the plane
	Eigen::Vector4f u(-1, -1, 1, 1);
	cout << u.transpose()*Plane << endl;

	return 0;

}
