# 3DVisionProject

## libraries

Ceres solver 1.11.0
Eigen 3.2.0
opencv 2.4.12
CImg.h as included in the src/ folder

## running project

1. Clone the project to a local folder. `git clone https://github.com/ryenelith/3DVisionProject localDirectory/`
2. Open the command window and navigate to project folder.
3. Open CMakeLists.txt and change the directory paths in the following lines

	# external libs
	find_package(Ceres REQUIRED)
	include_directories(${CERES_INCLUDE_DIRS})
	#include_directories("/home/ceres-solver-1.11.0/ceres-bin/lib")
	INCLUDE_DIRECTORIES ( "/usr/include/eigen3" )

	#SET(CERES_LIBRARIES,
	#/home/ceres-solver-1.11.0/ceres-bin/lib/libceres.a
	#)

   such that they point to your local Ceres and Eigen files.

4. Run `cmake CMakeLists.txt`
5. Make the project by running the generated makefile using `make`
6. Add all images to data/images folder 
7. Run project using `./latt_bal [PATH_TO]/images.txt data/points3667_10views.txt [PATH_TO]/model-1-cams.txt [PATH_TO]/K.txt` with [PATH_TO] begin the path to Andrea's data
