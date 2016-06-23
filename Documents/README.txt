# 3DVisionProject

## libraries

Ceres solver 1.11.0 	- http://ceres-solver.org/ceres-solver-1.11.0.tar.gz
Eigen 3.2.0 		- https://launchpad.net/ubuntu/+source/eigen3/3.2.0-8
opencv 2.4.12 		- https://github.com/Itseez/opencv/archive/2.4.12.zip
CImg.h as included in the src/ folder

## files

- src/	Folder with all source files. main.cpp contains main function, which loads our pre-grouped point groups and our pre-computed, selected lattices from file and applies bundle adjustment to the model

- data/	Folder that contains some subfolders, together with some files used in the program - relevant are points3667.txt and points3667_10views.txt, on which we ran our program when generating the final data (see file README_Notes.txt) for details

- data/distanceVectors	Folder that contains some of the generated results for the distance vectors in different configurations. data/distanceVectors/Trivial loss/distanceVectorPlot contains the MATLAB script and used result files to generate the plots shown in the report and also present in said folder

- data/grouping		Folder that contains the output of different runs of point grouping (different model size, different parameters etc). Relevant is the file outputPoints.txt, which was used when generating our results.

- data/savedLattices	Folder that contains saved versions of detected lattices. We manually selected a subset of them, that are loaded from file in the current version of the code. See src/main.cpp for details.	

## what, apart from the mentioned libraries, we did not code ourselves, but took from Andrea

- my_v3d_vrmlio.h (lattice and plane visualizers are ours though)
- camera.h
- 3dtools.cpp and 3dtools.h
- structs: TriangulatedPoint, PointMeasurement

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
