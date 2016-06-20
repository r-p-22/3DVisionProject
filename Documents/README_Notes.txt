

Useful notes:

0. Run the program with the following paramenters

	[images.txt] [points.txt] [cams.txt] [K.txt]

(with the more specific file names).

1. The lattice fitting works for each group independently.

2. The lattice fitting depends on RANSAC. The ransac is executed on a small sample of distorted points, so for some groups it can detect the wrong lattice. Rerunning the lattice fitting is required in that case. The thresholds in Planefitter.h are hand tuned for best performance already.

3. The file 3dtools.h contains the method compareSiftFronto(), which is called during the lattice fitting. It checks whether two points (the referecePoint and the pointToTest) have similar SIFT descriptors. The camera selection was done by hard-coding the ones we took into consideration for the lattice fitting. For larger collection of images, the method of the paper for computing the most frontoparallel view requires refinenement (issues with occlusions and large angles in z-direction).

4. The images folder must be in the path ${executable_dir}/data. This can be easily changed in the 3dtools.h file, in the function computeSIFT.

5. The points.txt files (containing the point information) used were:
	The lattices were generated with data/points3667.txt 
	For the augmented BA the data/points3667_10views.txt
All other input files (K.txt, images.txt,images folder, model-1-cams.txt) stay the same.

6. By the paper's method, if during the grouping, multiple points in different planes get grouped together (e.g. points in a window from perpendicular faces of a building) then only one lattice will be fitted. It will not be able to separate the planes accordingly.

7. Our data set was too small to establish a proper plane equality check, as needed for lattice consolidation. We therefore omitted the check and only evaluated the similarity of the basis vector sets of two lattices. In our case that was enough to get the desired result, because there were no two parallel faces with lattices on them. For bigger data sets, the plane check, at the moment outcommented in method calculateLatticeTransformation(...) in latticeClass.h, should be un-commented and refined.
