

Useful notes:

1. The lattice fitting works for each group independently.

2. The lattice fitting depends on RANSAC. The ransac is executed on a small sample of distorted points, so for some groups it can detect the wrong lattice. Rerunning the lattice fitting is required in that case. The thresholds in Planefitter.h are hand tuned for best performance already.

3. The file 3dtools.h contains the method compareSiftFronto(), which is called during the lattice fitting. It checks whether two points (the referecePoint and the pointToTest) have similar SIFT descriptors. The camera selection was done by hard-coding the ones we took into consideration for the lattice fitting. For larger collection of images, the method of the paper for computing the most frontoparallel view requires refinenement (issues with occlusions and large angles in z-direction).

4. The images folder must be in the path ${executable_dir}/data. This can be easily changed in the 3dtools.h file, in the function computeSIFT.

5. The points.txt files (containing the point information) used were: (*)
	The lattices were generated with data/points3667.txt 
	For the augmented BA the data/points3667_10views.txt
All other input files (K.txt, images.txt,images folder, model-1-cams.txt) stay the same.

6. By the paper's method, if during the grouping, multiple points in different planes get grouped together (e.g. points in a window from perpendicular faces of a building) then only one lattice will be fitted. It will not be able to separate the planes accordingly.

7. Our data set was too small to establish a proper plane equality check, as needed for lattice consolidation. We therefore omitted the check and only evaluated the similarity of the basis vector sets of two lattices. In our case that was enough to get the desired result, because there were no two parallel faces with lattices on them. For bigger data sets, the plane check, at the moment outcommented in method calculateLatticeTransformation(...) in latticeClass.h, should be un-commented and refined.

8. The meaning of the field "consolidationTransformation" of the class LatticeClass is explained in the file BasisVectorTransformations.pdf . All lattices in a consolidation class get the proper index to describe their basis vector configuration relative to a common reference, which makes them comparable to each other. This should help to understand the code parts related to transformations. In general, the transformations are needed to make sure the similarity that should be enforced by the bundle adjustment can be defined in the proper way. It is necessary because the choice of one of the corner points to represent the lattice corner and the naming of the two basisvectors is arbitrary.

9. While experimenting, we tried out an approach where instead of using an additional error term for basis vector deviations between the lattices of a consolidation group, we used the very same two basis vectors for all of them in the optimization. This condition turned out to be too harsh for getting good results, so instead we used separate basis vector pairs and added the additional error term. The code, however, still contains support also for the first approach, which is referred to as "rigid". Look for that keyword when inspecting the BundleOptimizer class to find the appropriate methods if you want to try out that approach by yourself.


(*)	The groups saved in points.txt were generated at an early state of the project. All tests for lattice fitting and evaluations of the augmented bundle adjustment were done based on this initial grouping. During the project the grouping was further developed, the new group points were not used though, in order to make the development of later stages of the pipeline independet of the changing grouping files. When rerunning the whole project, slightly different results are to be expected due to an improved grouping of points.
