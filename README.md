# 3DVisionProject

## running project
1. Clone the project to a local folder. `git clone https://github.com/ryenelith/3DVisionProject localDirectory/`
2. Open command window and navigate to project folder.
3. make the project by running the makefile using `make`
4. replace image_rep.txt with image.txt containing all image names (not just 5), also add all images to image/ folder 
5. run project using `./3DVisionProject data/images.txt data/points.txt data/grouping/outputPoints.txt data/grouping/outSiftFeaturesVector.txt`

Hope you get it working. It should 
- Load the sift features from a file (fast option),
- ouput the group members, 
- output 3d locations for given 3dpoint indexes for a reduced points.txt file (with less points for testing) and 
- run the planeFit function of Nektarios.

> Step 4 is important to make it work.

## Git

### Idea how it could work for our group

> Seems to be a good habit to have a separate local directory for each branch /topic

> GitHub Desktop makes your life easier

1. to open new branch: 
  - clone master to local directory 
  - make a branch using naming convention: NAME_branch-TOPIC

2. Edit in local directory
  - accept all changes in local directory on GitHub desktop you want to later add to master
  - sync button updates to show differences of current branch to master (doesn’t change local data)
  - update from master as you go to get most recent updates from pull requests from other co-workers (changes local data)

3. Sync changes of branch to master
  - Update from master to perhaps already solve conflicting cases locally.  
  - Open pull request of your branch
  - wait for pull request to be accepted by administrator of repo
  - master and branch are now identical -> delete branch and start new one for new topic
  - update from master (for meanwhile co-worker changes) 
  - continue working with local directory files
