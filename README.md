# 3DVisionProject

## running project
- Clone the project to a local folder. `git clone https://github.com/ryenelith/3DVisionProject localDirectory/`
- Open command window and navigate to project folder.
- make the project by running the makefile using `make`
- run project using `./3DVisionProject data/cams.txt data/points99.txt`

Hope you get it working. It should 
- ouput the group members, 
- output 3d locations for given 3dpoint indexes for a reduced points.txt file (with less points for testing) and 
- run the planeFit function of Nektarios.

## group organisation

An idea how we could organise our group repository:

* To make edits to the master, clone project to local directory. 
  * If you are working on two topics in parallel, clone project to two different local directories
  * add a branch and use following naming convention for new branches: NAME_branch-TOPIC (e.g. naming ryen_branch-siftStuff)
* To merge your changes with the master, make new pull request. Probably it makes sense to use the group meetings to look at what new features the different branches have and merge to the master.
* Once the master has been updated:
  * Continue working on your current branch / topic or
  * create a new branch for a new topic by making a new clone to a local directory e.g. overright previous directory if all changes from previous topic have been merged with master. Or clone to different directory if changes of your current topic are still being discussed or you are working on two topics in parallel.
