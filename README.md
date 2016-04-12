# 3DVisionProject

## running project
1. Clone the project to a local folder. `git clone https://github.com/ryenelith/3DVisionProject localDirectory/`
2. Open command window and navigate to project folder.
3. make the project by running the makefile using `make`
4. replace image_rep.txt with image.txt containing all image names (not just 5), also add all images to image/ folder 
5. run project using `./3DVisionProject data/images.txt data/points.txt`

Hope you get it working. It should 
- Load the sift features from a file (fast option),
- ouput the group members, 
- output 3d locations for given 3dpoint indexes for a reduced points.txt file (with less points for testing) and 
- run the planeFit function of Nektarios.

> This project requires an installation of openCV on your computer. If the makefile throws and error that's probably why.

> Step 4 is important to make it work when recomputing the sift features from images (because the points.txt file assumes they exist). If you choose to read the sift features from file you don't need the images/ folder or the images.txt file.

> For grouping, run initialise class with option to read from images `(argv,1)` once, then (now the siftFeatureVector.txt is writen) initialise with arguments `(argv,0)` to read sift features from file (faster). Every time the points.txt file is changed, initialise with `(argv,1)` once to update siftFeatureVector.txt.

## Git suggestions

> Download GitHub Desktop! It makes your life easier.

> All group members have write access to the project. This means everyone in the group can make branches, commit pull requests to the master, even edit the master directly online (not advisable...)

#### Instructions to work on project
0. First time downloading
  - clone master to local directory as mentioned in running project section

1. For each new topic (e.g. add class to detect lattices) open new branch from master: 
  - make a branch using naming convention: NAME_branch-TOPIC
  - make sure you are working on the branch and comparing to the master (not editing the master directly)
  - edit the files in your local directory and commit changes to your branch in gitHub desktop
  - when happy with the brach, start pull request to master
  - Switch to https://github.com/ryenelith/3DVisionProject (online) and commit to pull request

> Commiting a change to a branch and syncing with github will not change the master, but will change the branch for all group members. Github probably has a way to sort out issues when people work on the same branch and sync, but probably advisable not to work on other peoples branch...

> Pull requests to master should be discussed in the group. Especially when conflicts occur (due to two pull requests changing same file).

#### Workings of gitHub Desktop:

- commited to local directory changes to a branch on GitHub desktop will be merged with the master when you make a pull request
- `sync` button updates to show differences of current branch to master (doesn’t change local data)
- `update from master` as you go to get most recent updates from pull requests from other co-workers (changes local data, might require solving conflicts)
- Solving conflicts possible when
  - Using `update from master` before starting pull request or 
  - once two conflicting pullrequests have been uploaded (will just give a note -> go back to gitHub desktop and use `update from master`, manage conflicts)
- When a pull request of a branch to the master is commited, the branch is equal to the master and can be deleted (carefull with files included in gitignore)
- gitHub desktop keeps the files in your local directory linked with the branch. If you delete a branch and start a new branch from the master, any files included in the .gitignore (e.g. build directory etc.) will no longer be available in your local directory. 

> Keep copy of the files in .gitignore that you don't want to show up in changes to commit in the gitHub desktop, but that you plan to reuse in the next branch. In this project these are for example: /images folder and images.txt file and build folders.
