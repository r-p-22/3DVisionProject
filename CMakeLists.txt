CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(latt_bal)

#using c++11
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

# external libs
find_package(Ceres REQUIRED)
include_directories(${CERES_INCLUDE_DIRS})
#include_directories("/home/ceres-solver-1.11.0/ceres-bin/lib")
INCLUDE_DIRECTORIES ( "/usr/include/eigen3" )

#SET(CERES_LIBRARIES,
#/home/ceres-solver-1.11.0/ceres-bin/lib/libceres.a
#)

SET(SRC
src/3dtools.cpp
src/3dtools.h 
src/BundleOptimizer.cpp        
src/BundleOptimizer.h
src/camera.h                 
src/ceresReprojectionErrors.h 
src/CImg.h 
src/detectRepPoints.cpp 
src/detectRepPoints.h 
src/inputManager.h
src/latticeClass.h 
src/latticeDetector.cpp 
src/latticeDetector.h
src/latticeStruct.h
src/main.cpp        
src/my_v3d_vrmlio.h
src/planeFitter.cpp
src/planeFitter.h
)

SET(LINKFLAGS
pthread 
X11 
opencv_contrib 
opencv_core 
opencv_features2d opencv_highgui opencv_imgproc opencv_legacy opencv_nonfree 
)

ADD_EXECUTABLE(latt_bal ${SRC})
TARGET_LINK_LIBRARIES(latt_bal ${CERES_LIBRARIES} ${LINKFLAGS})
set (CMAKE_BUILD_TYPE Debug)

