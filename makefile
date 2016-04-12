CFLAGS=\
-I./Libraries/eigen \
-I./Libraries/openCV/include \
-I./Libraries/3dTools \
-I./Libraries/CImg \

LINKFLAGS= -lopencv_calib3d -lopencv_contrib -lopencv_core -lopencv_features2d -lopencv_flann -lopencv_gpu -lopencv_highgui -lopencv_imgproc -lopencv_legacy -lopencv_ml -lopencv_nonfree -lopencv_objdetect -lopencv_ocl -lopencv_photo -lopencv_stitching -lopencv_superres -lopencv_ts -lopencv_video -lopencv_videostab

all: main.o detectRepPoints.o planeFit.o latticeDetector.o
	g++ $(LINKFLAGS) -o 3DVisionProject main.o detectRepPoints.o planeFit.o latticeDetector.o

main.o: src/main.cpp src/detectRepPoints.h
	g++ $(CFLAGS) -c src/main.cpp -o main.o

detectRepPoints.o: src/detectRepPoints.cpp src/detectRepPoints.h
	g++ $(CFLAGS) -c src/detectRepPoints.cpp -o detectRepPoints.o

planeFit.o: src/planeFit.cpp src/main.h
	g++ $(CFLAGS) -c src/planeFit.cpp -o planeFit.o

latticeDetector.o: src/latticeDetector.h src/latticeDetector.cpp
	g++ $(CFLAGS) -c src/latticeDetector.cpp -o latticeDetector.o

clear:
	rm *.o 3DVisionProject
