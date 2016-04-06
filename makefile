CFLAGS=\
-I./Libraries/eigen \
-I./Libraries/CImg \

all: main.o detectRepPoints.o planeFit.o latticeDetector.o
	g++ -o 3DVisionProject main.o detectRepPoints.o planeFit.o latticeDetector.o

main.o: src/main.cpp src/main.h src/detectRepPoints.h src/latticeDetector.h
	g++ $(CFLAGS) -c src/main.cpp -o main.o

detectRepPoints.o: src/detectRepPoints.cpp src/detectRepPoints.h
	g++ $(CFLAGS) -c src/detectRepPoints.cpp -o detectRepPoints.o

planeFit.o: src/planeFit.cpp src/main.h
	g++ $(CFLAGS) -c src/planeFit.cpp -o planeFit.o

latticeDetector.o: src/latticeDetector.cpp src/latticeDetector.h
	g++ $(CFLAGS) -c src/latticeDetector.cpp -o planeFit.o

clear:
	rm *.o 3DVisionProject