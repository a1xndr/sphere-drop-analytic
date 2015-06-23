CC=g++
CFLAGS= -I. 
all: sphere-drop-analytical
sphere-drop-analytical: sphere-drop-analytical.o sphere-math.o
	$(CC) sphere-drop-analytical.o sphere-math.o -o sphere-drop-analytical -I. -std=c++11
sphere-drop-analytical.o : sphere-drop-analytical.cpp
	$(CC) -c sphere-drop-analytical.cpp  -std=c++11
sphere-math.o : sphere-math.cpp
	$(CC) -c sphere-math.cpp
clean:
	rm *.o
