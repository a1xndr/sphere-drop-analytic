CC=g++
CFLAGS= -I. 
all: sphere-drop-analytical
sphere-drop-analytical: sphere-drop-analytical.o sphere_math.o
	$(CC) sphere-drop-analytical.o sphere_math.o -o sphere-drop-analytical -I. -std=c++11
sphere-drop-analytical.o : sphere-drop-analytical.cpp
	$(CC) -c sphere-drop-analytical.cpp  -std=c++11
sphere_math.o : sphere_math.cpp
	$(CC) -c sphere_math.cpp
clean:
	rm *.o
