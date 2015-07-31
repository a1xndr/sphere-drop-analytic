/*
 * =====================================================================================
 *
 *       Filename:  main.cpp
 *
 *    Description:  Sphere Packing Algorithm
 *
 *        Version:  1.0
 *        Created:  03/09/2015 04:02:39 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexander Oleinik
 *   Organization:  
 *
 * =====================================================================================
 */

#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <iostream>
 
const int NUM_SPHERES = 2000000;
const int X_MAX = 10;
const int Y_MAX = 10;
const int Z_MAX = 10;
const float R_MIN = 0.2;
const float R_MAX = 2;

struct position {
	float x, y, z;
};
struct sphere {
	float radius;
	position pos;
};

struct sphere spheres[NUM_SPHERES];

float rand_range(float min, float max)
{
    return min + static_cast<float>( rand()/static_cast<double>(RAND_MAX) * (max - min) );
}


int main()
{
	srand(time(NULL));
	double volume=0;
	for(int i=0; i < NUM_SPHERES; i++)
	{
		bool placed=false;
		int tries = 0;
		
		while(!placed && tries<1000000)
		{
			tries++;
			bool fail = false;
			
			spheres[i].radius = rand_range(R_MIN, R_MAX);
			float x = rand_range(spheres[i].radius, X_MAX - spheres[i].radius);
			float y = rand_range(spheres[i].radius, Y_MAX - spheres[i].radius);
			float z = rand_range(spheres[i].radius, Z_MAX - spheres[i].radius);
			//spheres[i].radius = (rand() + R_MIN) % R_MAX;
			//int x = (rand() + spheres[i].radius) % (X_MAX-spheres[i].radius); 
			//int y = (rand() + spheres[i].radius) % (Y_MAX-spheres[i].radius); 
			//int z = (rand() + spheres[i].radius) % (Z_MAX-spheres[i].radius); 
			for(int j=0; j < i; j++)
			{
				if( pow(spheres[j].pos.x - x , 2) + 
				    pow(spheres[j].pos.y - y , 2) + 
				    pow(spheres[j].pos.z - z , 2) 
				    < pow(spheres[j].radius + spheres[i].radius , 2))
				{
					fail = true;
					//std::cout << "ERR::PLACEMENT // Placing sphere" << i << " failed. attempt: "<< tries << std::endl;
					break;
				}
			}
			if(!fail)
			{
				placed = true;
				spheres[i].pos.x = x;
				spheres[i].pos.y = y;
				spheres[i].pos.z = z;
			}
		}
		if(!placed)
		{
			//std::cout << "ERR::PLACEMENT // Placing sphere" << i << " failed" << std::endl;
			break;
		}
		volume = volume + (4/3) * 3.14159 * pow(spheres[i].radius, 3);
		/*  std::cout << "INF::PLACEMENT // Placed sphere " << i 
		    << " of radius " << spheres[i].radius << " at coordinate ( " 
		    << spheres[i].pos.x << " , " << spheres[i].pos.y << " , " 
		    << spheres[i].pos.z << " )" << std::endl; */
		std::cout <<  spheres[i].radius << " " << spheres[i].pos.x <<  " " <<spheres[i].pos.y <<  " " <<spheres[i].pos.z << std::endl; 

	}
	std::cout << "INF::VOLUME // Total Container Volume is " << X_MAX*Y_MAX*Z_MAX << std::endl;
	std::cout << "INF::VOLUME // Total Sphere Volume is " << volume << std::endl;
	std::cout << "INF::VOLUME // Occupied Volume Ratio is " << volume/(X_MAX*Y_MAX*Z_MAX) << std::endl;
	return 0;
}
