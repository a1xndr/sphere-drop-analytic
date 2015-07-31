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
const float R_STEP = 0.01;
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
// R X Y Z maximum index 
int check_intersection(float radius, float x, float y, float z, float i){
    for(int j=0; j < i; j++)
    {
	if( pow(spheres[j].pos.x - x , 2) + 
	    pow(spheres[j].pos.y - y , 2) + 
	    pow(spheres[j].pos.z - z , 2) 
	    < pow(spheres[j].radius + radius , 2))
	{
	    return j;
	}
    }
    return -1;
}

int main()
{
	//srand(time(NULL));
	srand(123456);
	double volume=0;
	for(int i=0; i < NUM_SPHERES; i++)
	{
		bool placed=false;
		int tries = 0;
		while(!placed && tries<10000000)
		{
			tries++;
			bool posfail = false;
			bool radiusfail = false;
			
			float radius = R_MIN;	
			float x = rand_range(R_MIN, X_MAX - R_MIN);
			float y = rand_range(R_MIN, Y_MAX - R_MIN);
			float z = rand_range(R_MIN, Z_MAX - R_MIN);
			if(check_intersection(radius,x,y,z,i)!=-1)
			{
			    posfail=true;
			}
			if(!posfail)
			{
			    while(radius+R_STEP<R_MAX){
				if(x+radius+R_STEP>X_MAX || y+radius+R_STEP > Y_MAX || z+radius+R_STEP >Z_MAX || x-radius-R_STEP<0 || y-radius-R_STEP<0 || z-radius-R_STEP<0){
				    radiusfail=true;
				    break;
				}
				int collided_sphere=check_intersection(radius+R_STEP,x,y,z,i);
				if(collided_sphere!=-1)
				{
				    float vec_x = x - spheres[collided_sphere].pos.x;
				    float vec_y = y - spheres[collided_sphere].pos.y;
				    float vec_z = z - spheres[collided_sphere].pos.z;
				    float vec_length = pow(pow(vec_x,2)+ pow(vec_y, 2)+ pow(vec_z, 2),0.5);
				    float shifted_x = x + vec_x * R_STEP;
				    float shifted_y = y + vec_y * R_STEP;
				    float shifted_z = z + vec_z * R_STEP;
				    if(check_intersection(radius+R_STEP,shifted_x,shifted_y,shifted_z,i)!=-1){
					radiusfail=true;
				    }
				    else{
					if(shifted_x+radius+R_STEP>X_MAX || shifted_y+radius+R_STEP > Y_MAX || shifted_z+radius+R_STEP >Z_MAX || 
						shifted_x-radius-R_STEP<0 || shifted_y-radius-R_STEP<0 || shifted_z-radius-R_STEP<0){
					    radiusfail=true;
					}
					else{
					    x = shifted_x;
					    y = shifted_y;
					    z = shifted_z;
					}
				    }
				}
				if(radiusfail)break;
				radius+=R_STEP;
			    }
			    placed = true;
			    spheres[i].radius = radius;
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
