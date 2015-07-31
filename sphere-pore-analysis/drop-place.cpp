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
 
//Compile-time constants
#define PI 3.14159265

//parameters
const int NUM_SPHERES = 2000000;
const int X_MAX = 10;
const int Y_MAX = 10;
const int Z_MAX = 10;
const float R_MIN = 0.2;
const float R_MAX = 2;
const float STEP = 0.01;

//Sphere definition object
struct position {
	float x, y, z;
};
struct sphere {
	float radius;
	position pos;
};

//Array of spheres to be filled
struct sphere spheres[NUM_SPHERES];


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  rand_range
 *  Description:  Choose a random float bewtween min and max. Make sure srand is called 
 *  with a seed before using this function.
 * =====================================================================================
 */
float rand_range(float min, float max)
{
    return min + static_cast<float>( rand()/static_cast<double>(RAND_MAX) * (max - min));
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  check_intersection
 *  Description:  Compares a theoretical sphere with radius=radius and center coordinate
 *  (x,y,z) against the first i spheres in the sphere array. Efficiency of this function
 *  is extrememly important and there is room for improvement as making a few addition
 *  subtraction checks could save 4, square operations. The pow() implementation is also
 *  slower than performing multiplications but is better for readibility. There are mult-
 *  iple similar distance calculations used in the program to which this also applies. 
 * =====================================================================================
 */
int check_intersection(float radius, float x, float y, float z, float i){
    for(int j=0; j < i; j++)
    {
	//Square of distance between two venters comparison against square of radius.
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

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  dot_product
 *  Description:  Dot product of the first n indices of two vectors
 * =====================================================================================
 */
float dot_product(float v[], float u[], int n)
{
    float result = 0.0;
    for (int i = 0; i < n; i++)
        result += v[i]*u[i];
    return result;
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  free_fall
 *  Description:  Free-fall is the default mode the currently simulated sphere adheres to
 *  and searches for the first collision along a straight vertical line. Includes a modi-
 *  fied check_intersection which check whether it is possible for the sphere in vertical
 *  freefall to collide with another by comparing x and y coordinates and picking the 
 *  collision with a sphere of largest z.
 * =====================================================================================
 */
int free_fall(sphere s, int i){
    
    std::cout   << "INF::MODE // Sphere: " << i << " entering free_fall mode" <<std::endl; 
    float z_max = 0;
    int c_index = -1;
    for(int j=0; j < i; j++)
    {
	if( pow(spheres[j].pos.x - s.pos.x , 2) + 
	    pow(spheres[j].pos.y - s.pos.y , 2) 
	    < pow(spheres[j].radius + s.radius , 2) && s.pos.z > spheres[j].pos.z)
	{
	    if(spheres[j].pos.z > z_max)
	    {
		z_max = spheres[j].pos.z;
		c_index = j;	
	    }
	}
    }
    //If there is a collision, move the sphere to the correct location. Else: move the 
    //sphere to the bottom of the container.
    if(c_index!=-1)s.pos.z= pow(pow(spheres[c_index].radius + s.radius , 2) 
		- pow(spheres[c_index].pos.x - s.pos.x , 2) - 
		pow(spheres[c_index].pos.y - s.pos.y , 2), 0.5);
    else s.pos.z = s.radius;
    return c_index;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  single_sphere_roll
 *  Description:  Roll on a single sphere until a second collision or until the polar
 *  angle is greater than pi/2 
 * =====================================================================================
 */
int single_sphere_roll(sphere s,int j, int i){
    std::cout   << "INF::MODE // Sphere: " << i << " entering single_sphere_roll mode" <<std::endl; 
    float r_traj = s.radius + spheres[j].radius;
    float phi = acos((s.pos.z-spheres[j].pos.z)/r_traj);
    float theta = acos((s.pos.y-spheres[j].pos.y)/(s.pos.x-spheres[j].pos.x));
    if(phi>0)
    {
	for(float i=phi; i<=PI/2; i+=STEP)
	{
	    float new_x, new_y, new_z;
	    new_x = s.pos.x + r_traj * cos(theta) * sin(i);
	    new_y = s.pos.y + r_traj * sin(theta) * sin(i);
	    new_z = s.pos.z + r_traj * cos(i);
	    int intersect = check_intersection(s.radius, new_x, new_y, new_z, i);
	    if(intersect != -1) return intersect;
	    else{
		s.pos.x = new_x;
		s.pos.y = new_y;
		s.pos.z = new_z;
	    }

	}
    }    
    else
    {
	for(float i=phi; i>=-PI/2; i-=STEP)
	{
	    float new_x, new_y, new_z;
	    new_x = s.pos.x + r_traj * cos(theta) * sin(i);
	    new_y = s.pos.y + r_traj * sin(theta) * sin(i);
	    new_z = s.pos.z + r_traj * cos(i);
	    int intersect = check_intersection(s.radius, new_x, new_y, new_z, i);
	    if(intersect != -1) return intersect;
	    else{
		s.pos.x = new_x;
		s.pos.y = new_y;
		s.pos.z = new_z;
	    }

	}
    }
    return -1;   
}
 
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  double_sphere_roll
 *  Description:  Roll sphere over two others until a third collision or until it reaches
 *  critical point and returns to free-fall mode. This step is very heavy on matrix oper-
 *  ations and is without a doubt the most expensive simulation mode.
 * =====================================================================================
 */
int double_sphere_roll(sphere s,int j,int k, int i){
    std::cout   << "INF::MODE // Sphere: " << i << " entering double_sphere_roll mode" << std::endl; 
    //Vector between two spheres of contact
    float jk[3] = { spheres[k].pos.x - spheres[j].pos.x,
		    spheres[k].pos.y - spheres[j].pos.y,
		    spheres[k].pos.z - spheres[j].pos.z};
    
    //Magnitude of the jk vector
    float jk_mag = pow(pow(jk[0],2) + pow(jk[1], 2) + pow(jk[2],2),0.5);
    
    //Normalized version of jk vector
    float jk_norm[3] = {jk[0]/jk_mag , jk[1]/jk_mag, jk[2]/jk_mag}; 
    
    //Vector between the simulated sphere and one of the contact spheres
    float js[3] = {s.pos.x - spheres[j].pos.x,
				s.pos.y - spheres[j].pos.y,
				s.pos.z - spheres[j].pos.z};

    //js dotted against jk
    float js_dot_jk = dot_product(js, jk, 3);
    
    //center of the circular trajectory of simulation
    position c_traj = {	spheres[j].pos.x + jk_norm[0]*js_dot_jk,
			spheres[j].pos.y + jk_norm[1]*js_dot_jk,
			spheres[j].pos.z + jk_norm[2]*js_dot_jk};
    //raduys if he trajectory
    float r_traj = pow( pow(spheres[j].radius + s.radius,2 ) 
			- pow(js_dot_jk,2), 0.5);
    //float phi = acos((s.pos.z-c_traj.z)/r_traj);
    //float theta = acos((s.pos.y-c_traj.y)/(s.pos.x-c_traj.x));

    //relative position of s where c_traj is 0
    position rel_pos = {s.pos.x - c_traj.x,
			s.pos.y - c_traj.y,
			s.pos.z - c_traj.z};

    bool cont = true;
    double angle;
    //step iteration/collision checking loop
    while(cont){
	angle+=STEP;

	//Apply Rotation Matrix described in 
	//http://inside.mines.edu/fs_home/gmurray/ArbitraryAxisRotation/  
	float new_x = c_traj.x + jk_norm[0]*(  jk_norm[0]*rel_pos.x + 
			    jk_norm[1]*rel_pos.y + 
			    jk_norm[2]*rel_pos.z)*(1-cos(angle)) 
		+ rel_pos.x * cos(angle) + 
		(-jk_norm[2]*rel_pos.y + jk_norm[1]*rel_pos.z)*(sin(angle));

	float new_y = c_traj.y + jk_norm[1]*(  jk_norm[0]*rel_pos.x + 
			    jk_norm[1]*rel_pos.y + 
			    jk_norm[2]*rel_pos.z)*(1-cos(angle)) 
		+ rel_pos.y * cos(angle) + 
		(jk_norm[2]*rel_pos.x - jk_norm[0]*rel_pos.z)*(sin(angle));

	float new_z = c_traj.z + jk_norm[2]*(  jk_norm[0]*rel_pos.x + 
			    jk_norm[1]*rel_pos.y + 
			    jk_norm[2]*rel_pos.z)*(1-cos(angle)) 
		+ rel_pos.z * cos(angle) + 
		(-jk_norm[1]*rel_pos.x + jk_norm[0]*rel_pos.y)*(sin(angle));

	//Make sure the sphere isnt trying to roll defying "gravity". This is done by 
	//making sure that the distance between X and Y of s and j/k isnt decreasing.
	if(pow(new_x-spheres[j].pos.x, 2)+pow(new_y-spheres[j].pos.y, 2) < pow(s.pos.x-spheres[j].pos.x,2)+pow(s.pos.y-spheres[j].pos.y,2))break;
	if(pow(new_x-spheres[k].pos.x, 2)+pow(new_y-spheres[k].pos.y, 2) < pow(s.pos.x-spheres[k].pos.x,2)+pow(s.pos.y-spheres[k].pos.y,2))break;

	//check that there are no intersections. 
	int intersect = check_intersection(s.radius, new_x, new_y, new_z, i);
	
	if(intersect !=-1) return intersect;
	//After all these transformations and checks, change the coordinates and continue the loop
	else{
	    s.pos.x = new_x;
	    s.pos.y = new_y;
	    s.pos.z = new_z;
	}
    }

    //back to freefall
    return 0;
}

int main()
{
	//define random seed
	//srand(time(NULL));
	srand(123456);

	//total volume occupied by spheres
	double volume=0;
	
	//Loop to attempt placmente of NUM_SPHERES spheres
	for(int i=0; i < NUM_SPHERES; i++)
	{	
		bool placed=false;
		int tries = 0;
		//In this simulation "tries" keeps a count of all the attempts to place a 
		//sphere without immediate intersection(into a valid location) at the top
		//of a container. If this fails many times, give up.
		while(!placed && tries<100000)
		{
			tries++;
			bool posfail = false;
			bool radiusfail = false;

			//Pick random radius and x,y	
			float radius = rand_range(R_MIN, R_MAX);;	
			float x = rand_range(radius, X_MAX - radius);
			float y = rand_range(radius, Y_MAX - radius);
			float z = Z_MAX - radius;

			//make sure there is no immediate intersection
			if(check_intersection(radius,x,y,z,i)!=-1)
			{
			    posfail=true;
			}
			if(!posfail)
			{
			    sphere s = {radius,{x,y,z}};
			    bool lodged = false;
			    int state = 0;
			    int contact[3]={-1,-1,-1};
			    int collision;

			    //Logic to iterate between modes by interpreting the
			    //return codes
			    while(!lodged)
			    {
				switch(state)
				{
				    case 0:
					contact[0]=free_fall(s, i);
					state++;
					break;
				    case 1:
					collision=single_sphere_roll(s, contact[0], i);
					if(collision!=-1){
					    state=2;
					    contact[1]=collision;
					}
					else{
					    contact[0]=-1;
					} 
					break;
				    case 2:
					collision=double_sphere_roll(s, contact[0],contact[1], i);
					if(collision!=-1){
					    state=3;
					    contact[2]=collision;
					}
					else{
					    contact[0]=-1;
					    contact[1]=-1;
					} 
					break;
				    case 3:
					lodged=true;
					break;
				}	    
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
			std::cout << "ERR::PLACEMENT // Placing sphere" << i << " failed" << std::endl;
			break;
		}
		volume = volume + (4/3) * 3.14159 * pow(spheres[i].radius, 3);
		/*  std::cout << "INF::PLACEMENT // Placed sphere " << i 
		    << " of radius " << spheres[i].radius << " at coordinate ( " 
		    << spheres[i].pos.x << " , " << spheres[i].pos.y << " , " 
		    << spheres[i].pos.z << " )" << std::endl; */
		std::cout   <<  spheres[i].radius << " " 
			    << spheres[i].pos.x 
			    <<  " " <<spheres[i].pos.y 
			    <<  " " <<spheres[i].pos.z << std::endl; 

	}
	std::cout   << "INF::VOLUME // Total Container Volume is " 
		    << X_MAX*Y_MAX*Z_MAX << std::endl;
	std::cout   << "INF::VOLUME // Total Sphere Volume is " 
		    << volume << std::endl;
	std::cout   << "INF::VOLUME // Occupied Volume Ratio is " 
		    << volume/(X_MAX*Y_MAX*Z_MAX) << std::endl;
	return 0;
}

