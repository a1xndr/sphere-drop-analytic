/*
 * =====================================================================================
 *
 *       Filename:  sphere_drop.cpp
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
#include <vector>
#include <iostream>
#include "sphere_math.hpp"
#include "sphere.hpp"
#include <random>


//parameters
const int NUM_SPHERES = 4;
const int X_MAX = 10;
const int Y_MAX = 10;
const int Z_MAX = 10;
const double R_MIN = 1;
const double R_MAX = 1;
const double STEP = 0.001;


//Array of spheres to be filled
struct sphere spheres[NUM_SPHERES];

int check_intersect(double radius, vec3 pos, double i, bool fastx, bool fasty, bool fastz){
    for(int j=i-1; j >=0; j--)
    {
    	if(fastx && abs(pos.x-spheres[j].pos.x) > radius + spheres[j].radius)continue;
    	if(fasty && abs(pos.y-spheres[j].pos.y) > radius + spheres[j].radius)continue;
    	if(fastz && abs(pos.z-spheres[j].pos.z) > radius + spheres[j].radius)continue;
    	//Square of distance between two centers comparison against square of radius.
    	if( distance(pos,spheres[j].pos)
	    < spheres[j].radius + radius - SMIDGE)return j;
    }
    return -1;
}

std::vector<int> check_intersects(double radius, vec3 pos, double i, bool fastx, bool fasty, bool fastz){
    std::vector<int> intersects;
	for(int j=i-1; j >=0; j--)
    {
    	if(fastx && abs(pos.x-spheres[j].pos.x) > radius + spheres[j].radius)continue;
    	if(fasty && abs(pos.y-spheres[j].pos.y) > radius + spheres[j].radius)continue;
    	if(fastz && abs(pos.z-spheres[j].pos.z) > radius + spheres[j].radius)continue;
    	//Square of distance between two centers comparison against square of radius.
    	if( distance(pos,spheres[j].pos)
	    < spheres[j].radius + radius - SMIDGE)intersects.push_back(j);
    }
    return intersects;
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
int free_fall(sphere *s, int i){
    std::cout << "Entering FF" <<std::endl;
    double z_max = 0;
    int c_index = -1;
    for(int j=i-1; j>=0; j--)
    {
    	if(s->pos.z < spheres[j].pos.z )continue;
    	if(abs(s->pos.x-spheres[j].pos.x) > s->radius + spheres[j].radius)continue;
    	if(abs(s->pos.y-spheres[j].pos.y) > s->radius + spheres[j].radius)continue;

		if( distance(vec3(s->pos.x,s->pos.y,0), vec3(spheres[j].pos.x,spheres[j].pos.y,0))
			< spheres[j].radius + s->radius - SMIDGE)
		{
			double z_col=spheres[j].pos.z +pow(pow(spheres[j].radius +  s->radius , 2)
					- pow(spheres[j].pos.x -  s->pos.x , 2) -
					pow(spheres[j].pos.y -  s->pos.y , 2), 0.5);
			if(z_col > z_max && s->pos.z>z_max)
			{
				z_max = z_col;
				c_index = j;
			}
		}
    }
    //If there is a collision, move the sphere to the correct location. Else: move the 
    //sphere to the bottom of the container.
    if(c_index!=-1 && s->pos.z>z_max){
    	s->pos.z= z_max;
        return c_index;
    std::cout << "Hit against: " << c_index <<std::endl;
    }
    else {s->pos.z = s->radius;}
    return -1;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  single_sphere_roll
 *  Description:  Roll on a single sphere until a second collision or until the polar
 *  angle is greater than pi. Analytical version
 * =====================================================================================
 */
int single_sphere_roll(sphere *s,int j, int i){
	if(i==3){
		std::cout << "Entering SSR" <<std::endl;
	}	//Limit the amount of spheres we run the math against to the radius of the two
	//touching spheres.
    double r_collision_volume = 2*s->radius + spheres[j].radius;
    std::vector<int> intersects = check_intersects(r_collision_volume,
    		spheres[j].pos, i ,1,1,1);
    int collision = -1;
    //Assemble the spherical coordinates that we will need to run the toroidal collision math
    double r_traj 	= s->radius + spheres[j].radius;
    double phi 		= asin((s->pos.z-spheres[j].pos.z)/r_traj);
    double theta 	= atan2((s->pos.y-spheres[j].pos.y),(s->pos.x-spheres[j].pos.x));
    vec3 u( cos(theta),	sin(theta), 0);
    vec3 v(-sin(theta), cos(theta), 0);
    vec3 w(0,0,1);

    double maxsol = 1;
    double newphi=0;
    double T = cos(0);
    for(int l=0; l<intersects.size(); l++){
    	if(intersects[l]==j || intersects[l]==i)continue;
    	//Set up the variables exactly as detailed in the paper.
    	double Wx 	= (spheres[intersects[l]].pos.x-spheres[j].pos.x)/(s->radius+spheres[j].radius);
    	double Wy 	= (spheres[intersects[l]].pos.y-spheres[j].pos.y)/(s->radius+spheres[j].radius);
    	double Wz 	= (spheres[intersects[l]].pos.z-spheres[j].pos.z)/(s->radius+spheres[j].radius);
    	double W 	= (spheres[intersects[l]].radius+s->radius)/(s->radius+spheres[j].radius);

    	double K1 	= 2.0*(Wx*cos(theta)+Wy*sin(theta));
    	double K2 	= 2.0*Wz;
    	double K3 	= pow(Wx,2) + pow(Wy,2) + pow(Wz,2) +1 - pow(W,2);

    	//Now solve the quadratic equation
    	double D = pow(2*K1*K3,2)-4*(pow(K1,2)+pow(K2,2))*(pow(K3,2)-pow(K2,2));
    	if(D<0)continue;
    	double T1 = (2*K1*K3 + pow(D,0.5))/(2*(pow(K1,2)+pow(K2,2)));
    	double T2 = (2*K1*K3 - pow(D,0.5))/(2*(pow(K1,2)+pow(K2,2)));
    	if(K2*(K3-K1*T1)<0)T1=T2;
    	if(K2*(K3-K1*T2)<0)T2=T1;
    	if(K2*(K3-K1*T2)<0)continue;

	double sol = std::min(T1, T2);
    	double ang1 = acos(T1);
    	double ang2 = acos(T2);
    	double ang;

	if(sol<maxsol && sol>0 && sol > cos(phi)){
	    maxsol = sol;
	    collision = intersects[l];
	}
//    	if(spheres[intersects[l]].pos.z<spheres[j].pos.z){
//    		ang=std::min(ang1, ang2);
//    	}
//    	else{
//    		ang=std::max(ang1, ang2);
//    	}
//    	if(ang>newphi && ang<phi)
//    	{
//    		collision = intersects[l];
//        	newphi=ang;
//    	}
//    	newphi = ang;
    }
    newphi = acos(maxsol); 
    std::cout <<"phi is:" <<newphi <<std::endl;
    std::cout << "hit against: " << collision <<std::endl;
    s->pos.x = spheres[j].pos.x + r_traj*cos(theta)*cos(newphi);
    s->pos.y = spheres[j].pos.y + r_traj*sin(theta)*cos(newphi);
    s->pos.z = spheres[j].pos.z + r_traj*sin(newphi);
    return collision;
}
 /* ===  FUNCTION  ======================================================================
 *         Name:  double_sphere_roll
 *  Description:  Roll sphere over two others until a third collision or until it reaches
 *  critical point and returns to free-fall mode. This step is very heavy on matrix oper-
 *  ations and is without a doubt the most expensive simulation mode.
 * =====================================================================================
 */
int double_sphere_roll(sphere *s,int j,int k, int i){
    // All notation follows the paper j=s2, k=s3 etc
    
    //We define a new coordinate basis u, v, w
    if(i==13){
	std::cout << "Culprit found" << std::endl;
    }
    //Vector between two spheres of contact
    vec3 c2c3 = spheres[k].pos - spheres[j].pos;   

    // Basis vector u: Normalized vector c2c3
    vec3 u(c2c3.x/c2c3.magnitude(), c2c3.y/c2c3.magnitude(), c2c3.z/c2c3.magnitude());
    u.equals(u.normalize());
    // Basis vector v: k-(k.u)u
    vec3 v( 0-(u.z)*u.x,
	    0-(u.z)*u.y,
	    1-(u.z)*u.z);
    v.equals(v.normalize());
    //Basis vector w: u x v
    vec3 w = cross_product(u,v);
    std::cout << "u is " << u.x << ", "<< u.y << " " << u.z<< std::endl;
    std::cout << "v is " << v.x << ", "<< v.y << " " << v.z<< std::endl;
    std::cout << "w is " << w.x << ", "<< w.y << " " << w.z<< std::endl;

    vec3 c2_omega = c2c3.scalar_multiply(
		dot_product(s->pos-spheres[j].pos,spheres[k].pos-spheres[j].pos)
		/pow(c2c3.magnitude(),2));
    vec3 c2c3_norm(c2c3.x/c2c3.magnitude(), c2c3.y/c2c3.magnitude(), c2c3.z/c2c3.magnitude());
    vec3 omega = c2c3_norm.scalar_multiply(dot_product(s->pos,c2c3_norm)); 
    omega = c2_omega + spheres[j].pos;
    vec3 omega_s = s->pos - omega;


    std::vector<int> intersects = check_intersects(omega_s.magnitude() + spheres[j].radius,
    		spheres[j].pos, i ,1,1,1);

    double beta_min = acos(w.x);
    double beta = acos(dot_product(omega_s,v)/omega_s.magnitude());
    if(dot_product(omega_s,w)<0) beta *= -1.;
    int sign_coef = 1;
    if(beta<0)sign_coef=-1;
    std::cout << "sign_coef is " << sign_coef << std::endl;
    std::cout << "beta is " << beta << std::endl;
    double T=-1;
    int collision_index =-1;
    int false_values=0;

    vec3 stform(dot_product(omega_s, u),dot_product(omega_s, v),dot_product(omega_s, w));
    std::cout << "s in our new coordinates is " << stform.x << ", "<< stform.y << " " << stform.z<< std::endl;
    /*-----------------------------------------------------------------------------
     *  Yet another attempt to define the horizon
     *-----------------------------------------------------------------------------*/
//    vec3 apex(omega.x, omega.y, omega.z+s->radius);
//    vec3 apex_pos = s->pos - apex;
//    int f = apex.z/apex_pos.z;
//    vec3 horizontal_dir = apex + apex_pos.scalar_multiply(f);
//    vec3 horizontal_pos = (final_dir.normalize()).scalar_multiply(s->radius);
//    vec3 horizonatal_pos_prime( 
//			u.x * horizontal_pos.x + 
//			u.y * horizontal_pos.y +
//			u.z * horizontal_pos.z, 
//			
//			v.x * horizontal_pos.x + 
//			v.y * horizontal_pos.y +
//			v.z * horizontal_pos.z, 
//
//			w.x * horizontal_pos.x + 
//			w.y * horizontal_pos.y +
//			w.z * horizontal_pos.z);
//    double beta_horizontal = acos(horizontal_pos_prime.z);
       
    for(int l=0; l<intersects.size(); l++)
    {
	//if(i==6) continue;
	if(intersects[l]==j || intersects[l]==k)
	{
	    false_values++;
	    continue;
	}
	//Set up the variables exactly as detailed in the paper.
	vec3 omega_c = spheres[intersects[l]].pos - omega;
	vec3 cprime(dot_product(omega_c, u),dot_product(omega_c, v),dot_product(omega_c, w));
	
	if(s->radius - pow(pow(pow( pow(cprime.y,2)+pow(cprime.z ,2) ,0.5) 
			    - s->radius - omega_s.magnitude(),2) + pow(cprime.x,2),0.5) 
			    + spheres[intersects[l]].radius <= SMIDGE){
	    false_values++;
	    continue;
		}
	double K1 	= cprime.y;
	double K2 	= cprime.z*sign_coef;
	double K3 	= (pow(cprime.x,2) + pow(cprime.y,2) + pow(cprime.z,2) 
		+ pow(omega_s.magnitude(),2) 
		- pow(s->radius+spheres[intersects[l]].radius,2))/(2*omega_s.magnitude());
    
	double D = pow(2*K1*K3,2)-4*(pow(K1,2)+pow(K2,2))*(pow(K3,2)-pow(K2,2));
	if(D<0)continue;
	double T1 = (2*K1*K3 + pow(D,0.5))/(2*(pow(K1,2)+pow(K2,2)));
	double T2 = (2*K1*K3 - pow(D,0.5))/(2*(pow(K1,2)+pow(K2,2)));
	std::cout   << "Found solutions with sphere: " << intersects[l] <<
			" with T1 corresponding to "<< T1 <<" and T2 to " <<
		       T2 << std::endl;
	std::cout << "Checking sphere " << intersects[l] << std::endl;	
    	double betaprime = sign_coef * acos(T1);
	vec3 sprime(0, omega_s.magnitude() * cos(betaprime), omega_s.magnitude()*sin(betaprime));
	vec3 omega_s_prime( u.x * sprime.x + 
			    v.x * sprime.y +
			    w.x * sprime.z, 
			    
			    u.y * sprime.x + 
			    v.y * sprime.y +
			    w.y * sprime.z, 

			    u.z * sprime.x + 
			    v.z * sprime.y +
			    w.z * sprime.z);
	vec3 s_prime=omega+omega_s_prime;

	std::cout <<"T1: "<< T1 <<" | "<< s_prime.x << " "<< s_prime.y << " " << s_prime.z<< std::endl;
    	betaprime = sign_coef * acos(T2);
	vec3 sprime2(0, omega_s.magnitude() * cos(betaprime), omega_s.magnitude()*sin(betaprime));
	vec3 omega_s_prime2( u.x * sprime2.x + 
			    v.x * sprime2.y +
			    w.x * sprime2.z, 
			    
			    u.y * sprime2.x + 
			    v.y * sprime2.y +
			    w.y * sprime2.z, 

			    u.z * sprime2.x + 
			    v.z * sprime2.y +
			    w.z * sprime2.z);
	vec3 s_prime2=omega+omega_s_prime2;
	std::cout <<"T2: "<< T2 <<" | "<< s_prime2.x << " "<< s_prime2.y << " " << s_prime2.z<< std::endl;
	//double sol = std::max(T1,T2);
	if(K2*(K3-K1*T1)<0)T1=T2;
	if(K2*(K3-K1*T2)<0)T2=T1;
	if(K2*(K3-K1*T2)<0)continue;
	double sol;
	
	/*-----------------------------------------------------------------------------
	 *  Old method that should be working
	 *-----------------------------------------------------------------------------*/
	/*if(abs(T1)<abs(T2)){
	    sol=T2;
	}	    
	else{
	    sol=T1;
	}	    
	if(abs(sol)<abs(T)){
	    collision_index=intersects[l];
	    T = sol;
	}*/	
	/*-----------------------------------------------------------------------------
	 *  Method described in paper
	 *-----------------------------------------------------------------------------*/
	sol = std::max(T1,T2);
	if(sol>T){
	    T=sol;
	    collision_index=intersects[l];
	}

    }
    if(T==1)
    {
	vec3 omega_s_prime = v.scalar_multiply(-omega_s.magnitude());
	s->pos.equals(omega+omega_s_prime);
	if(s->pos.z<=std::min(spheres[k].pos.z,spheres[j].pos.z))
	{
	    collision_index = -1;
	}
	else if(s->pos.z>spheres[k].pos.z){
	    collision_index = k;
	}
    }
    else
    {
	std::cout   << "We will use T=" << T <<
			" corresponding to "<< acos(T) << std::endl;
    	double betaprime = sign_coef * acos(T);
	//if(i==6)betaprime= betaprime;
	vec3 sprime(0, omega_s.magnitude() * cos(betaprime), omega_s.magnitude()*sin(betaprime));
	vec3 omega_s_prime( u.x * sprime.x + 
			    v.x * sprime.y +
			    w.x * sprime.z, 
			    
			    u.y * sprime.x + 
			    v.y * sprime.y +
			    w.y * sprime.z, 

			    u.z * sprime.x + 
			    v.z * sprime.y +
			    w.z * sprime.z);
       s->pos.equals(omega+omega_s_prime);	
    }

    if(s->pos.z-s->radius<0)
    {

    }
    return collision_index;
}

int main()
{
	//define random seed
	int time2 = time(NULL);
	srand(12345);
	//double ceil = Z_MAX/100.0;

	std::default_random_engine generator;
	std::exponential_distribution<double> distribution(35);
	//total volume occupied by spheres
	double volume=0;
	
	//Loop to attempt placement of NUM_SPHERES spheres
	for(int i=0; i < NUM_SPHERES; i++)
	{
		if(i==15000){

		}
		bool placed=false;

		//In this simulation "tries" keeps a count of all the attempts to place a 
		//sphere without immediate intersection(into a valid location) at the top
		//of a container. If this fails many times, give up.
		int tries = 0;
		while(!placed && tries<100000)
		{
			tries++;
			bool posfail = false;
			
			//Pick random radius and x,y
			double radius = distribution(generator);
			radius = rand_range(R_MIN, R_MAX);
			vec3 pos;
			pos.x = rand_range(radius, X_MAX - radius);
			pos.y = rand_range(radius, Y_MAX - radius);
			//pos.x = 5;
			//pos.y = 5;
			pos.z = Z_MAX - radius;
			if(pos.z-radius<0)posfail=true;
			//make sure there is no immediate intersection
			if(check_intersect(radius,pos,i,1,1,1)!=-1)
			{
			    posfail=true;
			}
			if(!posfail)
			{
			    sphere s;
			    s.radius = radius;
			    s.pos.equals(pos);

				/*if(i==0){
					s.radius = 2;
					s.pos.x=5;
					s.pos.y=5;
					s.pos.z=Z_MAX-s.radius;
				}
				if(i==1){
					s.radius = 2;
					s.pos.x=6;
					s.pos.y=6;
					s.pos.z=Z_MAX-s.radius;
				}
				if(i==2){
					s.radius = 1;
					s.pos.x=4;
					s.pos.y=4;
					s.pos.z=Z_MAX-s.radius;
				}
				if(i==3){
					s.radius = 1;
					s.pos.x=7;
					s.pos.y=8;
					s.pos.z=Z_MAX-s.radius;
				}
				if(i==4){
					s.radius = 1;
					s.pos.x=5.2;
					s.pos.y=5.3;
					s.pos.z=Z_MAX-s.radius;
				}
				if(i==5){
					s.radius = 1;
					s.pos.x=5.2;
					s.pos.y=5.3;
					s.pos.z=Z_MAX-s.radius;
				}
				if(i==6){
					s.radius = 1;
					s.pos.x=6;
					s.pos.y=6;
					s.pos.z=Z_MAX-s.radius;
				}
				  if(i==0){
					s.radius = 1;
					s.pos.x=1;
					s.pos.y=1;
					s.pos.z=Z_MAX-s.radius;
				}
				if(i==1){
					s.radius = 1;
					s.pos.x=1;
					s.pos.y=3;
					s.pos.z=Z_MAX-s.radius;
				}
				if(i==2){
					s.radius = 1;
					s.pos.x=1;
					s.pos.y=5;
					s.pos.z=Z_MAX-s.radius;
				}
				if(i==3){
					s.radius = 1;
					s.pos.x=1;
					s.pos.y=7;
					s.pos.z=Z_MAX-s.radius;
				}
				if(i==4){
					s.radius = 1;
					s.pos.x=3;
					s.pos.y=1;
					s.pos.z=Z_MAX-s.radius;
				}
				if(i==5){
					s.radius = 1;
					s.pos.x=3;
					s.pos.y=3;
					s.pos.z=Z_MAX-s.radius;
				}
				if(i==6){
					s.radius = 1;
					s.pos.x=3;
					s.pos.y=5;
					s.pos.z=Z_MAX-s.radius;
				}
				if(i==7){
					s.radius = 1;
					s.pos.x=3;
					s.pos.y=7;
					s.pos.z=Z_MAX-s.radius;
				}
				if(i==8){
					s.radius = 1;
					s.pos.x=5;
					s.pos.y=1;
					s.pos.z=Z_MAX-s.radius;
				}
				if(i==9){
					s.radius = 1;
					s.pos.x=5;
					s.pos.y=3;
					s.pos.z=Z_MAX-s.radius;
				}
				if(i==10){
					s.radius = 1;
					s.pos.x=5;
					s.pos.y=5;
					s.pos.z=Z_MAX-s.radius;
				}
				if(i==11){
					s.radius = 1;
					s.pos.x=5;
					s.pos.y=7;
					s.pos.z=Z_MAX-s.radius;
				}
				if(i==12){
					s.radius = 1;
					s.pos.x=3.1;
					s.pos.y=3.1;
					s.pos.z=Z_MAX-s.radius;
				}
				if(i==13){
					s.radius = 1;
					s.pos.x=2.9;
					s.pos.y=3.1;
					s.pos.z=Z_MAX-s.radius;
				}
				if(i==14){
					s.radius = 1;
					s.pos.x=2.9;
					s.pos.y=2.9;
					s.pos.z=Z_MAX-s.radius;
				}
				if(i==15){
					s.radius = 1;
					s.pos.x=3.1;
					s.pos.y=2.9;
					s.pos.z=Z_MAX-s.radius;
				}
				if(i==16){
					s.radius = 0.51;
					s.pos.x=3.1;
					s.pos.y=3;
					s.pos.z=Z_MAX-s.radius;
				}
				if(i==17){
					s.radius = 0.51;
					s.pos.x=3;
					s.pos.y=3.1;
					s.pos.z=Z_MAX-s.radius;
				}
				if(i==18){
					s.radius = 0.51;
					s.pos.x=2.9;
					s.pos.y=3;
					s.pos.z=Z_MAX-s.radius;
				}
				if(i==19){
					s.radius = 0.51;
					s.pos.x=3;
					s.pos.y=2.9;
					s.pos.z=Z_MAX-s.radius;
				}
				if(i==20){
					s.radius = 0.51;
					s.pos.x=3.1;
					s.pos.y=3;
					s.pos.z=Z_MAX-s.radius;
				}
				if(i==21){
					s.radius = 0.51;
					s.pos.x=3;
					s.pos.y=3.1;
					s.pos.z=Z_MAX-s.radius;
				}
				if(i==22){
					s.radius = 0.51;
					s.pos.x=2.9;
					s.pos.y=3;
					s.pos.z=Z_MAX-s.radius;
				}
				if(i==23){
					s.radius = 0.51;
					s.pos.x=3;
					s.pos.y=2.9;
					s.pos.z=Z_MAX-s.radius;
				}*/
				
			    bool lodged = false;
			    int state = 0;
			    int contact[3]={-1,-1,-1};
			    int collision;
			    //Logic to iterate between modes by interpreting the
			    //return codes
			    std::cout   << "NEW SPHERE: "
			    		<< i << std::endl;
			    std::cout   << "INIT:"
			    		<< s.pos.x << " "
						<< s.pos.y << " "
						<< s.pos.z<< std::endl;
			    while(!lodged)
			    {
				if(check_intersect(s.radius, s.pos, i, 1, 1, 1 )!=-1){
				    std::cout << "Something Bad Just Happened" << std::endl;
				}
				std::cout << "The contact array is currently: " << 
				    contact[0] << ", " <<
				    contact[1] << ", " <<
				    contact[2] << std::endl;
				switch(state)
				{

				    case 0:
				    collision=free_fall(&s, i);
				    std::cout   << "FF:"
				    		<< s.pos.x << " "
							<< s.pos.y << " "
							<< s.pos.z<< std::endl;
					if( s.pos.x-radius+SMIDGE<0  || 
					    s.pos.y-radius+SMIDGE<0  ||
					    s.pos.z-radius+SMIDGE<0  ||
					    s.pos.x+radius-SMIDGE>10 || 
					    s.pos.y+radius-SMIDGE>10 || 
					    s.pos.z+radius-SMIDGE>10){
						state=3;
						break;
					}
					if(collision==-1){
						state = 3;
					}
					else{
						state++;
						contact[0]=collision;
					}
					break;
				    case 1:
					collision=single_sphere_roll(&s, contact[0], i);
					std::cout   << "SSR:"
			    		<< s.pos.x << " "
						<< s.pos.y << " "
						<< s.pos.z<< std::endl;
					if( s.pos.x-radius+SMIDGE<0  || 
					    s.pos.y-radius+SMIDGE<0  ||
					    s.pos.z-radius+SMIDGE<0  ||
					    s.pos.x+radius-SMIDGE>10 || 
					    s.pos.y+radius-SMIDGE>10 || 
					    s.pos.z+radius-SMIDGE>10){	
							state=3;
						break;
					}
					if(collision!=-1){
					    state=2;
					    contact[1]=collision;
					}
					else{
						state = 0;
					    contact[0]=-1;
					}
					break;
				    case 2:
					collision=double_sphere_roll(&s, contact[0],contact[1], i);
					std::cout   << "DSR:"
				    		<< s.pos.x << " "
							<< s.pos.y << " "
							<< s.pos.z<< std::endl;
					if( s.pos.x-radius+SMIDGE<0  || 
					    s.pos.y-radius+SMIDGE<0  ||
					    s.pos.x+radius-SMIDGE>10 || 
					    s.pos.y+radius-SMIDGE>10){	
						state=3;
						break;
					}
					if(collision!=-1){
					    if(collision==contact[1]){
						state = 1;
						contact[0]=contact[1];
						contact[1]=-1;
					    }
					    else{
						state=3;
						contact[2]=collision;
					    }
					}
					else{
					    state=0;
					    contact[0]=-1;
					    contact[1]=-1;
					}
					break;
				    case 3:
				    	lodged=true;
					break;
				}
			    }
				if( s.pos.x-radius+SMIDGE<0  || 
				    s.pos.y-radius+SMIDGE<0  ||
				    s.pos.z-radius+SMIDGE<0  ||
				    s.pos.x+radius-SMIDGE>10 || 
				    s.pos.y+radius-SMIDGE>10 || 
				    s.pos.z+radius-SMIDGE>10){	
				state=3;
				break;
				}
				else{
					placed = true;
					spheres[i].radius = s.radius;
					spheres[i].pos.x = s.pos.x;
					spheres[i].pos.y = s.pos.y;
					spheres[i].pos.z = s.pos.z;
				}
			}
		}
		if(!placed)
		{
			//std::cout << "ERR::PLACEMENT // Placing sphere" << i << " failed" << std::endl;
		}
		else{
		volume += (4.0/3.0) * PI * pow(spheres[i].radius, 3.0);
		//if(volume/(ceil*X_MAX*Y_MAX)>=0.6)ceil+=Z_MAX/100.0;
		/*std::cout   << "INF::VOLUME // Occupied Volume Ratio is "
	    << volume/(X_MAX*Y_MAX*Z_MAX) << std::endl;//*/
		/*std::cout << "INF::PLACEMENT // Placed sphere " << i
		    << " of radius " << spheres[i].radius << " at coordinate ( "
		    << spheres[i].pos.x << " , " << spheres[i].pos.y << " , "
		    << spheres[i].pos.z << " )" << std::endl;*/
		std::cout   <<  spheres[i].radius << " "
			    << spheres[i].pos.x
			    <<  " " <<spheres[i].pos.y
			    <<  " " <<spheres[i].pos.z << std::endl << std::endl;/*
		std::cout   << "INF::VOLUME // Total Container Volume is "
			    << X_MAX*Y_MAX*Z_MAX << std::endl;
		std::cout   << "INF::VOLUME // Total Sphere Volume is "
	    << volume << std::endl;
		std::cout   << "INF::VOLUME // Occupied Volume Ratio is "
	    << volume/(X_MAX*Y_MAX*Z_MAX) << std::endl;//*/
		}
	}
	std::cout << "===================" << std::endl;
	for(int i=0; i<NUM_SPHERES; i++){
		std::cout   <<  spheres[i].radius << " "
			    << spheres[i].pos.x
			    <<  " " <<spheres[i].pos.y
			    <<  " " <<spheres[i].pos.z << std::endl;
	}
	/*std::cout   << "INF::VOLUME // Total Container Volume is "
		    << X_MAX*Y_MAX*Z_MAX << std::endl;
	std::cout   << "INF::VOLUME // Total Sphere Volume is "
		    << volume << std::endl;
	std::cout   << "INF::VOLUME // Occupied Volume Ratio is "
		    << volume/(X_MAX*Y_MAX*Z_MAX) << std::endl;*/
	return 0;
}

