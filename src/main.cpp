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
#include <string>
#include <fstream>
#include <sstream>
#include <ctime>
#include <iostream>
#include "sphere_math.hpp"
#include "sphere.hpp"
#include "distribution_selector.hpp"
#include <random>


//parameters
const int NUM_SPHERES = 10000000;
int X_MAX = 10;
int Y_MAX = 10;
int Z_MAX = 10;
double R_MIN = 0.1;
double R_MAX = 0.5;
double STEP = 0.001;


//Array of spheres to be filled
struct sphere spheres[NUM_SPHERES];

int check_intersect(double radius, vec3 pos, double i, bool fastx, bool fasty, bool fastz){
    for(int j=i-1; j >=0; j--)
    {
    	if(fastx && fabs(pos.x-spheres[j].pos.x) > radius + spheres[j].radius)continue;
    	if(fasty && fabs(pos.y-spheres[j].pos.y) > radius + spheres[j].radius)continue;
    	if(fastz && fabs(pos.z-spheres[j].pos.z) > radius + spheres[j].radius)continue;
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
    	if(fastx && fabs(pos.x-spheres[j].pos.x) > radius + spheres[j].radius)continue;
    	if(fasty && fabs(pos.y-spheres[j].pos.y) > radius + spheres[j].radius)continue;
    	if(fastz && fabs(pos.z-spheres[j].pos.z) > radius + spheres[j].radius)continue;
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
    	if(fabs(s->pos.x-spheres[j].pos.x) > s->radius + spheres[j].radius)continue;
    	if(fabs(s->pos.y-spheres[j].pos.y) > s->radius + spheres[j].radius)continue;
		if( distance(vec3(s->pos.x,s->pos.y,0), vec3(spheres[j].pos.x,spheres[j].pos.y,0))
			< spheres[j].radius + s->radius + SMIDGE)
		{
			double z_col=spheres[j].pos.z +pow(pow(spheres[j].radius +  s->radius , 2)
					- pow(spheres[j].pos.x -  s->pos.x , 2) -
					pow(spheres[j].pos.y -  s->pos.y , 2), 0.5);
			if(z_col > z_max && s->pos.z>z_col)
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
		std::cout << "Entering SSR" <<std::endl;
		//Limit the amount of spheres we run the math against to the radius of the two
	//touching spheres.
    //Assemble the spherical coordinates that we will need to run the toroidal collision math
    vec3 v_traj 	=(s->pos-spheres[j].pos);
    double r_traj 	=v_traj.magnitude();
    double d = sqrt(pow((s->pos.x-spheres[j].pos.x),2) + pow((s->pos.y-spheres[j].pos.y),2));
    double cosphi     = d/(r_traj);
    double sinphi     = (s->pos.z-spheres[j].pos.z)/(r_traj);
    double costheta   = (s->pos.x-spheres[j].pos.x)/d;
    double sintheta   = (s->pos.y-spheres[j].pos.y)/d;

    std::cout << cosphi << std::endl;
    std::cout << sinphi << std::endl;
    std::cout << sintheta << std::endl;
    std::cout << costheta << std::endl;
    /*-----------------------------------------------------------------------------
     *  Add checks for redundancies of the angles.
     *-----------------------------------------------------------------------------*/
    double phi 		= asin(sinphi);
    double theta 	= asin(sintheta);
    vec3 u( costheta,	sintheta, 0);
    vec3 v(-sintheta, costheta, 0);
    vec3 w(0,0,1);
   
    double r_collision_volume = (s->pos-spheres[j].pos).magnitude()+s->radius;
    int collision = -1;
    
    std::vector<int> intersects = check_intersects( r_collision_volume,
    		                                    spheres[j].pos, i ,1,1,1);

    double maxsol = 1;
    double minz = 0;
    std::cout<<"maxsol is " << maxsol <<std::endl;
    /*-----------------------------------------------------------------------------
     *  New phi and T can also be pi
     *-----------------------------------------------------------------------------*/
    double newphi=0;
    double T = cos(0);
    for(int l=0; l<intersects.size(); l++){
        std::cout << "intersect: " << intersects[l] <<std::endl;
        /*-----------------------------------------------------------------------------
         *  I only need to do this oncespheres[j].pos.z + r_traj*sin(newphi);

         *-----------------------------------------------------------------------------*/
        vec3 sj = spheres[intersects[l]].pos - spheres[j].pos;
        if(dot_product((vec3){v_traj.x,v_traj.y,0},(vec3){sj.x,sj.y,0})<0)continue;
        
        /*-----------------------------------------------------------------------------
         * This is most likely wrong. u.x * sj.x + u.y*sj.y...
         *-----------------------------------------------------------------------------*/
        vec3 sjprime(       u.x * sj.x + 
			    u.y * sj.y +
			    u.z * sj.z, 
			    
			    v.x * sj.x + 
			    v.y * sj.y +
			    v.z * sj.z, 

			    w.x * sj.x + 
			    w.y * sj.y +
			    w.z * sj.z);
        if(     spheres[intersects[l]].radius+
                s->radius-
                sqrt(pow(sqrt(pow(sjprime.x,2)+
                pow(sjprime.z,2))-r_traj,2)+
                pow(sjprime.y,2))
                <=0
        )continue;
        
        if(intersects[l]==j || intersects[l]==i)continue;
    	double Wx 	= (spheres[intersects[l]].pos.x-spheres[j].pos.x)/(s->radius+spheres[j].radius);
    	double Wy 	= (spheres[intersects[l]].pos.y-spheres[j].pos.y)/(s->radius+spheres[j].radius);
    	double Wz 	= (spheres[intersects[l]].pos.z-spheres[j].pos.z)/(s->radius+spheres[j].radius);
    	double W 	= (spheres[intersects[l]].radius+s->radius)/(s->radius+spheres[j].radius);

    	double K1 	= 2.0*(Wx*costheta+Wy*sintheta);
    	double K2 	= 2.0*Wz;
    	double K3 	= pow(Wx,2) + pow(Wy,2) + pow(Wz,2) +1 - pow(W,2);

        std::cout<<"     Wx:" << Wx << std::endl;
        std::cout<<"     Wy:" << Wy << std::endl;
        std::cout<<"     Wz:" << Wz << std::endl;
        std::cout<<"     W:" << W << std::endl;
        std::cout<<"     K1:" << K1 << std::endl;
        std::cout<<"     K2:" << K2 << std::endl;
        std::cout<<"     K3:" << K3 << std::endl;

        //Now solve the quadratic equation
    	double D = pow(2*K1*K3,2)-4*(pow(K1,2)+pow(K2,2))*(pow(K3,2)-pow(K2,2));
    	if(fabs(D)<0.0000000001)D=0;
        std::cout<<"     D is " << D <<std::endl;
        if(D<0) continue;
    	double T1 = (2*K1*K3 + sqrt(D))/(2*(pow(K1,2)+pow(K2,2)));
    	double T2 = (2*K1*K3 - sqrt(D))/(2*(pow(K1,2)+pow(K2,2)));
        std::cout<<"     T1 is " << T1 <<std::endl;
        std::cout<<"     T2 is " << T2 <<std::endl;
        std::cout << "     K2*(K3-K1*T1)="<<K2*(K3-K1*T1) <<std::endl; 
        std::cout << "     K2*(K3-K1*T2)="<<K2*(K3-K1*T2) <<std::endl; 
        if(K2*(K3-K1*T1)<0)
        {
            std::cout << "     changed" <<std::endl;
            T1=T2;
        }
    	if(K2*(K3-K1*T2)<0)
        {
            std::cout << "     changed" <<std::endl;
            T2=T1;
        }
    	if(K2*(K3-K1*T2)<0)continue;
        double sol=0;
        if(spheres[intersects[l]].pos.z-spheres[j].pos.z>=0)
        {
            sol=std::min(T1,T2);
        }
        else
        {
            sol=std::max(T1,T2);
            double phi_l = -asin((spheres[j].pos.z-spheres[intersects[l]].pos.z)/((spheres[j].pos-spheres[intersects[l]].pos).magnitude()));
//            if(spheres[j].pos.z>spheres[intersects[l]].pos.z+spheres[intersects[l]].radius){
 
            double relx = r_traj*costheta*sol;
            double rely = r_traj*sintheta*sol;
            if((vec3){relx,rely,0}.magnitude()<((vec3){v_traj.x,v_traj.y,0}).magnitude()){
            std::cout<<"     skipped due to both intersections being below z=0"<<std::endl;
                std::cout<<"     -acos(sol):" << -acos(sol)<<std::endl;
                std::cout<<"     phi_l: " <<phi_l<<std::endl;
                continue;
            }
        }
        std::cout<<"     sol is " << sol <<std::endl;
	double z =  r_traj*sin(acos(sol));
        std::cout<<"     z is " << z <<std::endl;
        /*-----------------------------------------------------------------------------
         *  We dont understand this
         *-----------------------------------------------------------------------------*/
	if(z>minz && z <(s->pos-spheres[j].pos).z ){
	    maxsol = sol;
	    minz =  z;
            std::cout<<"maxsol is " << maxsol <<std::endl;
	    collision = intersects[l];
	}
    }
    if(collision==-1)maxsol=1;
    newphi = acos(maxsol);
    
    std::cout <<"phi is:" <<newphi <<std::endl;
    std::cout << "hit against: " << collision <<std::endl;
    std::cout << "new z is:" << spheres[j].pos.z + r_traj*sin(newphi) << std::endl;

    /*  if(spheres[j].pos.z + r_traj*sin(newphi) - s->radius < 0 )
    {   
        collision==-1;
        newphi=asin((s->pos.z-spheres[j].pos.z)/(r_traj));
    }*/
    
    s->pos.x = spheres[j].pos.x + r_traj*costheta*maxsol;
    s->pos.y = spheres[j].pos.y + r_traj*sintheta*maxsol;
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
    std::cout << "Entering DSR" <<std::endl;

    int collision_index = -1;
    //Vector between two spheres of contact
    vec3 c2c3 = spheres[k].pos - spheres[j].pos; 

    //We define a new coordinate basis u, v, w
    // Basis vector u: Normalized vector c2c3
    vec3 u = c2c3.normalize();
    // Basis vector v: k-(k.u)u
    vec3 v( 0-(u.z)*u.x,
	    0-(u.z)*u.y,
	    1-(u.z)*u.z);
    v.equals(v.normalize());
    //Basis vector w: u x v
    vec3 w = cross_product(u,v);
    std::cout << "u: " << u.x << ", "<< u.y << " " << u.z<< std::endl;
    std::cout << "v: " << v.x << ", "<< v.y << " " << v.z<< std::endl;
    std::cout << "w: " << w.x << ", "<< w.y << " " << w.z<< std::endl;

    vec3 c2_omega = c2c3.scalar_multiply(
		                        dot_product(s->pos-spheres[j].pos,
                                                    spheres[k].pos-spheres[j].pos)
		                        /pow(c2c3.magnitude(),2));
    vec3 c2c3_norm  = c2c3.normalize();
    vec3 omega      = c2_omega + spheres[j].pos;
    vec3 omega_s    = s->pos - omega;

    // Make a list of spherers intersecting the collision volume centered at omega 
    std::vector<int> E4 = check_intersects(omega_s.magnitude() + spheres[j].radius,
    		omega , i ,1,1,1);

    std::vector<int> intersects;
    for(int l=0; l<E4.size(); l++)
    {
	std::cout   << "E4: " << E4[l] << std::endl;
        vec3 omega_c = spheres[E4[l]].pos - omega;
        vec3 cprime(dot_product(omega_c, u),dot_product(omega_c, v),dot_product(omega_c, w));
        std::cout <<"For " << E4[l] <<"our criterion is"  <<
            s->radius - pow(pow(pow(pow(cprime.y,2)+pow(cprime.z ,2) ,0.5) 
			    - s->radius - omega_s.magnitude(),2) + pow(cprime.x,2),0.5) 
			    + spheres[E4[l]].radius << std::endl;
        if(s->radius - pow(pow(pow(pow(cprime.y,2)+pow(cprime.z ,2) ,0.5) 
			    - s->radius - omega_s.magnitude(),2) + pow(cprime.x,2),0.5) 
			    + spheres[E4[l]].radius >= -0.2)
        {
            intersects.push_back(E4[l]);
        }
    }

    if(intersects.size()==0)
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
        return collision_index;
    }
    //Set up Beta
    double beta_min = acos(w.x);
    double beta = acos(dot_product(omega_s,v)/omega_s.magnitude());
    if(dot_product(omega_s,w)<0) beta *= -1.;
    int sign_coef = 1;
    if(beta<0)sign_coef=-1;
    double minang=beta;
    double maxang=PI;

    if(beta<0)
    {
        minang=PI;
        maxang=beta;
    }

    std::cout << "sign_coef is " << sign_coef << std::endl;
    std::cout << "beta is " << beta << std::endl;
    int false_values=0;

    double T=-1;
    for(int l=0; l<intersects.size(); l++)
    {
	std::cout   << "intersect: " << intersects[l] << std::endl;
	//if(i==6) continue;
	if(intersects[l]==j || intersects[l]==k)
	{
	    false_values++;
	    continue;
	}
	//Set up the variables exactly as detailed in the paper.
	vec3 omega_c = spheres[intersects[l]].pos - omega;
	vec3 cprime(dot_product(omega_c, u),dot_product(omega_c, v),dot_product(omega_c, w));
	
	double K1 	= cprime.y;
	double K2 	= cprime.z*sign_coef;
	double K3 	= (pow(cprime.x,2) + pow(cprime.y,2) + pow(cprime.z,2) 
		+ pow(omega_s.magnitude(),2) 
		- pow(s->radius+spheres[intersects[l]].radius,2))/(2*omega_s.magnitude());
    
	double D = pow(2*K1*K3,2)-4*(pow(K1,2)+pow(K2,2))*(pow(K3,2)-pow(K2,2));
    	if(fabs(D)<0.0000000001)D=0;
        std::cout<<"     D is " << D <<std::endl;
	if(D<0)continue;
	double T1 = (2*K1*K3 + pow(D,0.5))/(2*(pow(K1,2)+pow(K2,2)));
	double T2 = (2*K1*K3 - pow(D,0.5))/(2*(pow(K1,2)+pow(K2,2)));
	std::cout   << "     T1: "<< T1 << std::endl;
	std::cout   << "     T2: "<< T2 << std::endl;
	std::cout << "     Checking sphere " << intersects[l] << std::endl;	
    	double betaprime = sign_coef * acos(T1);
	vec3 sprime(0, omega_s.magnitude() * cos(betaprime), omega_s.magnitude()*sin(betaprime));
	vec3 omega_s_prime( u.x * sprime.x + 
			    u.y * sprime.y +
			    u.z * sprime.z, 
			    
			    v.x * sprime.x + 
			    v.y * sprime.y +
			    v.z * sprime.z, 

			    w.x * sprime.x + 
			    w.y * sprime.y +
			    w.z * sprime.z);
	vec3 s_prime=omega+omega_s_prime;

	std::cout <<"     T1: "<< T1 <<" | "<< s_prime.x << " "<< s_prime.y << " " << s_prime.z<< std::endl;
    	betaprime = sign_coef * acos(T2);
	vec3 sprime2(0, omega_s.magnitude() * cos(betaprime), omega_s.magnitude()*sin(betaprime));
	vec3 omega_s_prime2(u.x * sprime2.x + 
			    u.y * sprime2.y +
			    u.z * sprime2.z, 
			    
			    v.x * sprime2.x + 
			    v.y * sprime2.y +
			    v.z * sprime2.z, 

			    w.x * sprime2.x + 
			    w.y * sprime2.y +
			    w.z * sprime2.z);
	vec3 s_prime2=omega+omega_s_prime2;
	std::cout <<"     T2: "<< T2 <<" | "<< s_prime2.x << " "<< s_prime2.y << " " << s_prime2.z<< std::endl;
	//double sol = std::max(T1,T2);
	if(K2*(K3-K1*T1)<-SMIDGE)T1=T2;
	if(K2*(K3-K1*T2)<-SMIDGE)T2=T1;
	if(K2*(K3-K1*T2)<-SMIDGE){
            std::cout << "     skipping: K2*(K3-K1*T2) = " << K2*(K3-K1*T2) << std::endl;
            continue;
        }
	double sol;
	/*-----------------------------------------------------------------------------
	 *  Method described in paper
	 *-----------------------------------------------------------------------------*/
	sol = std::max(T1,T2);
	if(sol>T){  
            double betaprime = sign_coef * acos(sol);
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
            if(s->pos.z>(omega+omega_s_prime).z)
            {
                T=sol;
	        collision_index=intersects[l];
            }
	}

    }
	std::cout   << "We will use T=" << T <<
			" corresponding to "<< sign_coef*acos(T) <<" and sphere:"
                        << collision_index << std::endl;
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

    if(s->pos.z-s->radius<0)
    {

    }
    return collision_index;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  sphere_cage_gen
 *  Description:  Places several massive spheres at the edges of the container 
 *  to simulate flat walls
 * =====================================================================================
 */
int sphere_cage_gen(int num_spheres, double r)
{
    spheres[num_spheres]=(sphere){r, (vec3){-r, Y_MAX/2.0, Z_MAX/2.0}}; 
    num_spheres++;
    spheres[num_spheres]=(sphere){r, (vec3){X_MAX+r, Y_MAX/2.0, Z_MAX/2.0}}; 
    num_spheres++;
    spheres[num_spheres]=(sphere){r, (vec3){X_MAX/2.0, -r, Z_MAX/2.0}}; 
    num_spheres++;
    spheres[num_spheres]=(sphere){r, (vec3){X_MAX/2.0, Y_MAX + r, Z_MAX/2.0}}; 
    num_spheres++;
    spheres[num_spheres]=(sphere){r, (vec3){X_MAX/2.0, Y_MAX/2.0, -r}}; 
    num_spheres++;
    spheres[num_spheres]=(sphere){r, (vec3){X_MAX/2.0, Y_MAX/2.0, Z_MAX + r}};
    num_spheres++;
    return num_spheres;
}
int main(int argc, char* argv[])
{
        DistributionSelector distribution;
	std::default_random_engine generator;
        for(int i=0; i< argc; i++)
        {
            std::string arg = argv[i];
            //Help
            if(argc==1)
            {
                break;
            }
            else if(arg=="-h")
            {
                std::cout<< "Usage: " << argv[0] 
                    << " NUM_SPHERES" << " X_MAX" << " Y_MAX" << " Z_MAX" 
                    << " R_MIN" << " R_MAX" << " STEP" <<  std::endl;
                std::cout<< "No args implies hard-coded defaults." << std::endl;
                return 0;
            }
            else if(argc>=6) {
                //std::cout << i << ": " <<arg <<std::endl;
                switch(i)
                {
                    case 1:
                        //NUM_SPHERES = std::stoi(arg);
                        break;
                    case 2:
                            X_MAX = std::stoi(arg);
                        break;
                    case 3:
                        Y_MAX = std::stoi(arg);
                        break;
                    case 4:
                        Z_MAX = std::stoi(arg);
                        break;
                    case 5:
                        //What distribution should I use?
                            if(arg=="normal")
                            {
                                distribution.set(DistributionSelector::normal,
                                        std::stod(argv[i+1]), std::stod(argv[i+2]));
                                i+=2;
                            }
                            else if(arg=="gamma")
                            {
                                distribution.set(DistributionSelector::gamma,
                                        std::stod(argv[i+1]), std::stod(argv[i+2]));
                                i+=2;
                            }
                            else if(arg=="uniform")
                            {
                                distribution.set(DistributionSelector::uniform,
                                        std::stod(argv[i+1]), std::stod(argv[i+2]));
                                i+=2;
                            }
                            else if(arg=="exponential")
                            {
                                distribution.set(DistributionSelector::exponential,
                                        std::stod(argv[i+1]));
                                i+=2;
                            }
                            else if(arg=="chi_squared")
                            {
                                distribution.set(DistributionSelector::chi_squared,
                                        std::stod(argv[i+1]));
                                i+=2;
                            }
                            break;
                }
            }
            else if(argc!=9){
                std::cout << "Invalid number of arguments. " << std::endl;
                return 0;
            }
        }
        time_t now = time(0);
        tm *tm = localtime(&now);
        std::ostringstream oss;
        oss << "sda-coords-" 
            << 1900+tm->tm_year << "-"<< 1+tm->tm_mon << "-"<< tm->tm_mday 
            << "-"<< 1+tm->tm_hour << 1+tm->tm_min<<1+tm->tm_sec
            << "-" << argv[5] << "-" << argv[6] << "-" 
            << X_MAX <<"_" << Y_MAX <<"_"<< Z_MAX; 
        std::string fn = oss.str();
        std::cout << "fn is " << fn <<std::endl;
        std::cout << "argc is " << argc <<std::endl;
        std::cout << "R_MIN is " << R_MIN <<std::endl;
        std::cout << "R_MAX is " << R_MAX <<std::endl;
        std::ofstream out(fn);

        int sphere_count=0;
        sphere_count = sphere_cage_gen(sphere_count, 100000);
        int time2 = time(NULL);
	srand(time2); //define random seed


	double volume=0; //total volume occupied by spheres
	
	//Loop to attempt placement of NUM_SPHERES spheres
        for(int i=sphere_count; i < NUM_SPHERES; i++)
	{
		bool placed=false;

		//In this simulation "tries" keeps a count of all the attempts to place a 
		//sphere without immediate intersection(into a valid location) at the top
		//of a container. If this fails many times, give up.
		int tries = 0;
		double radius = -1;
		while(radius <= 0){
		radius = distribution.gen();
		}
                while(!placed && tries<1000000)
		{
			tries++;
			bool posfail = false;
			
			//Pick random radius and x,y
			vec3 pos;
			pos.x = rand_range(radius, X_MAX - radius);
			pos.y = rand_range(radius, Y_MAX - radius);
			pos.z = rand_range(radius, Z_MAX - radius);
			//pos.z = Z_MAX - radius;
			//pos.x = 5;
			//pos.y = 5;
			//pos.z = Z_MAX - radius;
			//if(pos.z-radius<0)posfail=true;
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
			    bool lodged = false;
			    int state = 0;
			    int contact[3]={-1,-1,-1};
			    int collision;
			    //Logic to iterate between modes by interpreting the
			    //return codes
	                    std::cout << "========================================================" << std::endl;
			    std::cout   << "NEW SPHERE: "
			    		<< i << std::endl;
			    std::cout   << "INIT:"
			    		<< s.pos.x << " "
						<< s.pos.y << " "
						<< s.pos.z<< std::endl;
			    while(!lodged)
			    {
				if(check_intersect(s.radius-SMIDGE, s.pos, i, 1, 1, 1 )!=-1){
				    std::cout << "Something Bad Just Happened" << std::endl;
				    std::cout << "Intersection of magnitude: " << 
                                        (spheres[check_intersect(s.radius-SMIDGE, s.pos, i, 1, 1, 1)].pos - s.pos).magnitude()
                                        -0.5*2<< std::endl;
				    std::cout << "Intersection with: "<< check_intersect(s.radius-SMIDGE, s.pos, i, 1, 1, 1) << std::endl; 
				}
			        if(s.pos.x-s.radius<0 || s.pos.y-s.radius<0 || s.pos.z-s.radius <0 || s.pos.x+s.radius>X_MAX || s.pos.y+s.radius>Y_MAX || s.pos.z+s.radius > Z_MAX)break;
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
					placed = true;
					spheres[i].radius = s.radius;
					spheres[i].pos.x = s.pos.x;
					spheres[i].pos.y = s.pos.y;
					spheres[i].pos.z = s.pos.z;
			}
		}
		if(!placed)
		{
			//std::cout << "ERR::PLACEMENT // Placing sphere" << i << " failed" << std::endl;
		}
                if(tries>=1000000)break;
		else{
                sphere_count++;
                volume += (4.0/3.0) * PI * pow(spheres[i].radius, 3.0);
		std::cout   <<  spheres[i].radius << " "
			    << spheres[i].pos.x
			    <<  " " <<spheres[i].pos.y
			    <<  " " <<spheres[i].pos.z << std::endl << std::endl;
                out << i << " "
                            << spheres[i].radius << " "
			    << spheres[i].pos.x
			    <<  " " <<spheres[i].pos.y
			    <<  " " <<spheres[i].pos.z << std::endl;
		}
	}
        out.close();
	std::cout << "===================" << std::endl;
	for(int i=0; i<sphere_count; i++){
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

