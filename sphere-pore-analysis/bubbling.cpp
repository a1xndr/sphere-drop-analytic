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
 *   Organization:··
 *
 * =====================================================================================
 */

#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "sphere_math.hpp"
#include "sphere.hpp"
#include "bubble.hpp"

const int NUM_SPHERES = 2000000;
const int NUM_BUBBLES = 2000000;
const int MAX_TRIES = 100000;
const int X_MAX = 10;
const int Y_MAX = 10;
const int Z_MAX = 10;
const double STEP_SIZE = 0.01;
const double R_STEP_SIZE = 0.001;




struct sphere spheres[NUM_SPHERES];
struct bubble bubbles[NUM_BUBBLES];


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  check_intersect
 *  Description:  Checks a sphere of radius radius amd position pos against the array of
 *  spheres and bubbles for intersections
 * =====================================================================================
 */
int check_intersect(double radius, vec3 pos, int i, int j, bool fastx, bool fasty, bool fastz){
    for(int k=i-1; k >=0; k--)
    {
        if(fastx && abs(pos.x-spheres[k].pos.x) > radius + spheres[k].radius)continue;
        if(fasty && abs(pos.y-spheres[k].pos.y) > radius + spheres[k].radius)continue;
        if(fastz && abs(pos.z-spheres[k].pos.z) > radius + spheres[k].radius)continue;
        //Square of distance between two centers comparison against square of radius.
        if( distance(pos,spheres[k].pos)
            < spheres[k].radius + radius - SMIDGE)return k;
    }
    for(int k=j-1; k >=0; k--)
    {
        if(fastx && abs(pos.x-bubbles[k].pos.x) > radius + bubbles[k].radius)continue;
        if(fasty && abs(pos.y-bubbles[k].pos.y) > radius + bubbles[k].radius)continue;
        if(fastz && abs(pos.z-bubbles[k].pos.z) > radius + bubbles[k].radius)continue;
        //Square of distance between two centers comparison against square of radius.
        if( distance(pos,bubbles[k].pos)
            < bubbles[k].radius + radius - SMIDGE)return NUM_SPHERES+k;
    }
    return -1;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  find_nearest
 *  Description:  Finds the nearest sphere to a sphere with radius radius and at position 
 *  pos. Checks against both the sphere and bubble arrays up to inidces i and j respectively
 * =====================================================================================
 */
/*/int find_nearest(double radius, vec3 pos, int i, int j)
{
    double min;        // Distance to closest sphere
    int min_index;  // Index of closest sphere
    for(int k=i-1; k>=0; k--)
    {
        double dist = distance(pos,spheres[k].pos) - spheres[k].radius - radius;
        if( dist < min )
        {
            min = dist;
            min_index = k;
        }
    }
    for(int k=j-1; k>=0; k--)
    {
        double dist = distance(pos,bubbles[k].pos) - bubbles[k].radius - radius;
        if( dist < min )
        {
            min = dist;
            min_index = k;
        }
    }
    return min_index;
}
*/

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  nearest_valid_sphere
 *  Description:  Look for a sphere or bubble that is closest to the bubble of radius:
 *  radius, position:pos with the condition that the contact point forms an acute angle
 *  with the normal.
 * =====================================================================================
 */
/*  int nearest_valid_sphere(double radius, vec3 pos, std::vector<vec3>* constraints, int i, int j)
{

    double min;        // Distance to closest sphere
    int min_index;  // Index of closest sphere
    for(int k=i-1; k>=0; k--)
    {
        double dist = distance(pos,spheres[k].pos) - spheres[k].radius - radius;
        if( dist < min )
        {
            bool valid=true;
            for(int l=0; l<constraints.size(); l++)
            {
                if(spheres[k].dot_product(constraints[l])
            }
            min = dist;
            min_index = k;
        }
    }
    for(int k=j-1; k>=0; k--)
    {
        double dist = distance(pos,bubbles[k].pos) - bubbles[k].radius - radius;
        if( dist < min )
        {
            min = dist;
            min_index = k;
        }
    }
    return min_index;
}

*/
int read_sphere_coords(){
    std::ifstream file("spheres.list");
    int count = 0;
    if(file.is_open())
    {
        std::string line;
        while(getline(file, line))
        {
            std::istringstream iss(line);
            double d;
            double params[4];
            int paramcount=0;;
            while (iss>>d){
                params[paramcount]=d;
                paramcount++;
            }
            spheres[count]=(sphere){params[0], (vec3){params[1],params[2],params[3]}};
            count++;
        }
    }
    return count;
}

int main()
{
    int sphere_count    = read_sphere_coords(); 
    int bubble_count    = 0;
    while(1)
    {
        int tries=0;
        bubble b;
        while(tries<MAX_TRIES)
        {
            tries++;
            b.pos.x = rand_range(0, X_MAX);
            b.pos.y = rand_range(0, Y_MAX);
            b.pos.z = rand_range(0, Z_MAX);
            b.radius=0;
            if(check_intersect(0, b.pos, sphere_count, bubble_count, 1, 1, 1)!=-1) continue;
            bool placed = false;
            std::vector<vec3> normals;
            std::vector<int> normal_indices;
            while(!placed)
            {
                vec3 trajectory{0,0,0};
                for(int l=0; l<normals.size(); l++)
                {
                    trajectory = trajectory + normals[l];
                }
                if(trajectory.x != 0 || trajectory.y != 0 || trajectory.z != 0)
                {
                    trajectory = trajectory.normalize();
                }
                b.pos = b.pos + trajectory.scalar_multiply(STEP_SIZE);
                bool grown = false;
                while(!grown)
                {
                    b.radius += R_STEP_SIZE;
                    int collision = check_intersect(b.radius, b.pos, sphere_count, bubble_count, 1, 1, 1);
                    
                    if(b.pos.x-b.radius < 0)
                    {
                        b.radius -= R_STEP_SIZE;
                        bool skip=false;
                        for(int l=0; l<normal_indices.size(); l++)
                        {
                            if(normal_indices[l]==-1){
                                skip=true;
                            }
                        }
                        if(!skip){
                            normals.push_back((vec3){1,0,0});
                            normal_indices.push_back(-1);
                        }
                        grown = true;
                    }
                    if(b.pos.x + b.radius > X_MAX){
                        b.radius -= R_STEP_SIZE;
                        bool skip=false;
                        for(int l=0; l<normal_indices.size(); l++)
                        {
                            if(normal_indices[l]==-2){
                                skip=true;
                            }
                        }
                        if(!skip){
                            normals.push_back((vec3){-1,0,0});
                            normal_indices.push_back(-2);
                        }
                        grown = true;
                    }
                    if(b.pos.y-b.radius < 0){
                        b.radius -= R_STEP_SIZE;
                        bool skip=false;
                        for(int l=0; l<normal_indices.size(); l++)
                        {
                            if(normal_indices[l]==-3){
                                skip=true;
                            }
                        }
                        if(!skip){
                            normals.push_back((vec3){0,1,0});
                            normal_indices.push_back(-3);
                        }
                        grown = true;
                    }
                    if(b.pos.y + b.radius > Y_MAX){
                        b.radius -= R_STEP_SIZE;
                        bool skip=false;
                        for(int l=0; l<normal_indices.size(); l++)
                        {
                            if(normal_indices[l]==-4){
                                skip=true;
                            }
                        }
                        if(!skip){
                            normals.push_back((vec3){0,-1,0});
                            normal_indices.push_back(-4);
                        }
                        grown = true;
                    }
                    if(b.pos.z-b.radius < 0){
                        b.radius -= R_STEP_SIZE;
                        bool skip=false;
                        for(int l=0; l<normal_indices.size(); l++)
                        {
                            if(normal_indices[l]==-5){
                                skip=true;
                            }
                        }
                        if(!skip){
                            normals.push_back((vec3){0,0,1});
                            normal_indices.push_back(-5);
                        }
                        grown = true;
                    }
                    if(b.pos.z + b.radius > Z_MAX){
                        b.radius -= R_STEP_SIZE;
                        bool skip=false;
                        for(int l=0; l<normal_indices.size(); l++)
                        {
                            if(normal_indices[l]==-6){
                                skip=true;
                            }
                        }
                        if(!skip){
                            normals.push_back((vec3){0,0,-1});
                            normal_indices.push_back(-6);
                        }
                        grown = true;
                    }
                    if(collision!=-1)
                    {
                        int col_bubble = false;
                        if(collision>=NUM_SPHERES)
                        {
                            bool skip = false;
                            for(int l=0; l<normal_indices.size(); l++)
                            {
                              if(normal_indices[l]==collision){
                                    normals[l]=((vec3)(b.pos - bubbles[collision-NUM_SPHERES].pos).normalize());
                                    skip=true;
                                }
                            }
                            if(!skip)
                            {
                                normal_indices.push_back(collision);
                                collision -=NUM_SPHERES;
                                normals.push_back((vec3)(b.pos - bubbles[collision].pos).normalize());
                            }
                        }
                        else
                        {
                            bool skip = false;
                            for(int l=0; l<normal_indices.size(); l++)
                            {
                                if(normal_indices[l]==collision){
                                    normals[l]=((vec3)(b.pos - spheres[collision].pos).normalize());
                                    skip=true;
                                }
                            }
                            if(!skip)
                            {
                                normal_indices.push_back(collision);
                                normals.push_back((vec3)(b.pos - spheres[collision].pos).normalize());
                            }
                        }
                        b.radius -= R_STEP_SIZE;
                        grown = true;
                    }
                /*for(int l=0; l<normals.size(); l++)
                {
                    //std::cout << l << " " << normals[l].x << " " << normals[l].y << " " << normals[l].z << std::endl;
                }*/
                }
                if(normals.size()>=3){
                    placed = true;
                    for(int l=0; l<normal_indices.size(); l++)
                    {
                        int index = normal_indices[l]-NUM_SPHERES;
                        if(index>=0){
                            b.neighboors.push_back(index);
                            bubbles[index].neighboors.push_back(bubble_count);
                        }
                    }
                }
            }
            bubbles[bubble_count]=b;
            std::cout << b.radius << " " << b.pos.x << " " << b.pos.y << " " << b.pos.z << std::endl;
 //           std::cout << "Neighboors: ";
            for(int l=0; l<bubbles[bubble_count].neighboors.size(); l++)
            {
                std::cout << bubbles[bubble_count].neighboors[l] <<" ";
            }
            std::cout << std::endl;
            break;
        }
        bubble_count++;
    }
}



