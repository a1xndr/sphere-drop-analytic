/*
 * =====================================================================================
 *
 *       Filename:  bubble.cpp
 *
 *    Description:  Splits a bubble distribution into individual pores
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




struct bubble bubbles[NUM_BUBBLES];
int bubble_map[NUM_BUBBLES]={-1};

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  check_intersect
 *  Description:  Checks a sphere of radius radius amd position pos against the array of
 *  spheres and bubbles for intersections
 * =====================================================================================
 */
int check_intersect(double radius, vec3 pos, int j, bool fastx, bool fasty, bool fastz){
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

int read_bubble_coords(){
    std::ifstream file("bubbles.list");
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
            bubbles[count].radius=params[0];
            bubbles[count].pos = (vec3){params[1],params[2],params[3]};
            bubble_map[count]=count;
            getline(file, line); 
            //std::cout << line <<std::endl;
            std::istringstream iss2(line);
            int i;
            while (iss2>>i){
                bubbles[count].neighboors.push_back(i);
                bubbles[i].neighboors.push_back(count);
            } 
            count++;
        }
    }
    return count;
}

void bubble_sort(bubble arr[], int size) {
    int count = 0;
    while(1)
    {
        count++;
        int switches=0;
        for(int i=0; i<size; i++)
        {
            if(bubble_map[i+1]==-1)break;
            if(arr[bubble_map[i]].radius < arr[bubble_map[i+1]].radius)
            {
                int tmp = bubble_map[i+1];
                bubble_map[i+1] = bubble_map[i];
                bubble_map[i] = tmp;
                switches++;
            }

        }
 //       std::cout<<"Sort iteration: " << count << " Swiwtches: " << switches <<std::endl;
        if(switches==0)break;
    }
}
int main()
{
    int bubblecount    = read_bubble_coords();     
    bubble_sort(bubbles,bubblecount);
 
    int maxid=0;
    for(int i=0; i<bubblecount; i++)
    {
        double maxradius=bubbles[bubble_map[i]].radius;
        int biggest=-1;
        for(int j=0; j<bubbles[bubble_map[i]].neighboors.size(); j++)
        {
            int neighboor = bubbles[bubble_map[i]].neighboors[j];
            //std::cout << neighboor << std::endl ;
            if(bubbles[neighboor].radius>=maxradius && bubbles[neighboor].id>-1)
            {
                maxradius=bubbles[neighboor].radius;
                biggest=neighboor;
            }
        }
        if(biggest>-1)
        {
            bubbles[bubble_map[i]].id=bubbles[biggest].id;
        }
        else{
            bubbles[bubble_map[i]].id=maxid;
            maxid++;
        }
    }
    for(int i=0; i<bubblecount; i++)
    {
        std::cout << bubbles[bubble_map[i]].radius << " " ;
        std::cout << bubbles[bubble_map[i]].pos.x << " " ;
        std::cout << bubbles[bubble_map[i]].pos.y << " " ;
        std::cout << bubbles[bubble_map[i]].pos.z << " " ;
        std::cout << bubbles[bubble_map[i]].id << std::endl ;
    }
}



