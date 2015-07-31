/*
 * =====================================================================================
 *
 *       Filename:  pore-statistics.cpp
 *
 *    Description:  Count pores, their volume etc
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

int MAX_PORES = 2000000;
const int MAX_TRIES = 100000;
const int X_MAX = 10;
const int Y_MAX = 10;
const int Z_MAX = 10;
const double STEP_SIZE = 0.01;
const double R_STEP_SIZE = 0.001;



struct pore {
    int id=0;
    vector<int> bubble_ids;
    double volume=0;
    double count=0;
}

vector<pore> pores;
int pores_map[MAX_PORES] = -1;


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  read_pores
 *  Description:  Reads pores from a file of bubbles. Need to make this filename indep
 * =====================================================================================
 */
int read_pores(){
    std::ifstream file("bubbles_with_pores.list");
    int count = 0;
    int num_pores=0; //The current number of registered pores;
    if(file.is_open())
    {
        std::string line;
        while(getline(file, line))
        {
            std::istringstream iss(line);
            double d;
            double params[5];
            int paramcount=0;;
            while (iss>>d)
            {
                params[paramcount]=d;
                paramcount++;
            }
            if(pores_map[params[4]]==-1)
            {
                pore p;
                p.id= params[4];
                p.bubble_ids.push_back(count);
                p.volume += (4/3)*PI*pow(params[0],3);
                p.count++;
                pores.push_back(p);
                pores_map[params[4]]=num_pores;
                num_pores++;

            }
            else 
            {
                pores[pores_map[params[4]]].bubble_ids.push_back(count);
                pores[pores_map[params[4]]].volume+= (4/3)*PI*pow(params[0],3);
                pores[pores_map[params[4]]].count++;
            }
            count++;
        }
    }
    return num_pores;
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



