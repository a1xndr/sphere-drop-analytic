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

const int MAX_PORES = 2000000;
const int MAX_TRIES = 100000;
const int X_MAX = 10;
const int Y_MAX = 10;
const int Z_MAX = 10;
const double STEP_SIZE = 0.01;
const double R_STEP_SIZE = 0.001;



struct pore {
    int id=0;
    std::vector<int> bubble_ids;
    double volume=0;
    double count=0;
};

std::vector<pore> pores;
int pores_map[MAX_PORES] = {-1};
double max_volume = 0;
double total_volume = 0;
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
            int paramcount=0;
            int poreid;
            while (iss>>d)
            {
                params[paramcount]=d;
                paramcount++;
                poreid=params[4];
                //std::cout << poreid <<"\n" ;
            }
            //std::cout << pores_map[poreid] <<"\n" ;
            if(pores_map[poreid]==-1)
            {
                pore p;
                p.id= poreid;
                p.bubble_ids.push_back(count);
                p.volume += (4/3)*PI*pow(params[0],3);
                total_volume+=(4/3)*PI*pow(params[0],3);
                if(p.volume>max_volume)max_volume=p.volume;
                p.count++;
                pores.push_back(p);
                pores_map[poreid]=num_pores;
                num_pores++;
                //std::cout << "incremented \n" ;
            }
            else 
            {
                pores[pores_map[poreid]].bubble_ids.push_back(count);
                pores[pores_map[poreid]].volume+= (4/3)*PI*pow(params[0],3);
                total_volume+=(4/3)*PI*pow(params[0],3);
                if(pores[pores_map[poreid]].volume>max_volume)max_volume=pores[pores_map[poreid]].volume;
                pores[pores_map[poreid]].count++;
            }
            count++;
        }
    }
    return num_pores;
}

int main()
{

    for(int i=0; i<MAX_PORES; i++)
    {   
        pores_map[i]=-1;
    }
    int porecount    = read_pores();
    double volint=max_volume/50;
    double volume_frac[50];
    for(int i=0; i<porecount; i++)
    {
        if(pores[i].count>2)
        {
        //std::cout << pores[i].id << " " ;
        //std::cout << pores[i].volume << " " ;
        //std::cout << pores[i].count << std::endl ;
        std::cout << (int)(pores[i].volume*50/max_volume) << std::endl ;
        volume_frac[(int)(pores[i].volume*50/max_volume)]+=pores[i].volume/total_volume;
        }
    }
    for(int i=0; i<50; i++)
    {
        std::cout << i*volint << " " << volume_frac[i] << std::endl;
    }
}



