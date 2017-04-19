#include <string>
#include <memory>
#include <vector>
#include <iostream>
#include <cassert>
#include <stack>
#include <utility>

#include "pugixml.hpp"
#include "Distribution.h"
#include "Reaction.h"
#include "Nuclide.h"
#include "Material.h"
#include "Surface.h"
#include "Cell.h"
#include "Source.h"
#include "Particle.h"
#include "Random.h"
#include "VarianceReduction.h"

using namespace std;

template< typename T >
std::shared_ptr< T > findByName( std::vector< std::shared_ptr< T > > vec, std::string name ) {
  for ( auto v : vec ) {
    if ( v->name() == name ) { return v; }
  }
  return nullptr;
}

int main(){
#include "Input.h"
int time;

//TRANSPORT////////////////////////////////////
for (int i=1; i<= nhist; i++)
{
    stack<particle> bank = src -> sample();
    
    while (! bank.empty())
    {  
        particle p = bank.top();   
        bank.pop();
            //cout <<"history: "<<i<<endl;
            //cout << "bank size: " << bank.size()<<endl;
            //cout << "particle weight: " << p.wgt()<<endl;
            //cout <<endl;     
      
        shared_ptr<cell> current_cell;
        //figure out which cell the particle is in
        for (auto c:cells)
        { 
            if (c-> testPoint(p.pos())){current_cell = c; break;}
        }
        while(p.alive())
        { 
            //sample distance to collision in the cell
            exponential_distribution expDist("CollDist ",current_cell ->macro_xs(&p));
            
            double dColl = expDist.sample();
            
            //sample closest surface, and distance to it
            pair <shared_ptr<surface>,double> intersect = current_cell -> surfaceIntersect(p.getRay());

            if (intersect.second < dColl) // Surface Crossing
            {
                double oldImp = current_cell -> getImportance(); // get old imp. before moving
                
                current_cell -> moveParticle(&p,intersect.second);
                intersect.first -> crossSurface(&p);
                
                for (auto c:cells)
                { 
                    if (c-> testPoint(p.pos())){current_cell = c; break;}
                }
                p.recordCell(current_cell);
                
                double impRatio =(current_cell -> getImportance())/oldImp; // get ratio of new to old imp
                VarianceReduction variance_reduction(current_cell, impRatio, &p, &bank);

            }       
            else //Collision
            {
                current_cell -> moveParticle(&p,dColl);
                current_cell -> sampleCollision(&p, &bank); //REACTION ROUTINE
            }
            time++;

        }

    }    
    
    for(auto e:estimators)
    {
        e -> endHistory();
    }
}
cout <<endl;
for(auto e:estimators)
{
    double relUnc = e -> report();
    //cout <<"rel_unc"<<relUnc<<endl;
    cout <<"Figure of Merit: "<< 1.0/(relUnc*relUnc*time) <<endl;
}
///////////////////////////////////////////////
//cout << "time in number of tracks: "<<time<<endl;
cout <<endl;

  return 0;
}


