#include <cmath>
#include <cassert>
#include <stack>

#include "Cell.h"
#include "VarianceReduction.h"

void VarianceReduction :: roulette(particle* p, double impRatio){

    //std::cout << "roulette"<< std::endl;
    if (Urand() < (impRatio))
    {
        p->adjustWeight(1.0/impRatio);
    }
    else{p->kill();}
};

void VarianceReduction :: split(particle* p, double impRatio, stack<particle>* bank){

    //std::cout << "split"<< std::endl;
    double n = floor(impRatio + Urand());
        for (int j = 1; j<=n; j++) // make this <= to n, not <n
        {
            particle q( p->pos(), p->dir() );
            q.recordCell( p->cellPointer() );
            q.adjustWeight(p->wgt()/n);
            bank->push( q );
        }
        p->kill();

};