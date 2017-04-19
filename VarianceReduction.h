#ifndef _VARIANCE_REDUCTION_HEADER_
#define _VARIANCE_REDUCTION_HEADER_

#include "Cell.h"

using namespace std;

class VarianceReduction {
    private:
    
    public:
    VarianceReduction( shared_ptr<cell> current_cell, double impRatio, particle* p, stack<particle>* bank){
        if (current_cell -> getImportance() == 0.0)
        {
            //std::cout << "kill"<< std::endl;
            p->kill();
        }
        else if (impRatio < 1.0)
        {
            roulette(p, impRatio);
        }
        else if (impRatio > 1.0)
        {
            split(p, impRatio, bank);
        }

    };
    ~VarianceReduction() {};
    
    void roulette(particle* p, double impRatio);
    void split(particle* p, double impRatio, stack<particle>* bank);

};
#endif