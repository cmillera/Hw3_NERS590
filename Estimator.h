#ifndef _ESTIMATOR_HEADER_
#define _ESTIMATOR_HEADER_

#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <cassert>
#include <typeinfo>

#include "Particle.h"
#include "Material.h"
#include "Reaction.h"

class estimator {
  private:
    std::string estimator_name;
  protected:
    unsigned long long nhist;
    unsigned long long time;
  public:
     estimator( std::string label ) : estimator_name(label) {};
    ~estimator() {};

    virtual std::string name() final { return estimator_name; };

    virtual void score( particle* ) = 0;
    virtual void score2( particle*, double, std::shared_ptr< material > )=0;
    template< typename T >
    void score( particle*, T ) { assert(false); };

    

    virtual void endHistory()       = 0;
    virtual double report()           = 0;
};

class single_valued_estimator : public estimator {
  private:

  protected:
    double tally_hist, tally_sum, tally_squared;
  public:
     using estimator::score;

     single_valued_estimator(std::string label ) : estimator(label) { 
       nhist         = 0;
       tally_hist    = 0.0;   
       tally_sum     = 0.0; 
       tally_squared = 0.0;
     };
    ~single_valued_estimator() {};

     virtual void endHistory()    final { 
       nhist++;
       tally_sum     += tally_hist;
       tally_squared += tally_hist * tally_hist;
       tally_hist = 0.0; }

     virtual void score( particle* ) = 0;
     virtual void score2( particle*, double, std::shared_ptr< material > )=0;

     virtual double report() =0;

};

class surface_current_estimator : public single_valued_estimator {
  private:

  public:
     surface_current_estimator( std::string label ) : single_valued_estimator(label) {};
    ~surface_current_estimator() {};

    double report() {
       double mean = tally_sum / nhist;
       double var  = ( tally_squared / nhist - mean*mean ) / nhist;
       std::cout << name() << "   " << mean << "   " << std::sqrt( var ) / mean << std::endl;
       return std::sqrt( var ) / mean;
     };
    void score( particle* );
    virtual void score2( particle*, double, std::shared_ptr< material > ){};
};

class counting_estimator : public estimator {
  private:
    int count_hist;
    std::vector< double > tally;
  public:
     counting_estimator( std::string label ) : estimator(label) { count_hist = 0; };
    ~counting_estimator() {};

    void score( particle* );
    virtual void score2( particle*, double, std::shared_ptr< material > ){};
    void endHistory();
    double report();
};


// absorption scalar flux estimator in a cell
class cell_pathLengthFluxCapture_estimator : public single_valued_estimator {
  private:

  public:
     cell_pathLengthFluxCapture_estimator( std::string label) : 
       single_valued_estimator(label) {};
    ~cell_pathLengthFluxCapture_estimator() {};

    void score2( particle*, double, std::shared_ptr< material > );
    virtual void score( particle* ) {};
    
    double report() {
       double mean = tally_sum / nhist;
       //std::cout<<tally_sum<<std::endl;
       //std::cout<<nhist<<std::endl;
       double var  = ( tally_squared / nhist - mean*mean ) / nhist;
       std::cout << name() << "   " << mean << "   " << std::sqrt( var ) / mean << std::endl;
       return std::sqrt( var ) / mean;
     };
    
};

// scatter scalar flux estimator in a cell
class cell_pathLengthFluxScatter_estimator : public estimator {
  private:
    double bin_width = 5.0;
    int num_bins = 20;
    double time_bin;
    std::vector< double > count_hist;
    std::vector< double > tally;
  public:
     cell_pathLengthFluxScatter_estimator( std::string label) : estimator(label)
     {
        nhist = 0;
        count_hist.resize(num_bins, 0.0);
        tally.resize(num_bins, 0.0);
     };
    ~cell_pathLengthFluxScatter_estimator() {};

    void score( particle* ) {};
    void score2( particle*, double, std::shared_ptr< material > );
    void endHistory();
    double report();
    
};

#endif
