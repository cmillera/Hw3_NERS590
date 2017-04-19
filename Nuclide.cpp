#include <vector>
#include <memory>
#include <cassert>

#include "Random.h"
#include "Nuclide.h"

// add a new reaction to the current nuclide
void nuclide::addReaction( std::shared_ptr< reaction > R ) { rxn.push_back( R ); }

// return the total microscopic cross section
double nuclide::total_xs(particle* p) {
  double xs = 0.0;
  for ( auto r : rxn ) 
  { 
    std::shared_ptr<distribution<double>> XSdist = r->xs(p);
    xs += XSdist -> sample2(p);
  }
  return xs;
}

// return the capture microscopic cross section
double nuclide::capture_xs(particle* p) {
  double cap_xs = 0.0;
  for ( auto r : rxn ) {             
    if (r->name() == "capture")
    {    
        std::shared_ptr<distribution<double>> XSdist = r->xs(p);
        cap_xs += XSdist -> sample2(p); 
    }
  }
  return cap_xs;
}

// return the capture microscopic cross section
double nuclide::scatter_xs(particle* p) {
  double scat_xs = 0.0;
  for ( auto r : rxn ) {             
    if (r->name() == "scatter")
    {    
        std::shared_ptr<distribution<double>> XSdist = r->xs(p);
        scat_xs += XSdist -> sample2(p); 
    }
  }
  return scat_xs;
}

// randomly sample a reaction type from this nuclide
std::shared_ptr< reaction > nuclide::sample_reaction(particle* p) {
  double u = total_xs(p) * Urand();
  double s = 0.0;
  for ( auto r : rxn ) {
    std::shared_ptr<distribution<double>> XSdist = r->xs(p);
    s += XSdist -> sample2(p);
    if ( s > u ) { return r; }
  }
  assert( false ); // should never reach here
  return nullptr;
}
