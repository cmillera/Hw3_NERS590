#include <vector>
#include <utility>
#include <memory>
#include <cassert>

#include "Random.h"
#include "Particle.h"
#include "Reaction.h"
#include "Nuclide.h"
#include "Material.h"
#include "Point.h"

// add a nuclide and atom fraction to the material
void material::addNuclide( std::shared_ptr< nuclide > N, double frac ) { 
  nuclides.push_back( std::make_pair( N, frac ) ); 
}

// private utility function that returns sum of atomic fraction * microscopic total xs
// multiply this by atomic density to get macroscopic cross section
// (this function is useful to accelerate sampling the nuclide)
double material::micro_xs(particle* p) {
  double xs = 0.0;
  for ( auto n : nuclides ) { 
    // first is pointer to nuclide, second is atomic fraction
    xs += n.first->total_xs(p) * n.second;
    //std::cout << n.first->total_xs(p) << " " << n.second << std::endl;
  }
  return xs;
}

//same as above but only capture
double material::micro_capture_xs(particle* p) {
  double capture_xs = 0.0;
  for ( auto n : nuclides ) { 
    // first is pointer to nuclide, second is atomic fraction
    capture_xs += n.first->capture_xs(p) * n.second;
  }
  return capture_xs;
}

double material::micro_scatter_xs(particle* p) {
  double scatter_xs = 0.0;
  for ( auto n : nuclides ) { 
    // first is pointer to nuclide, second is atomic fraction
    scatter_xs += n.first->scatter_xs(p) * n.second;
  }
  return scatter_xs;
}

// return the macroscopic cross section
double material::macro_xs(particle* p) {
  return atom_density() * micro_xs(p);
}

// same as above but capture
double material::macro_capture_xs(particle* p) {
  return atom_density() * micro_capture_xs(p);
}

double material::macro_scatter_xs(particle* p) {
  return atom_density() * micro_scatter_xs(p);
}

// randomly sample a nuclide based on total cross sections and atomic fractions
std::shared_ptr< nuclide > material::sample_nuclide(particle* p) {
  double u = micro_xs(p) * Urand();
  //point q =p->pos();
  //std::cout << q.x <<std::endl;
  //std::cout << q.y <<std::endl;
  //std::cout << q.z <<std::endl;
  //std::cout << material_name<<std::endl;
  //std::cout <<micro_xs(p)<<std::endl;
  double s = 0.0;
  for ( auto n : nuclides ) {
    // first is pointer to nuclide, second is atomic fraction
    s += n.first->total_xs(p) * n.second;
    if ( s > u ) { return n.first; }
  }
  assert( false ); // should never reach here
  return nullptr;
}

// function that samples an entire collision: sample nuclide, then its reaction, 
// and finally process that reaction with input pointers to the working particle p
// and the particle bank
void material::sample_collision( particle* p, std::stack<particle>* bank ) {
  // first sample nuclide
  std::shared_ptr< nuclide >  N = sample_nuclide(p);

  // now get the reaction
  std::shared_ptr< reaction > R = N->sample_reaction(p);

  // finally process the reaction
  R->sample( p, bank, N->nuc_A() );
}
