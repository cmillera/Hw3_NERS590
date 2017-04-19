#include "Reaction.h"
#include "Particle.h"
#include "Distribution.h"

void  capture_reaction::sample( particle* p, std::stack< particle >* bank, double A ) {
  // kill the particle and leave the bank unmodified
  p->kill();
  //std::cout << "capture"<< std::endl;
}

void  scatter_reaction::sample( particle* p, std::stack< particle >* bank, double A ) {
  // scatter the particle and leave the bank unmodified
  double mu0 = scatter_dist->sample();
  p->scatter( mu0, A );
  //std::cout << "scatter"<< std::endl;
}

void  fission_reaction::sample( particle* p, std::stack< particle >* bank, double A ) {
  // create random number of secondaries from multiplicity distributon and
  // push all but one of them into the bank, and set working particle to the last one
  // if no secondaries, kill the particle
  int n = multiplicity_dist->sample();
  if ( n <= 0 ) {
    p->kill();
  }
  else {
    // bank all but last particle (skips if n = 1)
    for ( int i = 0 ; i < (n - 1) ; i++ ) {
      particle q( p->pos(), isotropic->sample() );
      q.recordCell( p->cellPointer() );
      bank->push( q );
    }
    // set working particle to last one
    particle q( p->pos(), isotropic->sample() );
    q.recordCell( p->cellPointer() );
    *p = q;
  }
}
