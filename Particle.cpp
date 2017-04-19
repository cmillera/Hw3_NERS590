#include <cmath>
#include <limits>

#include "Particle.h"
#include "Random.h"
#include "Material.h"
#include "Cell.h"
#include "Nuclide.h"

// constructor for a new particle
particle::particle( point p, point d ) : p_pos(p), p_dir(d) {
  p_dir.normalize();
  exist = true;
  p_wgt = 1.0;
  p_cell = nullptr;
  p_erg = 1.0;
  p_time = 0.0;
}

// move the particle along its current trajectory
void particle::move( double s ) {
  p_pos.x += s * p_dir.x;
  p_pos.y += s * p_dir.y;
  p_pos.z += s * p_dir.z;
  
  //increment time, assuming non-relativistic neutron, in mks units, converted to ns
  // (2*1.602E-13)/1.660540E-27 [neutron mass, kg] = 1.9295e+14
  p_time += (( s / 100 )/pow((1.9295E14*p_erg),0.5))*1E9;
}

// scatter particle given input direction cosine cos_t0 = mu0
void particle::scatter( double cos_t0, double A ) {
    
    // sample a random azimuthal angle uniformly
    double azi = 2.0 * std::acos(-1.0) * Urand();
    double cos_azi = std::cos(azi);
    double sin_azi = std::sin(azi);

    // rotate the local particle coordinate system aligned along the incident direction
    // to the global problem (x,y,z) coordinate system 
    double sin_t  = std::sqrt( 1.0 - p_dir.z  * p_dir.z  );
    double sin_t0 = std::sqrt( 1.0 - cos_t0 * cos_t0 );

    point q;

    if ( sin_t > std::numeric_limits<double>::epsilon() * 1000.0 ) {
        double c = sin_t0 / sin_t;
        q.x = p_dir.x * cos_t0 + ( p_dir.x * p_dir.z * cos_azi - p_dir.y * sin_azi ) * c;
        q.y = p_dir.y * cos_t0 + ( p_dir.y * p_dir.z * cos_azi + p_dir.x * sin_azi ) * c;
        q.z = p_dir.z * cos_t0 - cos_azi * sin_t0 * sin_t;
    }
    else {
        // if incident direction along z, reorient axes to avoid division by zero
        sin_t  = std::sqrt( 1.0 -  p_dir.y * p_dir.y );
        double c = sin_t0 / sin_t;
        q.x = p_dir.x * cos_t0 + ( p_dir.x * p_dir.y * cos_azi + p_dir.z * sin_azi ) * c;
        q.y = p_dir.y * cos_t0 - cos_azi * sin_t0 * sin_t;
        q.z = p_dir.z * cos_t0 + ( p_dir.y * p_dir.z * cos_azi - p_dir.x * sin_azi ) * c;
    }

    p_dir.x = q.x;
    p_dir.y = q.y;
    p_dir.z = q.z;
  
  // ADJUST NEUTRON ENERGY

  //std::cout << "A= " <<A<< " mu= "<<cos_t0<<" E_initial= " << p_erg<<std::endl;
  double Erecoil = ((4.0 * A)/(( 1.0 + A )*( 1.0 + A )))* cos_t0 * cos_t0 * p_erg;
  p_erg = p_erg - Erecoil;
  //std::cout << "E_final= "<<p_erg<<std::endl;
   
}

// set the particles life flag to false
void particle::kill() {
  exist = false;
}

// set the particle's direction
void particle::setDirection( point p ) {
  p_dir.x = p.x;
  p_dir.y = p.y;
  p_dir.z = p.z;
}

// adjust the weight by a factor f
void particle::adjustWeight( double f ) {
  p_wgt *= f;
}

// set the cell pointer for efficiency
void particle::recordCell( std::shared_ptr< cell > cel ) {
  p_cell = cel;
}
