#ifndef _CELL_HEADER_
#define _CELL_HEADER_

#include <string>
#include <vector>
#include <utility>
#include <memory>

#include "Point.h"
#include "Surface.h"
#include "Material.h"
#include "Estimator.h"
#include "Cell.h"

class cell {
  private:
    std::string cell_name;

    std::vector< std::pair< std::shared_ptr< surface >, int > > surfaces;
    std::shared_ptr< material > cell_material;

    std::vector< std::shared_ptr< estimator > > cell_estimators;

    double importance;
  public:
     cell( std::string label ) : cell_name(label) { importance = 1.0; };
    ~cell() {};

    std::string name() { return cell_name; };

    void setMaterial( std::shared_ptr< material > M ) { cell_material = M; };
    std::shared_ptr< material > getMaterial() { return cell_material; }

    void setImportance( double imp ) { importance = imp; };
    double getImportance() { return importance; }

    void addSurface( std::shared_ptr< surface > S, int sense );

    void attachEstimator( std::shared_ptr< estimator > E ) { cell_estimators.push_back( E ); };
    
    //void cellEvent(particle* p, double track_length, std::shared_ptr< material > medium){
        //for ( auto e : cell_estimators ) {e->score2(p, track_length, medium);}
    //}

    //virtual void scoreEstimators( particle* p, double track_length, std::shared_ptr< material > medium) final {
    //void scoreEstimators( particle* p, double track_length, std::shared_ptr< material > medium) {
    //  for ( auto e : cell_estimators ) { e->score( p, track_length, medium); }
    // }
    
    bool testPoint( point p );
    std::pair< std::shared_ptr< surface >, double > surfaceIntersect( ray r );

    double macro_xs(particle* p) {
      if ( cell_material ) { return getMaterial()->macro_xs(p); }
      else { return 0.0; }
    };
    
    double macro_capture_xs(particle* p) {
      if ( cell_material ) { return getMaterial()->macro_capture_xs(p); }
      else { return 0.0; }
    };

    void moveParticle( particle* p, double s );
    void sampleCollision( particle* p, std::stack<particle>* bank );
//    std::shared_ptr< cell > nextCell( particle* p );
};

#endif
