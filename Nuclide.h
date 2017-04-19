#ifndef _NUCLIDE_HEADER_
#define _NUCLIDE_HEADER_

#include <vector>
#include <string>
#include <memory>

#include "Reaction.h"
#include "Particle.h"

class nuclide {
  private:
    std::string nuclide_name;
    double nuclide_A;
    std::vector< std::shared_ptr< reaction > > rxn;
  public:
     nuclide( std::string label, double a_num ) : nuclide_name(label), nuclide_A(a_num) {};
    ~nuclide() {};

    std::string name() { return nuclide_name; }
    std::vector< std::shared_ptr< reaction > > getReactions() {return rxn;} ;
    double nuc_A() {return nuclide_A;};

    void        addReaction( std::shared_ptr< reaction >);
    double      total_xs(particle* p);
    double      capture_xs(particle* p);
    double      scatter_xs(particle* p);

    std::shared_ptr< reaction > sample_reaction(particle* p);
};


#endif
