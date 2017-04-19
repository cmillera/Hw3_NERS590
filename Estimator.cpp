#include <cmath>
#include <iostream>

#include "Estimator.h"
#include "Material.h"
#include "Particle.h"

void surface_current_estimator::score( particle* p ) { tally_hist += p->wgt(); }

void counting_estimator::score( particle* p ) { count_hist++; }



void counting_estimator::endHistory() {
  if ( tally.size() < count_hist + 1 ) { tally.resize( count_hist + 1, 0.0 ); }
  tally[ count_hist ] += 1.0;
  nhist++;
  count_hist = 0;
}

double counting_estimator::report() {
  std::cout << name() << std::endl;
  double s1 = 0.0, s2 = 0.0;
  for ( int i = 0 ; i < tally.size() ; i++ ) {
    double p = tally[i] / nhist;
    std::cout << i << " " << p << "   " << std::sqrt( p * ( 1.0 - p ) / nhist ) / p <<  std::endl;
    s1 += p * i;
    s2 += p * i * i;
  }
  std::cout << "   mean = " << s1 << std::endl;
  std::cout << "   var  = " << s2 - s1*s1 << std::endl;
  return std::sqrt( s2-s1*s1 ) / s1;
}

void cell_pathLengthFluxCapture_estimator::score2( particle* p, double path_length, std::shared_ptr< material > medium ) {
  
  double capxs = medium -> macro_capture_xs(p);
    //std::cout <<"path length" <<path_length << std::endl;
    //std::cout << "capxs"<<capxs << std::endl;
    time++;
  tally_hist += p->wgt() * path_length*capxs;
  //std::cout << "tallied particle weight: " <<p->wgt()<<std::endl;
}

void cell_pathLengthFluxScatter_estimator::score2( particle* p, double path_length, std::shared_ptr< material > medium ) {
  
  time_bin = floor(p->time()/bin_width);
  //std::cout << "time " <<p->time()<<std::endl;
  //std::cout << "time_bin " <<time_bin<<std::endl;
  if (time_bin > num_bins)
  {
    p->kill();
  }
  else
  {
    double scatxs = medium -> macro_scatter_xs(p);
        //std::cout <<"path length " <<path_length << std::endl;
        //std::cout << "scatxs "<<scatxs << std::endl;
        //std::cout << "weight "<<p->wgt() << std::endl;
        time++;
    double val = p->wgt() * path_length*scatxs;
            //std::cout << "val "<<val << std::endl;
    count_hist[time_bin] = count_hist[time_bin] + val;
  }
}

void cell_pathLengthFluxScatter_estimator::endHistory() {
  
  for (int i=0; i<tally.size(); i++)
  {
    tally[i] = tally[i]+count_hist[i];
  }
  
  nhist++;
  std::fill(count_hist.begin(), count_hist.end(), 0);
}

double cell_pathLengthFluxScatter_estimator::report() {
  std::cout << name() << std::endl;
  std::cout << nhist << std::endl;
  double s1 = 0.0, s2 = 0.0;
  for ( int i = 0 ; i < tally.size() ; i++ ) {
    double p = tally[i] / nhist;
    std::cout << (i+1)*bin_width << " ns " << p << "   " << std::sqrt( p * ( 1.0 - p ) / nhist ) / p <<  std::endl;
    s1 += p * i;
    s2 += p * i * i;
  }
  std::cout << "   mean = " << s1 << std::endl;
  std::cout << "   var  = " << s2 - s1*s1 << std::endl;
  return std::sqrt( s2-s1*s1 ) / s1;
}