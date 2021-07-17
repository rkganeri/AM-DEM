#include <string>   // std::string
#include <fstream>  // std::flush
#include <iostream> // std::cerr, std::endl

#include "bins.hpp"

namespace amdem {

Bins::Bins(const GlobalSettings& gs) 
    : num_bins_x_(static_cast<int>(gs.length_/(4*gs.max_rad_))),
      num_bins_y_(static_cast<int>(gs.width_/(4*gs.max_rad_))),
      num_bins_z_(static_cast<int>(gs.height_/(4*gs.max_rad_))),
      num_particles_(gs.num_particles_),
      bins_("bins",num_bins_x_,num_bins_y_,num_bins_z_),
      linked_list_("linked_list",num_particles_)
{ /* everything done via initialization */ }



} // namespace amdem


