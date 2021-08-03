#include <string>   // std::string
#include <fstream>  // std::flush
#include <iostream> // std::cerr, std::endl

#include "bins.hpp"

namespace amdem {

Bins::Bins(const GlobalSettings& gs) 
    : num_bins_x_(static_cast<int>(gs.length_/(4*gs.max_rad_))),
      num_bins_y_(static_cast<int>(gs.width_/(4*gs.max_rad_))),
      num_bins_z_(static_cast<int>(gs.height_/(4*gs.max_rad_))),
      particle_bin_("particle_bin",gs.num_particles_),
      bins_("bins",num_bins_x_,num_bins_y_,num_bins_z_),
      linked_list_("linked_list",gs.num_particles_)
{ /* everything done via initialization */ }

void Bins::initBins() {
    // Kokkos views are initialized to 0 by default but we really need to set the bins_
    // and linked_list_ arrays to -1
    Kokkos::parallel_for("init_linked_list", linked_list_.extent(0), KOKKOS_LAMBDA(int i) {
        linked_list_(i) = -1;
    });

    typedef Kokkos::MDRangePolicy<Kokkos::Rank<3>> mdrange_policy3;
    Kokkos::parallel_for("init_bins", mdrange_policy3({0,0,0},{num_bins_x_,num_bins_y_,num_bins_z_}), 
            KOKKOS_LAMBDA(int i, int j, int k) {
        bins_(i,j,k) = -1;
    });
}




} // namespace amdem


