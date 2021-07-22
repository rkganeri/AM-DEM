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
      particle_bin_("particle_bin",gs.num_particles_),
      bins_("bins",num_bins_x_,num_bins_y_,num_bins_z_),
      linked_list_("linked_list",num_particles_)
{ /* everything done via initialization */ }

void Bins::init() {
    // Kokkos views are initialized to 0 by default but we really need to set the bins_
    // and linked_list_ arrays to -1
    Kokkos::parallel_for("init_linked_list", num_particles_, KOKKOS_LAMBDA(int i) {
        linked_list_(i) = -1;
    });

    typedef Kokkos::MDRangePolicy<Kokkos::Rank<3>> mdrange_policy3;
    Kokkos::parallel_for("init_bins", mdrange_policy3({0,0,0},{num_bins_x_,num_bins_y_,num_bins_z_}), 
            KOKKOS_LAMBDA(int i, int j, int k) {
        bins_(i,j,k) = -1;
    });
}


void Bins::setParticleBins(const std::unique_ptr<Particles>& particles, 
                           const GlobalSettings& gs) {

    // settings the bins easy to do in parallel
    const double bin_length = gs.length_ / num_bins_x_;
    const double bin_width = gs.width_ / num_bins_y_;
    const double bin_height = gs.height_ / num_bins_z_;
    // we need to create a pointer to the underlying data as we cannot capture a unique pointer
    // in the kokkos lambda 
    const Particles* particles_ptr = particles.get();

    // the below could probably be done just as fast without the loop collapse, but meh...
    Kokkos::parallel_for("set_particle_bins", num_particles_, KOKKOS_LAMBDA(int i) {
        // we use the c-style int conversion since we are in a parallel for loop which may need
        // to be compiled with nvcc if using GPUs
        particle_bin_(i,0) = (int) (particles_ptr->coordsnp1_(i,0)/bin_length);
        particle_bin_(i,1) = (int) (particles_ptr->coordsnp1_(i,1)/bin_width);
        particle_bin_(i,2) = (int) (particles_ptr->coordsnp1_(i,2)/bin_height);
    });


    // ok so to create the linked list we need to do it serially on the CPU. however, as we only
    // loop through num_particles once it's not actually that expensive
    Kokkos::View<int**>::HostMirror h_particle_bin = Kokkos::create_mirror_view(particle_bin_);
    Kokkos::View<int***>::HostMirror h_bins = Kokkos::create_mirror_view(bins_);
    Kokkos::View<int*>::HostMirror h_linked_list = Kokkos::create_mirror_view(linked_list_);
    Kokkos::deep_copy(h_particle_bin, particle_bin_);

    // kokkos views are initialized to 0, which is handy for setting the linked list
    for (int i=num_particles_-1; i>-1; i--) {
        h_linked_list(i) = h_bins(h_particle_bin(i,0), h_particle_bin(i,1), h_particle_bin(i,2));
        h_bins(h_particle_bin(i,0), h_particle_bin(i,1), h_particle_bin(i,2)) = i;
    }

    // copy back to device (if needed) and we're done
    Kokkos::deep_copy(bins_, h_bins);
    Kokkos::deep_copy(linked_list_, h_linked_list);
        
}


} // namespace amdem


