#ifndef AMDEM_BINS_HPP
#define AMDEM_BINS_HPP

#include <string>   
#include <memory>

#include "Kokkos_Core.hpp"

#include "global_settings.hpp"

namespace amdem {

class Bins {
    // This class stores the bins information related to where each particle is located

    public:
        const int num_bins_x_;
        const int num_bins_y_;
        const int num_bins_z_;

        Kokkos::View<int*[3]> particle_bin_; // (num_particles,3)
        Kokkos::View<int***> bins_;          // (num_bins_x,num_bins_y,num_bins_z)
        Kokkos::View<int*> linked_list_;     // (num_particles)

        // methods - again delete default constructor
        Bins(const GlobalSettings& gs);
        Bins() = delete;

        void initBins();

};

} // namespace amdem

#endif // AMDEM_BINS_HPP



