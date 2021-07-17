#ifndef AMDEM_BINS_HPP
#define AMDEM_BINS_HPP

#include <string>   

#include "Kokkos_Core.hpp"

#include "global_settings.hpp"
#include "particles.hpp"

namespace amdem {

class Bins {
    // This class stores the bins information related to where each particle is located

    public:
        const int num_bins_x_;
        const int num_bins_y_;
        const int num_bins_z_;
        const int num_particles_;

        Kokkos::View<double***> bins_;  // (num_particles)
        Kokkos::View<double*> linked_list_;

        // methods - again delete default constructor
        Bins(const GlobalSettings& gs);
        Bins() = delete;

        void setBins(std::unique_ptr<Particles>& particles);

};

} // namespace amdem

#endif // AMDEM_BINS_HPP



