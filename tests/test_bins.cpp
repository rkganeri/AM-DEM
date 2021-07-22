#include <string>   
#include <limits>
#include <float.h>

#include "gtest/gtest.h"
#include "Kokkos_Core.hpp"
#include "fmt/core.h"

#include "io_utils.hpp"
#include "global_settings.hpp"
#include "particles.hpp"
#include "bins.hpp"

// wrapping in the namespace so I can omit prepending to all classes/functions
namespace amdem {

// tests the binning algorithm
TEST(bins,setParticleBins) {
    
    // follow same setup as the particles test 
    int argc = 0;
    char **argv = nullptr;
    Kokkos::initialize(argc, argv);

    // wrap all below in braces to automaticall destroy kokkos views
    {

    const int num_particles = 100;
    const double mean_rad = 13.5e-06;
    const double stdev_rad = 4.0e-06;

    // global settings and particle init already tested in test_particles; here we'll just check binning
    amdem::GlobalSettings& global_settings = amdem::GlobalSettings::getInstance(num_particles,mean_rad,stdev_rad);

    auto particles = std::make_unique<amdem::Particles>(num_particles);
    int seed = 2371;
    particles->init(global_settings, seed);

    auto bins = std::make_unique<amdem::Bins>(global_settings);
    bins->setParticleBins(particles, global_settings);

    // num_bins_x = int(length/(4*max_rad)) = int(0.5/(4*23e-03)) = 5
    // num_bins_x = int(width/(4*max_rad)) = int(1.0/(4*23e-03)) = 10
    // num_bins_x = int(height/(4*max_rad)) = int(3.0/(4*23e-03)) = 32
    EXPECT_EQ(bins->num_bins_x_, 5);
    EXPECT_EQ(bins->num_bins_y_, 10);
    EXPECT_EQ(bins->num_bins_z_, 32);
    EXPECT_EQ(bins->num_particles_, num_particles);

    Kokkos::View<double*>::HostMirror h_particle_bin = Kokkos::create_mirror_view(bins->particle_bin_);
    Kokkos::View<double*>::HostMirror h_bins = Kokkos::create_mirror_view(bins->bins_);
    Kokkos::View<double*>::HostMirror h_linked_list = Kokkos::create_mirror_view(bins->linked_list_);
    Kokkos::deep_copy(h_particle_bins, bins->particle_bin_);
    Kokkos::deep_copy(h_bins, bins->bins_);
    Kokkos::deep_copy(h_linked_list, bins->linked_list_);

    // bin_length = length/num_bins_x = (0.5/5) = 0.1 mm
    // bin_width = width/num_bins_y = (1.0/10) = 0.1 mm
    // bin_height = height/num_bins_z = (3.0/32) = 0.09375 mm
    EXPECT_EQ(h_particle_bins(96,0), 3);
    EXPECT_EQ(h_particle_bins(96,1), 1);
    EXPECT_EQ(h_particle_bins(96,2), 29);

    // all particles should be in the top 1/3rd of the box at initialization
    int min_zbin = INT_MAX;
    int max_zbin = 0;
    for (int i=0; i<num_particles; i++) {
        if (h_particle_bins(i,2) < min_zbin) min_zbin = h_particle_bins(i,2);
        if (h_particle_bins(i,2) > max_zbin) max_zbin = h_particle_bins(i,2);
    }
    EXPECT_EQ(min_zbin, 21);
    EXPECT_EQ(max_zbin, 31);

    // TODO: test bins_ and linked_list_ arrays
    } // end wrapper to destroy views
    
    // finalize
    Kokkos::finalize();
}


} // namespace amdem




