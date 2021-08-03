#include <string>   
#include <limits>
#include <float.h>

#include "gtest/gtest.h"
#include "Kokkos_Core.hpp"
#include "fmt/core.h"

#include "io_utils.hpp"
#include "global_settings.hpp"
#include "particles.hpp"

// wrapping in the namespace so I can omit prepending to all classes/functions
namespace amdem {

// tests particle initialization
TEST(particles,initParticles) {
    
    // initialize Kokkos 
    int argc = 0;
    char **argv = nullptr;
    Kokkos::initialize(argc, argv);

    // wrap all below in braces to automaticall destroy kokkos views
    {

    const int num_particles = 100;
    const double mean_rad = 13.5e-06;
    const double stdev_rad = 4.0e-06;

    // may as well test the global settings initialization here while we're at it
    amdem::GlobalSettings& global_settings = amdem::GlobalSettings::getInstance(num_particles,mean_rad,stdev_rad);
    EXPECT_FLOAT_EQ(global_settings.rho_,7952);
    EXPECT_FLOAT_EQ(global_settings.height_,3.0e-03);
    EXPECT_FLOAT_EQ(global_settings.mean_rad_,mean_rad);
    EXPECT_FLOAT_EQ(global_settings.stdev_rad_,stdev_rad);
    EXPECT_EQ(global_settings.num_particles_,num_particles);

    auto particles = std::make_unique<amdem::Particles>(num_particles, global_settings);
    EXPECT_EQ(particles->radius_.extent(0), num_particles);
    EXPECT_EQ(particles->vn_.extent(1), 3);

    // init the particles with a given seed so the resulting radii and positions are deterministic
    int seed = 2371;
    particles->initParticles(global_settings, seed);

    // create host mirrors in case this is performed on the device
    Kokkos::View<double*>::HostMirror h_radius = Kokkos::create_mirror_view(particles->radius_);
    Kokkos::View<double*>::HostMirror h_volume = Kokkos::create_mirror_view(particles->volume_);
    Kokkos::View<double*>::HostMirror h_mass = Kokkos::create_mirror_view(particles->mass_);
    Kokkos::View<double**>::HostMirror h_coordsn = Kokkos::create_mirror_view(particles->coordsn_);
    Kokkos::deep_copy(h_radius, particles->radius_);
    Kokkos::deep_copy(h_volume, particles->volume_);
    Kokkos::deep_copy(h_mass, particles->mass_);
    Kokkos::deep_copy(h_coordsn, particles->coordsn_);

    double min_rad = DBL_MAX;
    double max_rad = 0.0;
    double min_z = DBL_MAX;
    double max_z = 0.0;
    for (int i=0; i<num_particles; i++) {
        if (h_radius(i) < min_rad) min_rad = h_radius(i);
        if (h_radius(i) > max_rad) max_rad = h_radius(i);
        if (h_coordsn(i,2) < min_z) min_z = h_coordsn(i,2);
        if (h_coordsn(i,2) > max_z) max_z = h_coordsn(i,2);
    }
    // ensure we stay within our set min/max radii
    EXPECT_FLOAT_EQ(min_rad, global_settings.min_rad_);
    EXPECT_FLOAT_EQ(max_rad, global_settings.max_rad_);
    // ensure all particles start in top 1/3rd of box
    EXPECT_TRUE(min_z > (global_settings.height_*2./3.0+global_settings.min_rad_));
    EXPECT_TRUE(max_z < (global_settings.height_-global_settings.min_rad_));

    // now check some of the deterministic results 
    EXPECT_FLOAT_EQ(h_radius(13), 1.8209942501009266e-05);
    EXPECT_FLOAT_EQ(h_coordsn(96,0), 0.0003729775063670094);
    EXPECT_FLOAT_EQ(h_coordsn(96,1), 0.00014433068736764189);
    EXPECT_FLOAT_EQ(h_coordsn(96,2), 0.0027633519);

    EXPECT_FLOAT_EQ(h_volume(13), 2.529381471E-14);
    EXPECT_FLOAT_EQ(h_mass(13), 2.529381471E-14*7952);

    } // end wrapper to destroy views
    
    // finalize
    Kokkos::finalize();
}


} // namespace amdem




