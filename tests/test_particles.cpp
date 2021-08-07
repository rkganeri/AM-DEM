#include <string>   
#include <limits>
#include <float.h>
#include <cmath>

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
    amdem::GlobalSettings& global_settings = 
        amdem::GlobalSettings::getInstance(num_particles,mean_rad,stdev_rad);
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


// tests particle initialization
TEST(particles,calcWallForce) {
    
    // initialize Kokkos 
    int argc = 0;
    char **argv = nullptr;
    Kokkos::initialize(argc, argv);

    // wrap all below in braces to automaticall destroy kokkos views
    {

    const int num_particles = 3;
    const double mean_rad = 13.5e-06;
    const double stdev_rad = 4.0e-06;

    // may as well test the global settings initialization here while we're at it
    amdem::GlobalSettings& global_settings = 
        amdem::GlobalSettings::getInstance(num_particles,mean_rad,stdev_rad);

    auto particles = std::make_unique<amdem::Particles>(num_particles, global_settings);

    // init the particles with a given seed so the resulting radii and positions are deterministic
    int seed = 2371;
    particles->initParticles(global_settings, seed);

    Kokkos::View<double*>::HostMirror h_radius = Kokkos::create_mirror_view(particles->radius_);
    Kokkos::View<double*>::HostMirror h_mass = Kokkos::create_mirror_view(particles->mass_);
    Kokkos::deep_copy(h_radius, particles->radius_);
    Kokkos::deep_copy(h_mass, particles->mass_);
    
    EXPECT_FLOAT_EQ(h_radius(0), 1.8026120e-05);
    EXPECT_FLOAT_EQ(h_radius(1), 1.6664266e-05);
    EXPECT_FLOAT_EQ(h_radius(2), 1.0623141e-05);

    EXPECT_FLOAT_EQ(h_mass(0), 1.9510652e-10);
    EXPECT_FLOAT_EQ(h_mass(1), 4/3.*M_PI*pow(1.6664266e-05,3)*7952); // 1.5414290e-10
    EXPECT_FLOAT_EQ(h_mass(2), 3.9932258e-11);

    // estar = 103.4963534964e9
    // rho = 7952
    // zeta = 0.1
    // mu_fric = 0.1

    Kokkos::View<double**>::HostMirror h_coordsn = Kokkos::create_mirror_view(particles->coordsn_);
    Kokkos::View<double**>::HostMirror h_vn = Kokkos::create_mirror_view(particles->vn_);
    Kokkos::deep_copy(h_coordsn, particles->coordsn_);
    Kokkos::deep_copy(h_vn, particles->vn_);

    // vel should be initialized to 0
    EXPECT_FLOAT_EQ(h_vn(0,0), 0.0);
    EXPECT_FLOAT_EQ(h_vn(1,0), 0.0);
    EXPECT_FLOAT_EQ(h_vn(2,0), 0.0);
    EXPECT_FLOAT_EQ(h_vn(2,1), 0.0);
    EXPECT_FLOAT_EQ(h_vn(2,2), 0.0);

    // manually set coords and velocity to test force calculations
    // particle 0 hits bottom wall
    h_coordsn(0,0) = 0.25e-03;
    h_coordsn(0,1) = 0.5e-03;
    h_coordsn(0,2) = 1.7526120e-05; // overlap is 0.5e-06 with z=0
    h_vn(0,0) = 0.1;
    h_vn(0,1) = 0.05;
    h_vn(0,2) = -0.2;
    // particle 1 hits 2 walls simultaneously
    h_coordsn(1,0) = 0.25e-03;
    h_coordsn(1,1) = 1.5664266e-05; // overlap 1.0e-06 with y=0
    h_coordsn(1,2) = 1.6164266e-05; // overlap is 0.5e-06 with z=0
    h_vn(1,0) = 0.1;
    h_vn(1,1) = -0.05;
    h_vn(1,2) = -0.2;
    // copy back to device (if needed)
    Kokkos::deep_copy(particles->coordsn_, h_coordsn);
    Kokkos::deep_copy(particles->vn_, h_vn);

    // set up views for storing force and calculate forces from walls
    Kokkos::View<double*[3]> psi_con("psi_con",particles->num_particles_);
    Kokkos::View<double*[3]> psi_fric("psi_fric",particles->num_particles_);
    // z-dir, top wall
    double wall_plane = global_settings.height_;
    int n_index = 2;  // normal index 
    int n_value = 1;  // normal value (e.g. normal = [0, 0, 1])
    particles->calcWallForce(psi_con, psi_fric, particles->coordsn_, particles->vn_, 
                             wall_plane, n_index, n_value);
    // z-dir, bottom wall
    wall_plane = 0.0;
    n_index = 2;  
    n_value = -1;  // normal value (e.g. normal = [0, 0, -1])
    particles->calcWallForce(psi_con, psi_fric, particles->coordsn_, particles->vn_, 
                             wall_plane, n_index, n_value);
    // x-dir, right wall
    wall_plane = global_settings.length_;
    n_index = 0;  // normal index 
    n_value = 1;  // normal value (e.g. normal = [1, 0, 0])
    particles->calcWallForce(psi_con, psi_fric, particles->coordsn_, particles->vn_, 
                             wall_plane, n_index, n_value);
    // x-dir, left wall
    wall_plane = 0.0;
    n_index = 0;  // normal index 
    n_value = -1;  // normal value (e.g. normal = [-1, 0, 0])
    particles->calcWallForce(psi_con, psi_fric, particles->coordsn_, particles->vn_, 
                             wall_plane, n_index, n_value);
    // y-dir, front wall
    wall_plane = global_settings.width_;
    n_index = 1;  // normal index 
    n_value = 1;  // normal value (e.g. normal = [0, 1, 0])
    particles->calcWallForce(psi_con, psi_fric, particles->coordsn_, particles->vn_, 
                             wall_plane, n_index, n_value);
    // y-dir, back wall
    wall_plane = 0.0;
    n_index = 1;  // normal index 
    n_value = -1;  // normal value (e.g. normal = [0, -1, 0])
    particles->calcWallForce(psi_con, psi_fric, particles->coordsn_, particles->vn_, 
                             wall_plane, n_index, n_value);

    // get views back on host so we can query them
    Kokkos::View<double**>::HostMirror h_psi_con = Kokkos::create_mirror_view(psi_con);
    Kokkos::View<double**>::HostMirror h_psi_fric = Kokkos::create_mirror_view(psi_fric);
    Kokkos::deep_copy(h_psi_con, psi_con);
    Kokkos::deep_copy(h_psi_fric, psi_fric);

    // I calculated these via hand to ensure we're getting the right answers 
    EXPECT_FLOAT_EQ(h_psi_con(0,0), 0.0);
    EXPECT_FLOAT_EQ(h_psi_con(0,1), 0.0);
    EXPECT_NEAR(h_psi_con(0,2), 0.206702392, 1.0e-06);

    EXPECT_FLOAT_EQ(h_psi_con(1,0), 0.0);
    EXPECT_NEAR(h_psi_con(1,1), 0.563208091, 1.0e-06);
    EXPECT_NEAR(h_psi_con(1,2), 0.198780606, 1.0e-06);

    EXPECT_NEAR(h_psi_fric(0,0), -0.03697605, 1.0e-06);
    EXPECT_NEAR(h_psi_fric(0,1), -0.01848802, 1.0e-06);
    EXPECT_FLOAT_EQ(h_psi_fric(0,2), 0);

    EXPECT_NEAR(h_psi_fric(1,0), -0.08593382, 1.0e-06);
    EXPECT_NEAR(h_psi_fric(1,1), 0.01777948, 1.0e-06);
    EXPECT_NEAR(h_psi_fric(1,2), 0.10074973, 1.0e-06);

    } // end wrapper to destroy views
    
    // finalize
    Kokkos::finalize();
}

} // namespace amdem




