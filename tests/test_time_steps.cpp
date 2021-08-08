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
#include "deposit_powder.hpp"

namespace amdem {

// tests RK4 time stepping scheme
TEST(time_steps,updateRK4) {
    // we will re-use the particle positions and forces calculated in test(particles,calcForces) 
    // to ensure the new coordinates and positions come out correctly given those forces

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
    Kokkos::View<double**>::HostMirror h_coordsn = Kokkos::create_mirror_view(particles->coordsn_);
    Kokkos::View<double**>::HostMirror h_vn = Kokkos::create_mirror_view(particles->vn_);
    Kokkos::deep_copy(h_radius, particles->radius_);
    Kokkos::deep_copy(h_coordsn, particles->coordsn_);
    Kokkos::deep_copy(h_vn, particles->vn_);

    // manually set coords and velocity to test force calculations
    h_coordsn(0,0) = 0.25e-03;
    h_coordsn(0,1) = 0.5e-03;
    h_coordsn(0,2) = 1.0e-03; // overlap is 0.5e-06 with z=0
    h_vn(0,0) = 0.1;
    h_vn(0,1) = 0.05;
    h_vn(0,2) = -0.2;
    h_coordsn(1,0) = 0.25e-03;
    h_coordsn(1,1) = 0.5e-03 + h_radius(0)+h_radius(1) - 1.0e-06; // overlap 1.0e-06 
    h_coordsn(1,2) = 1.0e-03; // overlap is 0.5e-06 with z=0
    h_vn(1,0) = 0.1;
    h_vn(1,1) = -0.05;
    h_vn(1,2) = -0.2;
    h_coordsn(2,0) = 0.25e-03 + h_radius(0)+h_radius(2) - 2.0e-06;  // overlap 2e-6
    h_coordsn(2,1) = 0.5e-03;
    h_coordsn(2,2) = 1.0e-03;
    h_vn(2,0) = 0;
    h_vn(2,1) = 0;
    h_vn(2,2) = 0;

    // copy back to device (if needed)
    Kokkos::deep_copy(particles->coordsn_, h_coordsn);
    Kokkos::deep_copy(particles->vn_, h_vn);

    // now set the bins
    particles->initBins();
    particles->setParticleBins(global_settings);

    // set up views for storing force and calculate forces from walls
    Kokkos::View<double*[3]> psi_tot("psi_tot",particles->num_particles_);

    particles->calcForces(psi_tot, global_settings, particles->coordsn_, particles->vn_);

    // we know the force values from test(particles,calcForces), now let's see if the position/vel
    // get updated correctly
    Kokkos::View<double*[3]> y1_pos("y1_pos",particles->num_particles_);
    Kokkos::View<double*[3]> y1_vel("y1_vel",particles->num_particles_);
    Kokkos::View<double*[3]> y2_pos("y2_pos",particles->num_particles_);
    Kokkos::View<double*[3]> y2_vel("y2_vel",particles->num_particles_);
    Kokkos::View<double*[3]> y3_pos("y3_pos",particles->num_particles_);
    Kokkos::View<double*[3]> y3_vel("y3_vel",particles->num_particles_);
    Kokkos::View<double*[3]> y4_pos("y4_pos",particles->num_particles_);
    Kokkos::View<double*[3]> y4_vel("y4_vel",particles->num_particles_);
    Kokkos::deep_copy(y1_pos,particles->coordsn_);
    Kokkos::deep_copy(y1_vel,particles->vn_);

    updateRK4SubStep(particles, y2_pos, y2_vel, y1_vel, psi_tot, global_settings.dt_, 1);

    Kokkos::View<double**>::HostMirror h_y2_pos = Kokkos::create_mirror_view(y2_pos);
    Kokkos::View<double**>::HostMirror h_y2_vel = Kokkos::create_mirror_view(y2_vel);
    Kokkos::deep_copy(h_y2_pos, y2_pos);
    Kokkos::deep_copy(h_y2_vel, y2_vel);

    EXPECT_NEAR(h_y2_vel(0,0), -1.2918698E+02, 1.0e-04);
    EXPECT_NEAR(h_y2_vel(0,1), -5.8234917E+01, 1.0e-05);
    EXPECT_NEAR(h_y2_vel(0,2), 2.4885358E+01, 1.0e-05);

    EXPECT_NEAR(h_y2_vel(1,0), 9.9999886E-02, 1.0e-08);
    EXPECT_NEAR(h_y2_vel(1,1), 6.5786235E+01, 1.0e-05);
    EXPECT_NEAR(h_y2_vel(1,2), -2.0000002E-01, 1.0e-07);

    EXPECT_NEAR(h_y2_pos(0,0), 2.5000250E-04, 1.0e-10);
    EXPECT_NEAR(h_y2_pos(0,1), 5.0000125E-04, 1.0e-10);
    EXPECT_NEAR(h_y2_pos(0,2), 9.9999500E-04, 1.0e-10);

    EXPECT_NEAR(h_y2_pos(2,0), 2.7664926E-04, 1.0e-10);
    EXPECT_NEAR(h_y2_pos(2,1), 5.0000000E-04, 1.0e-10);
    EXPECT_NEAR(h_y2_pos(2,2), 1.0000000E-03, 1.0e-09);

    // now let's do the next 2 steps to make sure we still get what we expect
    updateRK4SubStep(particles, y3_pos, y3_vel, y2_vel, psi_tot, global_settings.dt_, 2);
    updateRK4SubStep(particles, y4_pos, y4_vel, y3_vel, psi_tot, global_settings.dt_, 3);

    // finally copy these views back and query them
    Kokkos::View<double**>::HostMirror h_y3_pos = Kokkos::create_mirror_view(y3_pos);
    Kokkos::View<double**>::HostMirror h_y3_vel = Kokkos::create_mirror_view(y3_vel);
    Kokkos::deep_copy(h_y3_pos, y3_pos);
    Kokkos::deep_copy(h_y3_vel, y3_vel);
    Kokkos::View<double**>::HostMirror h_y4_pos = Kokkos::create_mirror_view(y4_pos);
    Kokkos::View<double**>::HostMirror h_y4_vel = Kokkos::create_mirror_view(y4_vel);
    Kokkos::deep_copy(h_y4_pos, y4_pos);
    Kokkos::deep_copy(h_y4_vel, y4_vel);

    EXPECT_NEAR(h_y3_vel(0,0), -1.2918698E+02, 1.0e-04);
    EXPECT_NEAR(h_y3_vel(1,1), 6.5786235E+01, 1.0e-05);
    EXPECT_NEAR(h_y3_vel(2,2), -1.2256549E+02, 1.0e-04);

    EXPECT_NEAR(h_y3_pos(0,1), 4.9854413E-04, 1.0e-10);
    EXPECT_NEAR(h_y3_pos(1,2), 9.9999500E-04, 1.0e-10);
    EXPECT_NEAR(h_y3_pos(2,0), 2.9244146E-04, 1.0e-10);

    EXPECT_NEAR(h_y4_vel(0,0), -2.5847395E+02, 1.0e-04);
    EXPECT_NEAR(h_y4_vel(1,1), 1.3162247E+02, 1.0e-04);
    EXPECT_NEAR(h_y4_vel(2,1), 6.1282747E+01, 1.0e-05);

    EXPECT_NEAR(h_y4_pos(0,2), 1.0012443E-03, 1.0e-09);
    EXPECT_NEAR(h_y4_pos(1,0), 2.5000500E-04, 1.0e-10);
    EXPECT_NEAR(h_y4_pos(2,2), 9.9387173E-04, 1.0e-10);


    } // end wrapper to destroy views
    
    // finalize
    Kokkos::finalize();
}

} // namespace amdem
