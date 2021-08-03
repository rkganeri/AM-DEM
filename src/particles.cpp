#define _USE_MATH_DEFINES
#include <string>   
#include <ctime>
#include <random>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <algorithm>

#include "Kokkos_Core.hpp"
#include "Kokkos_Random.hpp"
#include "fmt/core.h"

#include "particles.hpp"
#include "utilities.hpp"
#include "io_utils.hpp"
#include "terminate.hpp"

namespace amdem {

Particles::Particles(const int num_particles, const GlobalSettings& gs)
    : Bins(gs), // initialize base class
      num_particles_(num_particles),
      radius_("radius",num_particles),
      mass_("mass",num_particles),
      volume_("volume",num_particles),
      vn_("vn",num_particles),  // (num_particles,3)
      vnp1_("vnp1",num_particles),  // (num_particles,3)
      coordsn_("coordsn",num_particles),  // (num_particles,3)
      coordsnp1_("coordsnp1",num_particles)  // (num_particles,3)
{   }  

// initialize particle radii and locations (also calculate mass/volume)
void Particles::initParticles(const GlobalSettings& global_settings, int seed) {

    // unpack local vars
    const double length = global_settings.length_;
    const double width = global_settings.width_;
    const double height = global_settings.height_;
    const double mean_rad = global_settings.mean_rad_;
    const double stdev_rad = global_settings.stdev_rad_;
    const double min_rad = global_settings.min_rad_;
    const double max_rad = global_settings.max_rad_;

    // NOTE: while there are ways to generate the random radii and particle coordinates in parallel on a 
    // per thread basis using the Kokkos_Random.hpp header and the functor template in the utilities.hpp file,
    // it is actually a bit painful.  Since we only need to generate this once at initialization, I will opt
    // to just do this in serial for now using standard c++ template libraries, to make things a bit simpler
    // (and at minimal extra computational expense)

    // We will have to set all values on the host first and then copy them over to the device
    // We accomplish this using mirror views. Note that when compiled for CPU-only architectures, the
    // host mirror is merely a refence to the original data already residing on the host, no extra copy is created.
    Kokkos::View<double*>::HostMirror h_radius = Kokkos::create_mirror_view(this->radius_);
    Kokkos::View<double**>::HostMirror h_coordsn = Kokkos::create_mirror_view(this->coordsn_);

    // generate random normal distribution for particle radii
    std::default_random_engine generator;
    if (seed == 0) {
        generator.seed(std::chrono::system_clock::now().time_since_epoch().count());
    } else {
        generator.seed(seed);
    }
    std::normal_distribution<double> normal_distribution(mean_rad, stdev_rad);
    std::uniform_real_distribution<double> uniform_distribution(0.0, 1.0);   // between 0 and 1, we'll scale this later
    for (int i=0; i<num_particles_; i++) {
        double rad = normal_distribution(generator);
        // min and max are bounded
        if (rad < min_rad) {
            h_radius(i) = min_rad;
        } else if (rad > max_rad) {
            h_radius(i) = max_rad;
        } else {
            h_radius(i) = rad;
        }
        // generate starting coords evenly distributed in x/y and within top 1/3rd of box in z
        h_coordsn(i,0) = h_radius(i) + (length-2*h_radius(i))*uniform_distribution(generator);
        h_coordsn(i,1) = h_radius(i) + (width-2*h_radius(i))*uniform_distribution(generator);
        h_coordsn(i,2) = height - h_radius(i) - (1.0/3.0*height-2*h_radius(i))*uniform_distribution(generator);
    }

    // now check to see if particles overlap, if so we need to generate a new coordinate.
    // this can be expensive as it is an O(n!) check
    int i = num_particles_ - 2;
    int j; 
    double dist;
    while (i > -1) {
        j = i+1;
        while (j < num_particles_) {
            dist = utilities::norm2(h_coordsn(i,0)-h_coordsn(j,0), h_coordsn(i,1)-h_coordsn(j,1), 
                                    h_coordsn(i,2)-h_coordsn(j,2));

            // re-generate center coord if the particles overlap and reset j index as we need to check all again
            if (dist < (h_radius(i)+h_radius(j))) {
                h_coordsn(i,0) = h_radius(i) + (length-2*h_radius(i))*uniform_distribution(generator);
                h_coordsn(i,1) = h_radius(i) + (width-2*h_radius(i))*uniform_distribution(generator);
                h_coordsn(i,2) = height - h_radius(i) - (1.0/3.0*height-2*h_radius(i))*uniform_distribution(generator);
                j = i+1;
            } else {
                j+=1;
            }
        }
        i-=1;
    }

    // copy arrays onto device for remainder of calculations (note that fences are included within beginning + end of 
    // deep copy operations, and if we are doing host only memory then only references/pointers are being passed back and forth,
    // no actual deep copy occurs 
    Kokkos::deep_copy(radius_, h_radius);
    Kokkos::deep_copy(coordsn_, h_coordsn);
    // the below always incurs a deep copy because h_coordsn is not a mirror_view of this->coordsnp1_
    Kokkos::deep_copy(coordsnp1_, h_coordsn);


    // use our first actual kokkos loops now to calculate the mass and volume and set initial velocity
    // the Kokkos::parallel_for loop construct will default to OpenMP threading if compiled CPU-only
    // or create a Cuda/HIP device function call (through lambda functions) if compiled with Cuda or HIP.
    // the KOKKOS_LAMBDA macro captures by value, which is what we want since Kokkos Views are treated
    // as pointers
    double rho = global_settings.rho_;
    Kokkos::parallel_for("calc_mass_vol", num_particles_, KOKKOS_LAMBDA(int i) {
        volume_(i) = 4./3.*M_PI*pow(radius_(i),3.);
        mass_(i) = rho*volume_(i);
    });

    // finally lets set some parameters we use for Hertzian contact (these only need to be calculated once)
    // we assume all particles and the bounding walls are 316L SS with the same mat props
    double e = global_settings.youngs_mod_;
    double nu = global_settings.nu_;
    estar_ = (e*e) / (e*(1-pow(nu,2.0)) + e*(1-pow(nu,2.0)));
    zeta_ = 0.1;    // damping coefficient
    mu_fric_ = 0.1;  // friction coefficient


}


// update state data at the beginning of a new time step
void Particles::updateState() {

    // copy data from state np1 to state n
    Kokkos::deep_copy(coordsn_, coordsnp1_);
    Kokkos::deep_copy(vn_, vnp1_);

}
    


// set the bins based upon particle positions
void Particles::setParticleBins(const GlobalSettings& gs) {

    // settings the bins easy to do in parallel
    const double bin_length = gs.length_ / num_bins_x_;
    const double bin_width = gs.width_ / num_bins_y_;
    const double bin_height = gs.height_ / num_bins_z_;

    // the below could probably be done just as fast without the loop collapse, but meh...
    Kokkos::parallel_for("set_particle_bins", num_particles_, KOKKOS_LAMBDA(int i) {
        // we use the c-style int conversion since we are in a parallel for loop which may need
        // to be compiled with nvcc if using GPUs
        particle_bin_(i,0) = (int) (coordsnp1_(i,0)/bin_length);
        particle_bin_(i,1) = (int) (coordsnp1_(i,1)/bin_width);
        particle_bin_(i,2) = (int) (coordsnp1_(i,2)/bin_height);
    });


    // ok so to create the linked list we need to do it serially on the CPU. however, as we only
    // loop through num_particles once it's not actually that expensive
    Kokkos::View<int**>::HostMirror h_particle_bin = Kokkos::create_mirror_view(particle_bin_);
    Kokkos::View<int***>::HostMirror h_bins = Kokkos::create_mirror_view(bins_);
    Kokkos::View<int*>::HostMirror h_linked_list = Kokkos::create_mirror_view(linked_list_);
    Kokkos::deep_copy(h_particle_bin, particle_bin_);

    for (int i=num_particles_-1; i>-1; i--) {
        h_linked_list(i) = h_bins(h_particle_bin(i,0), h_particle_bin(i,1), h_particle_bin(i,2));
        h_bins(h_particle_bin(i,0), h_particle_bin(i,1), h_particle_bin(i,2)) = i;
    }

    // copy back to device (if needed) and we're done
    Kokkos::deep_copy(bins_, h_bins);
    Kokkos::deep_copy(linked_list_, h_linked_list);
        
}


// calculate particle-particle and particle-wall forces, which get stored in psi_tot
void Particles::calcForces(Kokkos::View<double**> psi_tot, const GlobalSettings& global_settings,
                           const Kokkos::View<double**> coords, const Kokkos::View<double**> vel) {
    
    // create separate views for each of the force components
    // N.B. kokkos views automatically get initialized to 0.0, as we desire
    Kokkos::View<double*[3]> psi_con("psi_con",num_particles_);
    Kokkos::View<double*[3]> psi_fric("psi_fric",num_particles_);
    Kokkos::View<double*[3]> psi_env("psi_env",num_particles_);

    // calculation of some parameters for hertzian contact


    // 1. calculate the particle-wall contact and frictional forces
    // z-dir, top wall
    double wall_plane = global_settings.height_;
    int n_index = 2;  // normal index 
    int n_value = 1;  // normal value (e.g. normal = [0, 0, 1])
    calcWallForce(psi_con, psi_fric, coords, vel, wall_plane, n_index, n_value);
    // z-dir, bottom wall
    wall_plane = 0.0;
    n_index = 2;  
    n_value = -1;  // normal value (e.g. normal = [0, 0, -1])
    calcWallForce(psi_con, psi_fric, coords, vel, wall_plane, n_index, n_value);

    // x-dir, right wall
    wall_plane = global_settings.length_;
    n_index = 0;  // normal index 
    n_value = 1;  // normal value (e.g. normal = [1, 0, 0])
    calcWallForce(psi_con, psi_fric, coords, vel, wall_plane, n_index, n_value);
    // x-dir, left wall
    wall_plane = 0.0;
    n_index = 0;  // normal index 
    n_value = -1;  // normal value (e.g. normal = [-1, 0, 0])
    calcWallForce(psi_con, psi_fric, coords, vel, wall_plane, n_index, n_value);

    // y-dir, front wall
    wall_plane = global_settings.width_;
    n_index = 1;  // normal index 
    n_value = 1;  // normal value (e.g. normal = [0, 1, 0])
    calcWallForce(psi_con, psi_fric, coords, vel, wall_plane, n_index, n_value);
    // y-dir, back wall
    wall_plane = 0.0;
    n_index = 1;  // normal index 
    n_value = -1;  // normal value (e.g. normal = [0, -1, 0])
    calcWallForce(psi_con, psi_fric, coords, vel, wall_plane, n_index, n_value);

    // 2. calculate particle-particle contact and frictional forces 
    // (eqns 2 - 14 in the paper)
    Kokkos::parallel_for("calc_part_part_force", num_particles_, KOKKOS_LAMBDA (int i) {
        // only loop through the neighboring bins when checking for contact
        int h0 = particle_bin_(i,0);
        int h1 = particle_bin_(i,1);
        int h2 = particle_bin_(i,2);
        // grr, C doesn't have a built-in max/min function so i had to write my own
        int min_bin0 = utilities::max(h0-1,0);
        int max_bin0 = utilities::min(h0+2,num_bins_x_);
        int min_bin1 = utilities::max(h1-1,0);
        int max_bin1 = utilities::min(h1+2,num_bins_y_);
        int min_bin2 = utilities::max(h2-1,0);
        int max_bin2 = utilities::min(h2+2,num_bins_z_);
        int j;
        double normal[3], vdiff[3], psi_con_part[3], v_tan0[3], v_tan1[3], v0[3], v1[3];
        double dist, rstar, mstar, delta, delta_dot, damp_coef, vel0_mag, vel1_mag;
        for (int x=min_bin0; x<max_bin0; x++) {
            for (int y=min_bin1; y<max_bin1; y++) {
                for (int z=min_bin2; z<max_bin2; z++) {
                    j = bins_(x,y,z);
                    while ((j >= 0) and (j < num_particles_)) {
                        normal[0] = coords(j,0) - coords(i,0);
                        normal[1] = coords(j,1) - coords(i,1);
                        normal[2] = coords(j,2) - coords(i,2);
                        dist = utilities::norm2(normal, 3);

                        if ( (dist < (radius_(i)+radius_(j))) and (j != i) ) {
                            for (int index=0; index<3; index++) {
                                // normalize normal vector
                                normal[index] = normal[index] / dist;
                                // calculate velocity difference between particles
                                vdiff[index] = vel(j,index) - vel(i,index);
                            }

                            delta = abs(dist - (radius_(i)+radius_(j)));
                            delta_dot = utilities::dotProduct(vdiff, normal, 3);

                            rstar = radius_(i)*radius_(j) / (radius_(i)+radius_(j));
                            mstar = mass_(i)*mass_(j) / (mass_(i)+mass_(j));

                            // now calculate the contact force for this particle-pair
                            damp_coef = 2.0*zeta_*sqrt(2.0*estar_*mstar)*pow(rstar*delta,0.25);
                            for (int index=0; index<3; index++) {
                                psi_con_part[index] = -4.0/3.0*sqrt(rstar)*estar_*pow(delta,1.5)*normal[index]
                                                      + damp_coef*delta_dot*normal[index];

                                psi_con(i,index) += psi_con_part[index];
                                
                                // store these velocities for friction calcs
                                v0[index] = vel(i,index);
                                v1[index] = vel(j,index);
                            }

                            // apply friction force if there is a difference in the particles' tangential vel
                            vel0_mag = utilities::dotProduct(v0, normal, 3);
                            vel1_mag = -utilities::dotProduct(v1, normal, 3);
                            for (int index=0; index<3; index++) {
                                v_tan0[index] = v0[index] - vel0_mag*normal[index];
                                v_tan1[index] = v1[index] - vel1_mag*(-normal[index]);
                                vdiff[index] = v_tan1[index] - v_tan0[index];  // overload this guy
                            }

                            double v_tan_norm = utilities::norm2(vdiff, 3);
                            if (v_tan_norm > 0.0) {
                                double psi_con_norm = utilities::norm2(psi_con_part, 3);
                                psi_fric(i,0) += mu_fric_*psi_con_norm*vdiff[0] / v_tan_norm;
                                psi_fric(i,1) += mu_fric_*psi_con_norm*vdiff[1] / v_tan_norm;
                                psi_fric(i,2) += mu_fric_*psi_con_norm*vdiff[2] / v_tan_norm;
                            }

                        }   // end if

                        j = linked_list_(j);    // go to next particle within this bin in linked list

                    }   // end while

                }   // end bin loop
            }
        }
    });


    // 3.  use Stoke's law to calculate drag (eqn 15 in paper)
    const double viscosity_ar = 2.23e-5; // [kg/m-s or Pa-s] at NTP
    //const double rho_ar = 1.66;  // [kg/m^3] at NTP (20 C, 1 atm) for argon atmosphere
    Kokkos::parallel_for("calc_drag", num_particles_, KOKKOS_LAMBDA (int i) {
        double coef =  -6.0 * M_PI * viscosity_ar * radius_(i);
        psi_env(i,0) = coef * vel(i,0);
        psi_env(i,1) = coef * vel(i,1);
        psi_env(i,2) = coef * vel(i,2);
    });

    // 4. sum up all forces and add gravity contribution (eqns 1, 16 in paper)
    const double g = 9.81;
    Kokkos::parallel_for("sum_forces", num_particles_, KOKKOS_LAMBDA (int i) {
        psi_tot(i,0) = psi_con(i,0) + psi_fric(i,0) + psi_env(i,0);
        psi_tot(i,1) = psi_con(i,1) + psi_fric(i,1) + psi_env(i,1);
        // grav force gets added to z-component
        psi_tot(i,2) = psi_con(i,2) + psi_fric(i,2) + psi_env(i,2) - mass_(i)*g;
    });


}


void Particles::calcWallForce(Kokkos::View<double**> psi_con, Kokkos::View<double**> psi_fric, 
                              const Kokkos::View<double**> coords, const Kokkos::View<double**> vel,
                              const double wall_plane, const int n_index, const int n_value) {

    // quick error/sanity checks
    if ((n_index < 0) or (n_index > 2)) {
        terminateError(fmt::format("Normal index must be 0, 1, or 2 in Particles::calcWallForce, current value is {}",n_index));
    }
    if (abs(n_value) != 1) {
        terminateError(fmt::format("Normal value must be +/- 1 in Particles::calcWallForce, current value is {}",n_value));
    }

    // calculate particle-wall interaction forces (eqns 2 - 14 in paper)
    Kokkos::parallel_for("calc_wall_force", num_particles_, KOKKOS_LAMBDA (int i) {
        // only calculate force if particle overlaps with the plane of the wall
        double dist = abs(coords(i,n_index) - wall_plane);
        if (dist < radius_(i)) {

            double normal[3] = {0};
            normal[n_index] = n_value;
            double v[3];
            v[0] = vel(i,0);
            v[1] = vel(i,1);
            v[2] = vel(i,2);

            double delta_wall = abs(dist - radius_(i));
            double delta_wall_dot = -1.*utilities::dotProduct(v, normal, 3);
            double rstar = radius_(i);
            double mstar = mass_(i);

            double damp_coef = 2.0*zeta_*sqrt(2.0*estar_*mstar)*pow(rstar*delta_wall,0.25);

            double psi_con_wall[3] = {0};
            psi_con_wall[n_index] = -4.0/3.0*sqrt(rstar)*estar_*pow(delta_wall,1.5)*normal[n_index]
                                    + damp_coef*delta_wall_dot*normal[n_index];

            psi_con(i,n_index) += psi_con_wall[n_index];

            // apply friction force if there is tengential velocity
            double v_tan[3];
            double vel_norm_mag = utilities::dotProduct(v, normal, 3);
            v_tan[0] = v[0] - vel_norm_mag*normal[0];
            v_tan[1] = v[1] - vel_norm_mag*normal[1];
            v_tan[2] = v[2] - vel_norm_mag*normal[2];

            if (utilities::norm2(v_tan,3) > 0.0) {
                double psi_con_wall_norm = utilities::norm2(psi_con_wall, 3);
                double v_tan_norm = utilities::norm2(v_tan, 3);
                psi_fric(i,0) += -mu_fric_*psi_con_wall_norm*v_tan[0]/v_tan_norm;
                psi_fric(i,1) += -mu_fric_*psi_con_wall_norm*v_tan[1]/v_tan_norm;
                psi_fric(i,2) += -mu_fric_*psi_con_wall_norm*v_tan[2]/v_tan_norm;
            }
        }

    });

}


} // namespace amdem

