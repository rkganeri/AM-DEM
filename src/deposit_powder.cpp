#include <string>  
#include <memory>
#include <cmath>

#include "Kokkos_Core.hpp"
#include "fmt/core.h"

#include "deposit_powder.hpp"
#include "io_utils.hpp"
#include "terminate.hpp"

namespace amdem {

bool depositPowder(std::unique_ptr<Particles>& particles, 
                   const GlobalSettings& global_settings) {

    // we are using an explicit RK-4 time stepping scheme, so our step size is fixed
    double current_time = global_settings.t_start_;
    const double end_time = global_settings.t_end_;
    const double dt = global_settings.dt_;
    int num_time_steps = static_cast<int>(ceil((end_time-current_time)/dt));
    const int plot_step_freq = static_cast<int>(global_settings.plot_time_freq_/dt);
    // with a max speed of roughly 0.2 m/s, each bin having dimesnion of roughly 1e-4,
    // and a time step size of 5e-08, it takes roughly 10k time steps for a particle to 
    // move through a bin... thus with a safety factor of 2 we update only every 5000 steps
    const int rebin_freq = 5000;

    // Kokkos views we create for the RK-4 time stepping
    Kokkos::View<double*[3]> y1_pos("y1_pos",particles->num_particles_);
    Kokkos::View<double*[3]> y1_vel("y1_vel",particles->num_particles_);
    Kokkos::View<double*[3]> psi1("psi1",particles->num_particles_);
    Kokkos::View<double*[3]> y2_pos("y2_pos",particles->num_particles_);
    Kokkos::View<double*[3]> y2_vel("y2_vel",particles->num_particles_);
    Kokkos::View<double*[3]> psi2("psi2",particles->num_particles_);
    Kokkos::View<double*[3]> y3_pos("y3_pos",particles->num_particles_);
    Kokkos::View<double*[3]> y3_vel("y3_vel",particles->num_particles_);
    Kokkos::View<double*[3]> psi3("psi3",particles->num_particles_);
    Kokkos::View<double*[3]> y4_pos("y4_pos",particles->num_particles_);
    Kokkos::View<double*[3]> y4_vel("y4_vel",particles->num_particles_);
    Kokkos::View<double*[3]> psi4("psi4",particles->num_particles_);
    // some more initialization/definitions best done outside the time loop
    typedef Kokkos::MDRangePolicy<Kokkos::Rank<2>> mdrange_policy2;
    Particles* particles_ptr = nullptr;

    // create initial plot state
    plotState(particles, current_time, 0);
    
    // time loop where it all happens!
    for (int istep=1; istep<=num_time_steps; istep++) {

        current_time += dt;
        // copy state data (position and velocity) from state np1 to state n
        particles->updateState();

        // update our bins
        if ((istep % rebin_freq == 0) or (istep == 1)) {
            print(fmt::format("Starting Step {} now, time = {:14.7e}", istep, current_time)
            particles->initBins();
            particles->setParticleBins(global_settings);
        }

        // RK4 step 1
        // We perform an explicit RK-4 time stepping scheme (see Eqn 29 in 2016 paper)
        Kokkos::deep_copy(y1_pos,particles->coordsn_);
        Kokkos::deep_copy(y1_vel,particles->vn_);
        // calculate the particle forces (stored in psi1)
        particles->calcForces(psi1, global_settings, y1_pos, y1_vel);
        // calculate updated position and velocity sub-increments
        updateRK4SubStep(particles, y2_pos, y2_vel, y1_vel, psi1, dt, 1);

        // RK4 step 2
        particles->calcForces(psi2, global_settings, y2_pos, y2_vel);
        updateRK4SubStep(particles, y3_pos, y3_vel, y2_vel, psi2, dt, 2);

        // RK4 step 3
        particles->calcForces(psi3, global_settings, y3_pos, y3_vel);
        updateRK4SubStep(particles, y4_pos, y4_vel, y3_vel, psi3, dt, 3);

        // RK4 step 4 (final step)
        particles->calcForces(psi4, global_settings, y4_pos, y4_vel);
        // the last step is different as we now do a full time step using all intermediary states
        // unique pointers can't be passed into kokkos lambdas so grab a pointer of the object itself
        particles_ptr = particles.get();
        Kokkos::parallel_for("update_state", mdrange_policy2({0,0},{particles_ptr->num_particles_,3}),
                KOKKOS_LAMBDA(int i, int j) {

            // recall a = vdot = F/m 
            particles_ptr->vnp1_(i,j) = particles_ptr->vn_(i,j) + dt/6.0*((psi1(i,j)+2.0*psi2(i,j)
                                        +2.0*psi3(i,j)+psi4(i,j))/particles_ptr->mass_(i));

            particles_ptr->coordsnp1_(i,j) = particles_ptr->coordsn_(i,j) + dt/6.0*(y1_vel(i,j)
                                             +2.0*y2_vel(i,j)+2.0*y3_vel(i,j)+y4_vel(i,j));

        });
        particles_ptr = nullptr; // probably unnecessary but just in case


        // plot results as specified
        if (istep % plot_step_freq == 0) plotState(particles, current_time, istep);

    }

    return true;

}

void updateRK4SubStep(const std::unique_ptr<Particles>& particles,
                      Kokkos::View<double**> y_pos,
                      Kokkos::View<double**> y_vel,
                      const Kokkos::View<double**> y_vel_prev,
                      const Kokkos::View<double**> psi,
                      const double dt, const int sub_step) {

    // Sets the new y_pos and y_vel values for each RK4 time step sub-increment
    // (only valid for sub-steps 1-3)

    const int num_particles = particles->num_particles_;
    double tstep = 0.0;
    if ((sub_step == 1) or (sub_step == 2)) {
        tstep = dt/2.0;
    } else if (sub_step == 3) {
        tstep = dt;
    } else {
        terminateError(fmt::format("sub_step should only be 1, 2, or 3 in updateRK4SubStep. "
            "Current sub_step value is {}.", sub_step));
    }
            
    // we can' access info stored within a unique pointer within kokkos lambdas so we
    // create a pointer to the underlying object
    Particles* particles_ptr = particles.get();

    typedef Kokkos::MDRangePolicy<Kokkos::Rank<2>> mdrange_policy2;
    Kokkos::parallel_for("update_RK4_substep", mdrange_policy2({0,0},{num_particles,3}), 
            KOKKOS_LAMBDA(int i, int j) {

        // recall a = F/m = psi/mass
        y_vel(i,j) = particles_ptr->vn_(i,j) + tstep*psi(i,j)/particles_ptr->mass_(i);
        y_pos(i,j) = particles_ptr->coordsn_(i,j) + tstep*y_vel_prev(i,j);

    });

}

} // namespace amdem



