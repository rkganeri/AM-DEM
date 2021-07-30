#ifndef AMDEM_PARTICLES_HPP
#define AMDEM_PARTICLES_HPP

#include <string>   
#include <memory>

#include "Kokkos_Core.hpp"

#include "global_settings.hpp"
#include "bins.hpp"

namespace amdem {

class Particles: public Bins {
    // This class stores the force and displacement data for each of the DEM particles.
    // Kokkos views are used to allow for multi-dimensional data access, whereby Kokkos
    // determines the optimal data layout based upon the compilation settings (e.g. LayoutRight
    // for HostSpace or LayoutLeft for CudaSpace)

    private:
        // some data used in Hertzian contact calculations (these get set in init method)
        double estar_;
        double zeta_;
        double mu_fric_;

    public:
        // make the data public so it's easier to access where we need it
        const int num_particles_;

        Kokkos::View<double*> radius_;  // (num_particles)
        Kokkos::View<double*> mass_;
        Kokkos::View<double*> volume_;

        Kokkos::View<double*[3]> vn_;           // velocity at step n - (num_particles,3)
        Kokkos::View<double*[3]> vnp1_;         // at step n plus 1
        Kokkos::View<double*[3]> coordsn_;      // position at step n
        Kokkos::View<double*[3]> coordsnp1_;    // position at step np1

        Kokkos::View<double*[3]> psi_tot_;      // particle forces

        // Methods below:
        // delete the default constructor
        Particles(const int num_particles, const GlobalSettings& gs);
        Particles() = delete;

        void initParticles(const GlobalSettings& global_settings, int seed=0);

        void setParticleBins(const GlobalSettings& global_settings);


        // views are treated as pointers so we capture them by value
        //void calcForces(const std::unique_ptr<Bins>& bins, const GlobalSettings& global_settings,
        void calcForces(const GlobalSettings& global_settings,
                        const Kokkos::View<double**> coords, const Kokkos::View<double**> vel);

        void calcWallForce(Kokkos::View<double**> psi_con, Kokkos::View<double**> psi_fric, 
                           const Kokkos::View<double**> coords, const Kokkos::View<double**> vel,
                           const double wall_plane, const int n_index, const int n_value);

};

} // namespace amdem

#endif // AMDEM_PARTICLES_HPP


