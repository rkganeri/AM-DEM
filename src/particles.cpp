#define _USE_MATH_DEFINES
#include <string>   
#include <ctime>
#include <random>
#include <chrono>
#include <cmath>

#include "Kokkos_Core.hpp"
#include "Kokkos_Random.hpp"

#include "particles.hpp"
#include "utilities.hpp"

namespace amdem {

Particles::Particles(const int num_particles)
    : num_particles_(num_particles),
      radius_("radius",num_particles),
      mass_("mass",num_particles),
      volume_("volume",num_particles),
      vn_("vn",num_particles),  // (num_particles,3)
      vnp1_("vnp1",num_particles),  // (num_particles,3)
      coordsn_("coordsn",num_particles),  // (num_particles,3)
      coordsnp1_("coordsnp1",num_particles),  // (num_particles,3)
{   }  

// our static method to actually instantiate this singleton-esque object
Particles& Particles::GetInstance(const int num_particles) {
    static Particles* particles = new Particles(num_particles);
    return *particles;
}

void Particles::initLocation(GlobalSettings global_settings) {

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

    // generate random normal distribution for particle raddi
    std::default_random_engine generator;
    generator.seed(std::chrono::system_clock::now().time_since_epoch().count());
    std::normal_distribution<double> normal_distribution(mean_rad, stdev_rad);
    std::uniform_distribition<double> uniform_distribution(0.0, 1.0);   // between 0 and 1, we'll scale this later
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
        h_coordsn(i,2) = height - h_radius(i) - 1.0/3.0*height*uniform_distribution(generator);
    }

    // now check to see if particles overlap, if so we need to generate a new coordinate.
    // this can be expensive as it is an O(n!) check
    int i = num_particles_ - 2;
    int j; 
    double xd, yd, zd, dist;
    while (i > -1) {
        j = i+1;
        while (j < num_particles_) {
            xd = h_coords(i,0) - h_coords(j,0);
            yd = h_coords(i,1) - h_coords(j,1);
            zd = h_coords(i,2) - h_coords(j,2);
            dist = sqrt(pow(xd,2.) + pow(yd,2.) + pow(zd,2.));
            // re-generate center coord if the particles overlap and reset j index as we need to check all again
            if (dist < (h_radius(i)+h_radius(j)) {
                h_coordsn(i,0) = h_radius(i) + (length-2*h_radius(i))*uniform_distribution(generator);
                h_coordsn(i,1) = h_radius(i) + (width-2*h_radius(i))*uniform_distribution(generator);
                h_coordsn(i,2) = height - h_radius(i) - 1.0/3.0*height*uniform_distribution(generator);
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

    // use our first actual kokkos loops now to calculate the mass and volume and set initial velocity
    // the Kokkos::parallel_for loop construct will default to OpenMP threading if compiled CPU-only
    // or create a Cuda/HIP device function call (through lambda functions) if compiled with Cuda or HIP.
    // the lambda defaults to capture by value, which is what we want since Kokkos Views are treated
    // as pointers
    double rho = global_settings.rho_;
    Kokkos::parallel_for("calc_mass_vol", num_particles_, KOKKOS_LAMBDA(int i) {
        volume_(i) = 4./3.*M_PI*pow(radius_(i),3.);
        mass_(i) = rho*volume_(i);
        vn_(i,0) = 0.0;
        vn_(i,1) = 0.0;
        vn_(i,2) = 0.0;
    });


    /*
    // first task is to generate random numbers. We do this in parallel using
    // built-in Kokkos helper classes, for details see 
    // https://github.com/kokkos/kokkos/blob/master/example/tutorial/Algorithms/01_random_numbers/random_numbers.cpp
    static_cast<uint64_t> rand_seed = time(NULL)
    Kokkos::Random_XorShift64_Pool<> rand_pool64(rand_seed);

    Kokkos::View<uint64_t*> rand_nums("rand_nums",num_particles_*3);

    // the Kokkos::parallel_for loop construct will default to OpenMP threading if compiled CPU-only
    // or create a Cuda/HIP device function call (through the GenerateRandom functor) if compiled with Cuda or HIP
    Kokkos::parallel_for("generate_random_nums_radii",num_particles_,
                         utilities::GenerateRandom<Kokkos::Random_XorShift64_Pool<> >(
                         rand_nums, rand_pool64, 3));
    */

}



} // namespace amdem

