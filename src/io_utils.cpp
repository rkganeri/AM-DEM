#include <string>   // std::string
#include <fstream>  // std::flush
#include <iostream> // std::cerr, std::endl

#include "fmt/core.h"
#include "Kokkos_Core.hpp"
#include "lean_vtk.hpp"

#include "io_utils.hpp"

namespace amdem {

// simple helper function to print message on its own line
void print(const std::string& msg) {
    std::cout << msg << std::endl << std::flush;
}

void plotState(const std::unique_ptr<Particles>& particles, const double current_time,
               const int step_num) {

    // 1. get kokkos views on host
    Kokkos::View<double*>::HostMirror h_radius = Kokkos::create_mirror_view(particles->radius_);
    Kokkos::View<double**>::HostMirror h_coordsnp1 = Kokkos::create_mirror_view(particles->coordsnp1_);
    Kokkos::View<double**>::HostMirror h_vnp1 = Kokkos::create_mirror_view(particles->vnp1_);
    Kokkos::deep_copy(h_radius, particles->radius_);
    Kokkos::deep_copy(h_coordsnp1, particles->coordsnp1_);
    Kokkos::deep_copy(h_vnp1, particles->vnp1_);
    
    // 2. copy views into std::vector... this is annoying but to use the lean-vtk lib, you need to pass
    // in arrays in vector format... if this turns into a large time sink might be worth writing my own
    // vtk writer
    const int dim = 3;
    const int num_particles = particles->num_particles_;
    std::vector<double> radius_vec(num_particles);
    std::vector<double> coords_vec(num_particles*dim);
    std::vector<double> vel_vec(num_particles*dim);
    // TODO: do this using CPU-only OpenMP
    for (int i=0; i<particles->num_particles_; i++) {
        radius_vec[i] = h_radius(i);
        coords_vec[3*i+0] = h_coordsnp1(i,0);
        coords_vec[3*i+1] = h_coordsnp1(i,1);
        coords_vec[3*i+2] = h_coordsnp1(i,2);
        vel_vec[3*i+0] = h_vnp1(i,0);
        vel_vec[3*i+1] = h_vnp1(i,1);
        vel_vec[3*i+2] = h_vnp1(i,2);
    }

    // 3. use lean-vtk lib to write out plot state, following example here:
    // https://github.com/mmorse1217/lean-vtk/blob/master/tests/test_lean_vtk.cpp
    std::string step_string = std::to_string(step_num);
    int num_prepend_zeros = 9 - step_string.length();
    step_string.insert(0,num_prepend_zeros,'0');
    std::string filename = fmt::format("powder_dep_step{}.vtu",step_string);
    
    leanvtk::VTUWriter writer;
    writer.add_scalar_field("radius", radius_vec);
    writer.add_vector_field("position", coords_vec, dim);
    writer.add_vector_field("velocity", vel_vec, dim);
    writer.write_point_cloud(filename, dim, coords_vec);

}

} // namespace amdem


