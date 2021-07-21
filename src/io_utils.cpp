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
    // 2. copy views into std::vector

    // 3. use lean-vtk lib to write out plot state, following example here:
    // https://github.com/mmorse1217/lean-vtk/blob/master/tests/test_lean_vtk.cpp
    std::string step_string = std::to_string(step_num);
    int num_prepend_zeros = 9 - step_string.length();
    step_string.insert(0,num_prepend_zeros,'0');
    std::string filename = fmt::format("powder_dep_step{}.vtu",step_string);
    leanvtk::VTUWriter writer;
    //writer.add_scalar_field("temperature", field_vec);
    //writer.write_volume_mesh(filename, dim, num_elem_nodes, node_coords_vec, elem_conn_vec);

}

} // namespace amdem


