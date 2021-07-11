// Insert Copyright stuff here

#include <memory>
#include <string>

#include "CLI11.hpp"
#include "Kokkos_Core.hpp"
#include "fmt/core.h"

#include "io_utils.hpp"
#include "terminate.hpp"

int main(int argc, char *argv[]) {

    // parse command line
    CLI::App app{"AMDEM"};
    
    double radius = 0.0;
    app.add_option("-r, --radius", radius, "Radius of DEM particles")
        ->required()->check(CLI::PositiveNumber);

    int num_particles = 0;
    app.add_option("-n, --numparticles", num_particles, "Number of DEM particles")
        ->required()->check(CLI::PositiveNumber);

    bool verbose = false;
    app.add_flag("-v, --verbose", verbose, "Enable verbose mode");

    // TODO: set up input deck parsing (maybe yaml input decks), but not too many parameters
    // so easy enough to just pass through as arguments for now

    // macro defined in CLI11.hpp, see CLI11 Readme docs
    CLI11_PARSE(app, argc, argv);

    if (verbose) {
        amdem::printMessage("Echoing command line inputs: ");
        amdem::printMessage(fmt::format("   x_range = {} ",x_range));
        amdem::printMessage(fmt::format("   y_range = {} ",y_range));
        amdem::printMessage(fmt::format("   z_range = {} ",z_range));
        amdem::printMessage(fmt::format("   radius = {} ",radius));
        amdem::printMessage(fmt::format("   num_particles = {} ",num_particles));
    }

    // initialize Kokkos
    Kokkos::initialize(argc, argv);

    // create our global settings object
    amdem::GlobalSettings global_settings(num_particles,radius);
    





    amdem::terminateNormal();


}
