// Insert Copyright stuff here

#include <memory>
#include <string>

#include "CLI11.hpp"
#include "Kokkos_Core.hpp"
#include "fmt/core.h"

#include "io_utils.hpp"
#include "terminate.hpp"
#include "global_settings.hpp"

int main(int argc, char *argv[]) {

    // parse command line
    CLI::App app{"AMDEM"};
    
    double mean_rad = 13.5e-06;
    app.add_option("-r, --radius", mean_rad, "Mean radius of DEM particles")
        ->check(CLI::PositiveNumber);

    double stdev_rad = 4.0e-06;
    app.add_option("-s, --stdev", stdev_rad, "Std. dev of radius of DEM particles")
        ->check(CLI::PositiveNumber);

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
        amdem::printMessage(fmt::format("   mean radius = {} ",mean_rad));
        amdem::printMessage(fmt::format("   std. deviation of radius = {} ",stdev_rad));
        amdem::printMessage(fmt::format("   num_particles = {} ",num_particles));
    }

    // initialize Kokkos
    Kokkos::initialize(argc, argv);

    // create our global settings object
    amdem::GlobalSettings global_settings(num_particles,mean_rad,stdev_rad);
    





    amdem::terminateNormal();


}
