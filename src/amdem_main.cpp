// Insert Copyright stuff here

#include <memory>
#include <string>

#include "CLI11.hpp"
#include "Kokkos_Core.hpp"
#include "fmt/core.h"

#include "io_utils.hpp"
#include "terminate.hpp"
#include "global_settings.hpp"
#include "particles.hpp"
#include "bins.hpp"
#include "deposit_powder.hpp"

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

    // TODO: set up input deck parsing (maybe yaml input decks) and re-configure global settings class

    // macro defined in CLI11.hpp, see CLI11 Readme docs
    CLI11_PARSE(app, argc, argv);

    if (verbose) {
        amdem::print("Echoing command line inputs: ");
        amdem::print(fmt::format("   mean radius = {} ",mean_rad));
        amdem::print(fmt::format("   std. deviation of radius = {} ",stdev_rad));
        amdem::print(fmt::format("   num_particles = {} ",num_particles));
    }

    // initialize Kokkos
    Kokkos::initialize(argc, argv);
    
    // we start initializing Kokkos views in the various routines below and they need to be destroyed 
    // (handled automaticallly once out of scope) before we terminate so wrapping everything below in braces
    {

    // create our global settings singleton object (does it need to be a singleton? prob not, but fun to try it out)
    amdem::GlobalSettings& global_settings = amdem::GlobalSettings::getInstance(num_particles,mean_rad,stdev_rad);

    // wrap the particles object in a unique pointer to ensure we don't accidentally
    // make copies, this object contains most of the data we use in the DEM calcs
    auto particles = std::make_unique<amdem::Particles>(num_particles, global_settings);

    particles->initParticles(global_settings);

    // main function where it all happens
    bool deposit_successful = amdem::depositPowder(particles, global_settings);
    if (not deposit_successful) amdem::terminateError("Failure during powder deposition");


    } // end wrapper brace to destroy kokkos views


    amdem::terminateNormal();


}
