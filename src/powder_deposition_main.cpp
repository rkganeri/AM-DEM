// Insert Copyright stuff here
)
#include <memory>
#include <string>

#include "CLI11/CLI11.hpp"
#include "fmt/fmt.hpp"
#include "Kokkos_Core.hpp"

#include "io_utils.hpp"
#include "terminate.hpp"

int main(int argc, char *argv[]) {

    // parse command line
    CLI::App app{"Powder_deposition_DEM"};
    
    double x_range = 0.0;
    app.add_option("-x, --xrange", x_range, "X-range of bounding box")
        ->required()->check(CLI::PositiveNumber);

    double y_range = 0.0;
    app.add_option("-y, --yrange", y_range, "Y-range of bounding box")
        ->required()->check(CLI::PositiveNumber);

    double z_range = 0.0;
    app.add_option("-z, --zrange", z_range, "Z-range of bounding box")
        ->required()->check(CLI::PositiveNumber);

    double radius = 0.0;
    app.add_option("-r, --radius", radius, "Radius of DEM particles")
        ->required()->check(CLI::PositiveNumber);

    int num_particles = 0;
    app.add_option("-n, --numparticles", radius, "Number of DEM particles")
        ->required()->check(CLI::PositiveNumber);

    // TODO: set up input deck parsing (maybe yaml input decks), but not too many parameters
    // so easy enough to just pass through as arguments for now

    // macro defined in CLI11.hpp, see CLI11 Readme docs
    CLI11_PARSE(app, argc, argv);


    // initialize Kokkos
    Kokkos::initialize(argc, argv);





    powderDep::terminateNormal();


}
