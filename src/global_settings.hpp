#ifndef AMDEM_GLOBAL_SETTINGS_HPP
#define AMDEM_GLOBAL_SETTINGS_HPP

#include <string>   
#include <fstream>


namespace amdem {

// This is a bit hacky but it's easiest just to store these global settings in an object.
// We'll do the actual setting of params in the cpp file so if we change settings only that file must
// be re-compiled.  If actual input deck parsing is added down the line we may re-visit this.
class GlobalSettings {
    // N.B. while presumably only 1 instance of Global setting will ever exist, it seems unnecessary
    // to enforce this to be a singleton class so I won't do that at this time.
    public:
        // domain dimensions
        const double length_;
        const double width_;
        const double height_;

        // particle size distribution
        const double mean_rad_;
        const double stdev_rad_;
        const double min_rad_;
        const double max_rad_;
        const double num_particles_;

        // simulation controls
        const double t_start_;
        const double t_end_;
        const double dt_;

        // material properties
        const double rho_;
        const double nu_;
        const double youngs_mod_;

        // declare constructor
        GlobalSettings(const int num_particles, const int mean_rad, const int stdev_rad);
};

} // namespace amdem

#endif // AMDEM_GLOBAL_SETTINGS_HPP

