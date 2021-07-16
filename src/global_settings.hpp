#ifndef AMDEM_GLOBAL_SETTINGS_HPP
#define AMDEM_GLOBAL_SETTINGS_HPP

#include <string>   
#include <fstream>


namespace amdem {

// This is a bit hacky but it's easiest just to store these global settings in an object.
// We'll do the actual setting of params in the cpp file so if we change settings only that file must
// be re-compiled.  If actual input deck parsing is added down the line we may re-visit this.
class GlobalSettings {
    // N.B. we make this a singleton class (for now) as it stores most of the global settings used
    // throughout the calculations and we don't want to accidentally allow a copy or alteration of it
    // (singleton probably unnecessary, but fun to try out ... )
    private:
        // make the constructor private and delete the default one
        GlobalSettings() = delete;
        GlobalSettings(const int num_particles, const int mean_rad, const int stdev_rad);
        // delete copy operator
        GlobalSettings(GlobalSettings& ) = delete;
        GlobalSettings& operator=(GlobalSettings& ) = delete;
        // delete move operator
        GlobalSettings(GlobalSettings&& ) = delete;
        GlobalSettings& operator=(GlobalSettings&& ) = delete;
    
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
        static GlobalSettings& getInstance(const int num_particles, const int mean_rad, 
                                           const int stdev_rad);
};

} // namespace amdem

#endif // AMDEM_GLOBAL_SETTINGS_HPP

