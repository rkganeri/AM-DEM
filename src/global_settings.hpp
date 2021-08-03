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
        GlobalSettings(const int num_particles, const double mean_rad, const double stdev_rad);
        // delete copy operator
        GlobalSettings(GlobalSettings& ) = delete;
        GlobalSettings& operator=(GlobalSettings& ) = delete;
        // delete move operator
        GlobalSettings(GlobalSettings&& ) = delete;
        GlobalSettings& operator=(GlobalSettings&& ) = delete;
    
    public:
        // domain dimensions
        const double length_ = 0.5e-03; // m
        const double width_ = 1.0e-03;
        const double height_ = 3.0e-03;

        // particle size distribution
        const double mean_rad_;
        const double stdev_rad_;
        const double min_rad_ = 4.0e-06;    
        const double max_rad_ = 23.0e-06;
        const int num_particles_;

        // simulation controls
        const double t_start_ = 0.0; // s
        const double t_end_ = 0.1;
        const double dt_ = 5.0e-08;
        const double plot_time_freq_ = 1.0e-03;

        // material properties
        const double rho_ = 7952.0; // kg/m^3
        const double nu_ = 0.26;    
        const double youngs_mod_ = 193.0e09; // Pa   

        // declare constructor
        static GlobalSettings& getInstance(const int num_particles, const double mean_rad, 
                                           const double stdev_rad);
};

} // namespace amdem

#endif // AMDEM_GLOBAL_SETTINGS_HPP

