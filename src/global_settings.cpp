#include <string>   // std::string
#include <fstream>  // std::flush
#include <iostream> // std::cerr, std::endl

#include "global_settings.hpp"

namespace amdem {

// Here we set all the simulation parameters.  This is the only
// section that needs to be edited (and re-compiled) if changing settings.
// Note that the units are standard SI (m-kg-s)
// N.B. These settings are consistent with the original Fortran code PowderDepositionKhairallah.F90
// which was used in the paper https://doi.org/10.1007/s10035-016-0626-0
GlobalSettings::GlobalSettings(const int num_particles, const double mean_rad, const double stdev_rad) 
    : mean_rad_(mean_rad),
      stdev_rad_(stdev_rad),
      num_particles_(num_particles)
{ /* everything else set in hpp file */ }

// our static method to actually instantiate this singleton object
GlobalSettings& GlobalSettings::getInstance(const int num_particles, const double mean_rad, 
                                            const double stdev_rad) {
    static GlobalSettings* global_settings = new GlobalSettings(num_particles, mean_rad, stdev_rad);
    return *global_settings;
}


} // namespace amdem

