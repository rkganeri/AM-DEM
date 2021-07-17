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
    : length_(0.5e-03),       // m
      width_(1.0e-03),
      height_(3.0e-03),
      mean_rad_(mean_rad),
      stdev_rad_(stdev_rad),
      min_rad_(4.0e-06), 
      max_rad_(23.0e-06),
      num_particles_(num_particles),
      t_start_(0.0),          // s
      t_end_(0.05),
      dt_(2.0e-08),
      rho_(7952.),            // kg/m^3
      nu_(0.26),
      youngs_mod_(193.e9)    // Pa
{ /* everything done via initialization */ }

// our static method to actually instantiate this singleton object
GlobalSettings& GlobalSettings::getInstance(const int num_particles, const double mean_rad, 
                                            const double stdev_rad) {
    static GlobalSettings* global_settings = new GlobalSettings(num_particles, mean_rad, stdev_rad);
    return *global_settings;
}


} // namespace amdem

