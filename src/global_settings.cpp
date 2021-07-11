#include <string>   // std::string
#include <fstream>  // std::flush
#include <iostream> // std::cerr, std::endl

#include "global_settings.hpp"

namespace amdem {

// Here we set all the simulation parameters.  This is the only
// section that needs to be edited (and re-compiled) if changing settings.
// Note that the units are standard SI (m-kg-s)
GlobalSettings::GlobalSettings(const int num_particles, const int mean_rad) :
    length_(0.4e-03),       // m
    width_(0.8e-03),
    height_(3.5e-03),
    mean_rad_(mean_rad),
    min_rad_(mean_rad-0.1*mean_rad), // particle rad varies within +/- 10% of mean
    max_rad_(mean_rad+0.1*mean_rad),
    num_particles_(num_particles),
    t_start_(0.0),          // s
    t_end_(0.08),
    dt_(5.0e-08),
    rho_(7952.),            // kg/m^3
    youngs_mod_(193.e9),    // Pa
    nu_(0.26)
{ /* everything done via initialization */ }

} // namespace amdem

