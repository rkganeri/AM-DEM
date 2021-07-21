#ifndef AMDEM_IO_UTILS_HPP
#define AMDEM_IO_UTILS_HPP

#include <string>   
#include <fstream>

#include "particles.hpp"


namespace amdem {

void print(const std::string& msg);

void plotState(const std::unique_ptr<Particles>& particles, const double current_time,
               const int step_num);

} // namespace amdem

#endif // AMDEM_IO_UTILS_HPP

