#ifndef AMDEM_DEPOSIT_POWDER_HPP
#define AMDEM_DEPOSIT_POWDER_HPP

#include <string>  
#include <memory>

#include "global_settings.hpp"
#include "particles.hpp"
#include "bins.hpp"

namespace amdem {

bool depositPowder(std::unique_ptr<Particles>& particles,
                   const GlobalSettings& global_settings);

void updateRK4SubStep(const std::unique_ptr<Particles>& particles,
                      Kokkos::View<double**> y_pos,
                      Kokkos::View<double**> y_vel,
                      const Kokkos::View<double**> y_vel_prev,
                      const Kokkos::View<double**> psi,
                      const double dt, const int sub_step);
} // namespace amdem

#endif // AMDEM_DEPOSIT_POWDER_HPP




