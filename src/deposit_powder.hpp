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

} // namespace amdem

#endif // AMDEM_DEPOSIT_POWDER_HPP




