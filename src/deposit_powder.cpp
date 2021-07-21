#include <string>  
#include <memory>
#include <cmath>

#include "Kokkos_Core.hpp"
#include "fmt/core.h"

#include "deposit_powder.hpp"
#include "io_utils.hpp"
#include "terminate.hpp"

namespace amdem {

bool depositPowder(std::unique_ptr<Particles>& particles, std::unique_ptr<Bins>& bins,
                   const GlobalSettings& global_settings) {

    // we are using an explicit RK-4 time stepping scheme, so our step size is fixed
    double current_time = global_settings.t_start_;
    const double end_time = global_settings.t_end_;
    const double dt = global_settings.dt_;
    const int num_time_steps = static_cast<int>(ceil((end_time-current_time) / dt));

    // create initial plot state
    plotState(particles, current_time, 0);

    for (int istep=1; istep<=num_time_steps; istep++) {



    }

    return true;

}

} // namespace amdem



