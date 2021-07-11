#include <string>   // std::string
#include <fstream>  // std::flush
#include <iostream> // std::cerr, std::endl
#include <stdlib.h>     /* exit, EXIT_FAILURE */

#include "Kokkos_Core.hpp"      

#include "terminate.hpp"

namespace amdem {

    // function to be called if there is an error
    void terminateError(const std::string& msg) {

        Kokkos::finalize();

        std::cerr << std::endl << std::flush;
        std::cerr << "ERROR: " << msg << std::endl << std::flush;
        std::cerr << std::endl << std::flush;

        exit(EXIT_FAILURE);
    }

    // normal termination (called at end of main)
    void terminateNormal(const std::string& msg) {

        Kokkos::finalize();

        std::cout << std::endl << std::flush;
        std::cout << msg << std::endl << std::flush;
        std::cout << std::endl << std::flush;
    }

} // namespace amdem

