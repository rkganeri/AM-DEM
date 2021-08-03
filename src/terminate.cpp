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

        std::cout << "\n ERROR: " << msg << "\n" << std::endl;

        exit(EXIT_FAILURE);
    }

    // normal termination (called at end of main)
    void terminateNormal(const std::string& msg) {

        Kokkos::finalize();

        std::cout << "\n" << msg << "\n" << std::endl;
    }

} // namespace amdem

