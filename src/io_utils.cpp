#include <string>   // std::string
#include <fstream>  // std::flush
#include <iostream> // std::cerr, std::endl

#include "io_utils.hpp"

namespace amdem {

    // simple helper function to print message on its own line
    void print(const std::string& msg) {
        std::cout << msg << std::endl << std::flush;
    }

} // namespace amdem


