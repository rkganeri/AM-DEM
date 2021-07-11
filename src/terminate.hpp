#ifndef AMDEM_TERMINATE_HPP
#define AMDEM_TERMINATE_HPP

namespace amdem {

    void terminateError(const std::string& msg="quitting...");

    void terminateNormal(const std::string& msg="Program finished, normal termination.");

} // namespace amdem

#endif // AMDEM_TERMINATE_HPP

