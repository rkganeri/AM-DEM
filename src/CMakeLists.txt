set(amdem-lib_source_files
    ${amdem-lib_source_files}
    ${CMAKE_CURRENT_LIST_DIR}/io_utils.cpp
    ${CMAKE_CURRENT_LIST_DIR}/terminate.cpp
    ${CMAKE_CURRENT_LIST_DIR}/global_settings.cpp
    ${CMAKE_CURRENT_LIST_DIR}/particles.cpp
    ${CMAKE_CURRENT_LIST_DIR}/bins.cpp
    PARENT_SCOPE
   )

set(amdem-exe_source_files
    ${CMAKE_CURRENT_LIST_DIR}/amdem_main.cpp
    PARENT_SCOPE
   )

# Add header files
set(amdem-lib_header_files
    ${amdem-lib_header_files}
    ${CMAKE_CURRENT_LIST_DIR}/io_utils.hpp
    ${CMAKE_CURRENT_LIST_DIR}/terminate.hpp
    ${CMAKE_CURRENT_LIST_DIR}/global_settings.hpp
    ${CMAKE_CURRENT_LIST_DIR}/particles.hpp
    ${CMAKE_CURRENT_LIST_DIR}/bins.hpp
    ${CMAKE_CURRENT_LIST_DIR}/utilities.hpp
    PARENT_SCOPE
   )

#install(
#    FILES
#        ${CMAKE_CURRENT_LIST_DIR}/io_utils.hpp
#        ${CMAKE_CURRENT_LIST_DIR}/terminate.hpp
#    DESTINATION
#        ${PROJECT_SOURCE_DIR}/src/include
#   )
