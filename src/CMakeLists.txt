
set(amdem-lib_source_files
    ${amdem-lib_source_files}
    ${CMAKE_CURRENT_LIST_DIR}/io_utils.cpp
    ${CMAKE_CURRENT_LIST_DIR}/terminate.cpp
    PARENT_SCOPE
   )

set(amdem-exe_source_files
    ${CMAKE_CURRENT_LIST_DIR}/powder_deposition_main.cpp
    PARENT_SCOPE
   )

# Add header files
set(amdem-lib_header_files
    ${amdem-lib_header_files}
    ${CMAKE_CURRENT_LIST_DIR}/io_utils.hpp
    ${CMAKE_CURRENT_LIST_DIR}/terminate.hpp
    PARENT_SCOPE
   )

#install(
#    FILES
#        ${CMAKE_CURRENT_LIST_DIR}/io_utils.hpp
#        ${CMAKE_CURRENT_LIST_DIR}/terminate.hpp
#    DESTINATION
#        ${PROJECT_SOURCE_DIR}/src/include
#   )
