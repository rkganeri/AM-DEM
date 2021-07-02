
set(DEM-LIB_SOURCE_FILES
    ${DIABLO-MTI-LIB_SOURCE_FILES}
    ${CMAKE_CURRENT_LIST_DIR}/io_utils.cpp
    ${CMAKE_CURRENT_LIST_DIR}/terminate.cpp
    PARENT_SCOPE
   )

set(DEM-EXE_SOURCE_FILES
    ${CMAKE_CURRENT_LIST_DIR}/powder_deposition_main.cpp
    PARENT_SCOPE
   )

# Add header files
set(DEM-LIB_HEADER_FILES
    ${DEM-LIB_HEADER_FILES}
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
