file(GLOB SRCS *.F90 *.f90 *.c)

add_compile_definitions(solver)
if(USE_NETCDF)
    add_compile_definitions(enable_netcdf)
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

add_executable(axisem ${SRCS})

if(NOT SERIAL)
    target_link_libraries(axisem MPI::MPI_Fortran)
endif()

add_custom_command(
    OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/gitversion.h
    COMMAND sh -c "echo \"character(len=*), parameter :: gitversion = \\\"`git rev-parse HEAD`\\\"\" > ${CMAKE_CURRENT_SOURCE_DIR}/gitversion.h"
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    VERBATIM
)
add_custom_target(solver_gitversion DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/gitversion.h)
add_dependencies(axisem solver_gitversion)

if(USE_NETCDF)
    target_link_libraries(axisem ${NETCDF_LIBRARIES})
endif()

add_subdirectory("UTILS")


