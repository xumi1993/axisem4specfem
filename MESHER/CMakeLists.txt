
file(GLOB SRCS *.F90 *.f90)

if(USE_NETCDF)
    add_compile_definitions(unc)
endif()

# Output binary to project bin directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

add_executable(xmesh ${SRCS})

# Generate gitversion.h which is required for compilation
add_custom_command(
    OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/gitversion.h
    COMMAND sh -c "echo \"character(len=*), parameter :: gitversion = \\\"`git rev-parse HEAD`\\\"\" > ${CMAKE_CURRENT_SOURCE_DIR}/gitversion.h"
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    VERBATIM
)
# Ensure gitversion.h is generated before compiling sources that might include it
# Since F90 files often dont Include headers but modules, check if any files include 'gitversion.h'.
# Actually the Makefile generated it. 
# We need to make sure the custom command is executed.
add_custom_target(mesher_gitversion DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/gitversion.h)
add_dependencies(xmesh mesher_gitversion)

if(USE_NETCDF)
    target_link_libraries(xmesh ${NETCDF_LIBRARIES})
endif()


