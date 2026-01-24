cmake_minimum_required(VERSION 3.10)

project(axisem LANGUAGES C Fortran)

# Default options
option(USE_NETCDF "Enable NetCDF" ON)
option(USE_PAR_NETCDF "Use parallel NetCDF" OFF)
option(SERIAL "Compile without MPI" OFF)
option(INCLUDE_MPI "Include mpif.h" OFF)

# Compiler flags
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native -fexternal-blas")
set(CMAKE_C_FLAGS_RELEASE "-O3 -march=native")
set(CMAKE_Fortran_FLAGS_DEBUG "-g -traceback")
set(CMAKE_C_FLAGS_DEBUG "-g")

if(SERIAL)
    add_compile_definitions(serial)
    if(NOT WIN32)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pthread")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -pthread")
    endif()
else()
    find_package(MPI REQUIRED)
    if(MPI_FOUND)
        link_libraries(MPI::MPI_Fortran)
    endif()
endif()

if(INCLUDE_MPI)
    add_compile_definitions(include_mpi)
endif()

if(USE_NETCDF)
    if(USE_PAR_NETCDF)
        add_compile_definitions(enable_parallel_netcdf)
    endif()
endif()

# Find NetCDF
if(USE_NETCDF)
    find_program(NF_CONFIG NAMES nf-config)
    
    if(NF_CONFIG)
        message(STATUS "Found nf-config: ${NF_CONFIG}")
        execute_process(COMMAND ${NF_CONFIG} --fflags OUTPUT_VARIABLE NETCDF_FFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
        execute_process(COMMAND ${NF_CONFIG} --flibs OUTPUT_VARIABLE NETCDF_LDFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
        
        # Handle Includes
        string(REGEX MATCHALL "-I[^ ]+" NETCDF_INCLUDES_RAW "${NETCDF_FFLAGS}")
        foreach(INC ${NETCDF_INCLUDES_RAW})
             string(SUBSTRING "${INC}" 2 -1 DIR)
             include_directories("${DIR}")
        endforeach()
        
        # Use flags directly
        set(NETCDF_LIBRARIES ${NETCDF_LDFLAGS})
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${NETCDF_FFLAGS}")
        
    else()
        # Fallback: Check environment variables often set by modules
        if(DEFINED ENV{NETCDF_ROOT})
             set(NETCDF_PATH $ENV{NETCDF_ROOT})
        elseif(DEFINED ENV{NETCDF_DIR})
             set(NETCDF_PATH $ENV{NETCDF_DIR})
        elseif(DEFINED ENV{NETCDF_PATH})
             set(NETCDF_PATH $ENV{NETCDF_PATH})
        endif()

        if(NETCDF_PATH)
            message(STATUS "Using NetCDF path from environment: ${NETCDF_PATH}")
            include_directories(${NETCDF_PATH}/include)
            link_directories(${NETCDF_PATH}/lib)
            set(NETCDF_LIBRARIES netcdff)
        else()
             # If no explicit path, assume compiler can find it or it's in standard paths.
             # We rely on user environment variables (CPATH, LIBRARY_PATH) implicitly if we don't set anything
             # or we can try to find the library.
             find_library(NETCDF_LIB NAMES netcdff)
             if(NETCDF_LIB)
                message(STATUS "Found NetCDF library: ${NETCDF_LIB}")
                set(NETCDF_LIBRARIES ${NETCDF_LIB})
                find_path(NETCDF_INCLUDE_DIR netcdf.mod HINTS ENV CPATH ENV C_INCLUDE_PATH ENV INCLUDE)
                if(NETCDF_INCLUDE_DIR)
                    include_directories(${NETCDF_INCLUDE_DIR})
                endif()
             else()
                message(STATUS "NetCDF library not found via find_library, assuming -lnetcdff works")
                set(NETCDF_LIBRARIES netcdff)
             endif()
        endif()
    endif()
endif()

add_subdirectory(MESHER)
add_subdirectory(SOLVER)
