if(USE_NETCDF)
    add_compile_definitions(unc)
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

add_executable(xpost_processing post_processing.F90 nc_postroutines.F90)
add_executable(xfield_transform field_transform.F90)

if(USE_NETCDF)
    target_link_libraries(xpost_processing ${NETCDF_LIBRARIES})
    target_link_libraries(xfield_transform ${NETCDF_LIBRARIES})
endif()


