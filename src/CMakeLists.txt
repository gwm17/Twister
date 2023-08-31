add_executable(Twister)

target_include_directories(Twister PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/.
    ${HDF5_INCLUDE_DIRS},
    ${GSL_INCLUDE_DIRS}
)

target_sources(Twister PRIVATE
    main.cpp
    core/Cluster.h
    core/GuessReader.h
    core/NuclearMap.h
    core/ParticleID.h
)

target_link_libraries(Twister ${HDF5_CXX_LIBRARIES} GSL::gsl)

set_target_properties(Twister PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${TWR_BINARY_DIR})