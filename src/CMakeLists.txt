add_executable(Twister)

target_include_directories(Twister PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/.
    ${GSL_INCLUDE_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}/../vendor/nlohmann_json/include
    ${CMAKE_CURRENT_SOURCE_DIR}/../vendor/catima/include
    ${CMAKE_CURRENT_SOURCE_DIR}/../vendor/csv-parser/single_include
    ${CMAKE_CURRENT_SOURCE}/../vendor/HighFive/include
)

target_sources(Twister PRIVATE
    main.cpp
    core/Cluster.h
    core/Cluster.cpp
    core/Guess.h
    core/GuessReader.h
    core/GuessReader.cpp
    core/NuclearMap.h
    core/NuclearMap.cpp
    core/ParticleID.h
    core/ParticleID.cpp
    core/HDFReader.h
    core/HDFReader.cpp
    core/Target.h
    core/Target.cpp
    core/Precision.h
    core/Solver.h
    core/Solver.cpp
    core/Constants.h
)

target_link_libraries(Twister ${HDF5_CXX_LIBRARIES} GSL::gsl nlohmann_json::nlohmann_json catima::catima csv HighFive)

set_target_properties(Twister PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${TWR_BINARY_DIR})