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
    Core/Cluster.h
    Core/Cluster.cpp
    Core/Guess.h
    Core/GuessReader.h
    Core/GuessReader.cpp
    Core/NuclearMap.h
    Core/NuclearMap.cpp
    Core/ParticleID.h
    Core/ParticleID.cpp
    Core/HDFReader.h
    Core/HDFReader.cpp
    Core/Target.h
    Core/Target.cpp
    Core/Precision.h
    Core/Solver.h
    Core/Solver.cpp
    Core/Constants.h
    Utils/Timer.h
)

target_link_libraries(Twister ${HDF5_CXX_LIBRARIES} GSL::gsl nlohmann_json::nlohmann_json catima::catima csv HighFive)

set_target_properties(Twister PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${TWR_BINARY_DIR})