cmake_minimum_required(VERSION 3.16)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED 20)

if(NOT CMAKE_BUILD_TYPE STREQUAL "Release")
    set(CMAKE_BUILD_TYPE "Debug")
    message("Building debug")
else()
    message("Building release")
endif()

project(Twister)

set(TWR_BINARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/bin)
set(CSV_CXX_STANDARD 20)

find_package(GSL REQUIRED)

set(HIGHFIVE_USE_BOOST Off)

add_subdirectory(vendor/nlohmann_json)
add_subdirectory(vendor/catima)
add_subdirectory(vendor/csv-parser)
add_subdirectory(vendor/HighFive)
add_subdirectory(src)