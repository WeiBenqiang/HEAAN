
cmake_minimum_required(VERSION 3.0)

project(TestMAX)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -Wall")

set(CMAKE_BUILD_TYPE Debug)

link_directories (${CMAKE_CURRENT_SOURCE_DIR}/../lib/) 

add_executable(ex2 example2.cpp)

target_link_libraries(ex2 HEAAN ntl pthread )