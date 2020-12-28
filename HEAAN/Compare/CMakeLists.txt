link_directories (${CMAKE_CURRENT_SOURCE_DIR}/../lib/) 
add_executable(Comp Comp.cpp)
target_link_libraries(Comp HEAAN ntl pthread )