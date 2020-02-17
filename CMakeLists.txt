# minimum version required for target_compile_features
cmake_minimum_required(VERSION 3.15)
# set project information
project(Propores VERSION 2.0 LANGUAGES CXX)
# add the source code and test sub-directories
add_subdirectory(src)
add_subdirectory(tests)
# create the propores executable
add_executable(propores src/propores.cpp)
# link the executable with the source code library
target_link_libraries(propores propores_source)
# make sure that the code gets compiled with C++17
target_compile_features(propores PUBLIC cxx_std_17)
# setup the installation
install (TARGETS propores RUNTIME DESTINATION ${CMAKE_SOURCE_DIR})
