# minimum version required for target_compile_features
cmake_minimum_required(VERSION 3.1)

# use the newest version between 3.1 and 3.15
if(${CMAKE_VERSION} VERSION_LESS 3.15)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
else()
    cmake_policy(VERSION 3.15)
endif()
# set project information
project(Propores VERSION 2.0 LANGUAGES CXX)
# set release type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
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
