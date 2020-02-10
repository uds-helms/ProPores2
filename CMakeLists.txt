cmake_minimum_required(VERSION 3.1)

if(${CMAKE_VERSION} VERSION_LESS 3.15)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
else()
    cmake_policy(VERSION 3.15)
endif()

project(Propores VERSION 2.0 LANGUAGES CXX)

# TODO
add_executable(propores src/propores.cpp)

add_subdirectory(src)
add_subdirectory(tests)

target_link_libraries(propores propores_source)
target_compile_features(propores PUBLIC cxx_std_17)
target_compile_features(propores_source PUBLIC cxx_std_17)

install (TARGETS propores RUNTIME DESTINATION bin)