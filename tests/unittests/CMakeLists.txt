# include the source code directory that is supposed to be tested
include_directories(${CMAKE_CURRENT_LIST_DIR}/../../src)
# include the google test and test file directories
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR} ${CMAKE_CURRENT_LIST_DIR}/test_files)
# create the test executable and set the debugger property for visual studio
add_executable(propores_tests test_basics.cpp test_pore_ID.cpp test_grid.cpp test_vector.cpp test_atom.cpp test_reader.cpp test_angles.cpp test_enums.cpp)
set_property(TARGET propores_tests PROPERTY VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")
# link the test executable with google test and the propores source library
target_link_libraries(propores_tests gtest gtest_main)
target_link_libraries(propores_tests propores_source)
# make sure that the code gets compiled with C++17
target_compile_features(propores_tests PUBLIC cxx_std_17)
# setup the installation
install (TARGETS propores_tests RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR})