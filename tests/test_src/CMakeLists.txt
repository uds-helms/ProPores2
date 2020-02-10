include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR}, ${CMAKE_CURRENT_LIST_DIR}/test_files)
include_directories(${CMAKE_CURRENT_LIST_DIR}/../../src)
add_executable(propores_tests test_basics.cpp test_pore_ID.cpp test_grid.cpp test_vector.cpp test_atom.cpp test_reader.cpp test_angles.cpp test_enums.cpp)
set_property(TARGET propores_tests PROPERTY VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")

target_link_libraries(propores_tests gtest gtest_main)
target_link_libraries(propores_tests propores_source)
target_compile_features(propores_tests PUBLIC cxx_std_17)