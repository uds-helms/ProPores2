set(CMAKE_CXX_STANDARD 17)
set(SOURCE_FILES pore_id.cpp basics.h vector.h grid.h atom.h reader.h atom.cpp basics.cpp pore_id.h reader.cpp
        trace.cpp trace.h settings.h gate_open.cpp gate_open.h gate.h angles.cpp angles.h enums.cpp enums.h)

add_library(propores_source STATIC ${SOURCE_FILES})