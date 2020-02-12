if exist "build\" rmdir /Q /S build
cmake -S . -B build
cmake --build build --config Release --target install
