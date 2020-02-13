DIR="build/"
if [ -d "$DIR" ]; then
  rmdir -v -rf delete "$DIR"
  mkdir "$DIR"
fi
cmake -S . -B build
cmake --build build --config Release --target install
chmod +x propores