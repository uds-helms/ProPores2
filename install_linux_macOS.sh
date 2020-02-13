DIR="build/"
if [ -d "$DIR" ]; then
  rm -rf "$DIR"
  mkdir "$DIR"
fi
cmake -S . -B build
cmake --build build --config Release --target install
chmod +x propores