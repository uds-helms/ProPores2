DIR="build/"
if [ -d "$DIR" ]; then
  cd "$DIR"
  rmdir -v -rf delete
  cd ..
  mkdir "$DIR"
fi
cmake -S . -B build
cmake --build build --config Release --target install
chmod +x propores