DIR="build/"
if [ -d "$DIR" ]; then
  rm -rf "$DIR"
  mkdir "$DIR"
fi
cmake -S . -B "$DIR"
cmake --build "$DIR" --config Release --target install
chmod +x propores