DIR="build/"
if [ -d "$DIR" ]; then
  rm -rf "$DIR"
fi
mkdir "$DIR"
cmake -S . -B "$DIR"
cmake --build "$DIR" --config Release --target install
chmod +x propores