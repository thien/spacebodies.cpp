# create build folder if it doesn't exist.
if [ ! -d build ]; then
  mkdir -p build
fi
g++ -O3 --std=c++11 spacebodies.cpp -o build/spacebodies 
# ./build/spacebodies 100.0 2.0 1.0 3.0 4.0 1.0 2.0 0.2 \ 1.0 2.0 4.0 1.0 1.0 0.1 2.0 \ 1.0 0.1 0.1 0.1 0.1 1.0 2.0
# ./spacebodies 100.0   2.0 0.0 0.0 0.0 1.0 0.0 0.2 \ 0.0 0.0 0.0 0.0 0.0 0.0 2.0 \ 1.0 0.0 0.0 0.0 0.0 0.0 2.0