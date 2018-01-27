# create build folder if it doesn't exist.
if [ ! -d build ]; then
  mkdir -p build
fi
# compile spacebodies code
g++ -O3 --std=c++11 spacebodies.cpp -o build/spacebodies

# run spacebodies; note that it comes in the form of px py pz vx vy vz m

# general spacebodies
# ./build/spacebodies 100.0 2.0 1.0 3.0 4.0 1.0 2.0 0.2 \ 1.0 2.0 4.0 1.0 1.0 0.1 2.0 \ 1.0 0.1 0.1 0.1 0.1 1.0 2.0

# collisions
./build/spacebodies 2.2  0.0000001 0.0 0.0 -111.1 0.0 0.0 5.0 \ -0.0000001 0.0 0.0 111.1 0.01 0.0 5.0
