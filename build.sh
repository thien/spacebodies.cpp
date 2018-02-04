# create build folder if it doesn't exist.
if [ ! -d build ]; then
  mkdir -p build
fi
# compile spacebodies code
g++ -O3 --std=c++11 spacebodies.cpp -o build/spacebodies

# run spacebodies; note that it comes in the form of px py pz vx vy vz m

# general spacebodies
./build/spacebodies 1.0 2.0 1.0 3.0 4.0 1.0 2.0 0.2 \ 1.0 2.0 4.0 1.0 1.0 0.1 2.0 \ 1.0 0.1 0.1 0.1 0.1 1.0 2.0

# collisions
# ./build/spacebodies 0.002  \ 0.00001 0.0 0.0 -0.10001 0.0 0.0 5.0 \ -0.00001 0.0 0.0 0.1 0.0 0.0 5.0

# collisions
# ./build/spacebodies 4.05 \ 10 0.0 0.0 -2.5 0.0 0.0 5.0 \ -10 0.0 0.0 +2.5 0.0 0.0 5.0
