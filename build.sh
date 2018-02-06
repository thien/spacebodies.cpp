# create build folder if it doesn't exist.
if [ ! -d build ]; then
  mkdir -p build
fi
# compile spacebodies code
g++ -O3 --std=c++11 spacebodies.cpp -o build/spacebodies





# ./build/spacebodies 

# run spacebodies; note that it comes in the form of px py pz vx vy vz m

# general spacebodies
# ./build/spacebodies 1.0 2.0 1.0 3.0 4.0 1.0 2.0 0.2 \ 1.0 2.0 4.0 1.0 1.0 0.1 2.0 \ 1.0 0.1 0.1 0.1 0.1 1.0 2.0

# collisions
# ./build/spacebodies 0.002  \ 0.00001 0.0 0.0 -0.10001 0.0 0.0 5.0 \ -0.00001 0.0 0.0 0.1 0.0 0.0 5.0

# collisions (floating point positions)
# ./build/spacebodies 0.4201101 \ 0.1 0.1 0.1 -2 -2 -2 0.00000000001 \ -1 -1 -1 +2 +2 +2 0.00000000001
./build/spacebodies 0.1001 \ 0.2 0.2 0.2 -1.0 -1.0 -1.0 0.0000000000001 \ 0 0 0 +1.0 +1.0 +1.0 0.0000000000001

# collisions (whole)
# ./build/spacebodies 0.401 \ 1 0.0 0.0 -2.5 0.0 0.0 100.0 \ -1 0.0 0.0 +2.5 0.0 0.0 100.0

