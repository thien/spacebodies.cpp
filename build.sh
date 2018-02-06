# create build folder if it doesn't exist.
if [ ! -d build ]; then
  mkdir -p build
fi
# compile spacebodies code
g++ -O3 --std=c++11 spacebodies.cpp -o build/spacebodies


# Run random bodies
./build/spacebodies 

# run spacebodies; note that it comes in the form of px py pz vx vy vz m

# general spacebodies
# ./build/spacebodies 1.0 2.0 1.0 3.0 4.0 1.0 2.0 0.2 \ 1.0 2.0 4.0 1.0 1.0 0.1 2.0 \ 1.0 0.1 0.1 0.1 0.1 1.0 2.0

# collisions
# ./build/spacebodies 0.002  \ 0.00001 0.0 0.0 -0.10001 0.0 0.0 5.0 \ -0.00001 0.0 0.0 0.1 0.0 0.0 5.0

# collisions (floating point positions)
# ./build/spacebodies 0.4201101 \ 0.523 0.523 0.523 -2.5 -2.5 -2.5 100.0 \ -1.061 -1.061 -1.061 +2.5 +2.5 +2.5 100.2
# ./build/spacebodies 0.4 \ 0.22300000000001 0.22300000000001 0.22300000000001 -2.500000000001 -2.500000000001 -2.500000000001 0.1 \ -0.06100000000001 -0.06100000000001 -0.06100000000001 +2.500000000001 +2.500000000001 +2.500000000001 0.1
# ./build/spacebodies 1 0 0 0 0 0 0 0.0000000000000001 0 0 -0.000001 0 0 0 0.0000000000000001
# collisions (whole)
# ./build/spacebodies 0.401 \ 1 0.0 0.0 -2.5 0.0 0.0 100.0 \ -1 0.0 0.0 +2.5 0.0 0.0 100.0