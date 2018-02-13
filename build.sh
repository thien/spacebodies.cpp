# create build folder if it doesn't exist.
if [ ! -d build ]; then
  mkdir -p build
fi

if [ ! -d paraview ]; then
  mkdir -p paraview
fi
# compile spacebodies code (Parallel)
# g++ -O3 --std=c++11 -fopenmp spacebodies.cpp -o build/spacebodies

# compile spacebodies code (Series)
g++ -O3 --std=c++11 spacebodies.cpp -o build/spacebodies 


# Run random bodies
./build/spacebodies 

# run spacebodies; note that it comes in the form of px py pz vx vy vz m

# general spacebodies
# ./build/spacebodies 1.0 2.0 1.0 3.0 4.0 1.0 2.0 0.2 \ 1.0 2.0 4.0 1.0 1.0 0.1 2.0 \ 1.0 0.1 0.1 0.1 0.1 1.0 2.0

# collisions
# ./build/spacebodies 0.002  \ 0.00001 0.0 0.0 -0.10001 0.0 0.0 5.0 \ -0.00001 0.0 0.0 0.1 0.0 0.0 5.0

# collisions (floating point positions)
# ./build/spacebodies 0.27502 \ 0.1 0.1 0.1 -2 -2 -2 0.00000000001 \ -1 -1 -1 +2 +2 +2 0.00000000001
# ./build/spacebodies 0.15 \ 0.1 0.1 0.1 +1.0 +1.0 +1.0 0.0000000000001 \ 0.21 0.21 0.21 -1.0 -1.0 -1.0 0.0000000000001 

# collisions (whole)
# ./build/spacebodies 1.0 \ 1 0.0 0.0 -10 0.0 0.0 100.0 \ -1 0.0 0.0 +10 0.0 0.0 100.0

