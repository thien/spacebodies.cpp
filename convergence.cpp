// g++ -O3 --std=c++11 convergence.cpp -o convergence

#include <iostream>
#include <math.h>


const double yo = 4.0;

double t        = 0.0;
double y        =  yo;
double lambda   = -1.0;
double dt       = 0.256;

// explicit euler
int main(void){
  // simulate till t = 2.0
  while (t < 2.0){
    t = t + dt;
    y = (1.0 + lambda * dt) * y;

    // calculate error
    double error = y - yo * exp(lambda * t);

    std::cout << "t = " << t << ",y=" << y << ", e=" << error << std::endl;
    // change dt if necessary (h refinement)
  }
}