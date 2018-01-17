//
// Explicit Euler for Numerical Algorithms course
//
// g++ Euler.cpp -o Euler
//
// (C) Tobias Weinzierl
//
#include <iostream>


double F(double f, double t) {
 return 0.0;
}


int main() {
  double t      =   0.0;
  double f      =   0.0;
  double dt     = 0.001;

  std::cout << "t, f" << std::endl;
  while (t<2.0) {                          // simulate until t=2.0
    t = t + dt;
    f = f + dt * F(f,t);
    std::cout << t << ", " << f << std::endl;
  }

  return 0;
}
