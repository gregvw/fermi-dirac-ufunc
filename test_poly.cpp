#include "horner.hpp"
#include <vector>
#include <iostream>
#include <iomanip>

typedef long double Scalar;

int main( int argc, char *argv[] ) {

  const unsigned int m = 11;
  const unsigned int n = 2;

  Scalar a[n+1] = {3,2,1}; // p(x) = 3*x^2 + 2*x + 1

  std::vector<Scalar> x(m);
  std::vector<Scalar> p(m);

  for( int i=0; i<m; ++i ) {
    x[i] = static_cast<Scalar>(2*i)/(m-1);
    p[i] = horner(x[i],a,2);
    std::cout << std::setw(16) << x[i] << std::setw(16) << p[i] << std::endl;
  }
  return 0;

}
