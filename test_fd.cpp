#include "fermidirac.hpp"
#include <vector>
#include <iostream>
#include <iomanip>

typedef double Scalar;

int main( int argc, char *argv[] ) {

  int m = 11;

  FermiDirac<Scalar> fd;

  std::vector<Scalar> x(m);
  std::vector<Scalar> f(m);

  for(int i=0; i<m; ++i) {

    Scalar t = static_cast<Scalar>((m-1-i))/(m-1);
    x[i] = -6*t + 10*(1-t); 
    f[i] = fd.fermi_half(x[i]);
    std::cout << std::setw(16) << x[i] << std::setw(16) << f[i] << std::endl;

  }
   
  

  return 0;
}
