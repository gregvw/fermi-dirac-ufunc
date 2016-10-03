#ifndef HORNER_HPP
#define HORNER_HPP

#include <cmath>

/*  Evaluate a polynomial of the form
 *
 *  p(x) = a[0]*x^n + ... + a[n-1]*x + a[n]
 *
 *  using Horner's method with fused multiply-accumulate 
 *
 *  Supported types: float, double, long double
 */

template<typename T>
T horner(const T &x, const T * const a, unsigned int n) {

  T p = a[0];
  for( unsigned int i=1; i<=n; ++i ) {
    p = std::fma(p,x,a[i]);
  }
  return p;
} 

#endif // HORNER_HPP
