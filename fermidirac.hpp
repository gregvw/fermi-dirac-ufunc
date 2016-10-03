#ifndef FERMI_DIRAC_HPP
#define FERMI_DIRAC_HPP

#include "horner.hpp"

/*  Approximate evaluation of Fermi-Diract integrals j=-1/2 and j=1/2
 *
 *  Reference: F.G. Lether, "Analytical Expansion and Numerical Approximation
 *  of the Fermi-Dirac Integrals Fj(x) of Order j=-1/2 and j=1/2", Journal of
 *  Scientific Computing, Dec 2000, Vol 15, No 4, pp/ 479--497
 */


template<typename T>
class FermiDirac {

private:

  T pi;
  T pi2;
  T sqrt2pi; 
  T* a;
  T* b;

  static const long numer_a[7]; // Expansion coefficients of j=1/2 approximation
  static const long denom_a[7];
  static const long numer_b[6]; // Expansion coefficients of j=-1/2 approximation
  static const long denom_b[6]; 

  T beta( const T &x, unsigned int k ) {
    T y = pi*(2*k-1);
    return std::sqrt(x*x + y*y);
  }

  // Terms which appear in the j=1/2 series
  T function_pos( const T &x, unsigned int k ) {
    return std::sqrt(function_f(x,k)-x);
  }

  // Terms which appear in the j=-1/2 series
  T function_neg( const T &x, unsigned int k ) {
    return -0.5*function_pos(x,k)/function_f(x,k);
  }


public:

  FermiDirac() : pi(4*std::atan(1)), pi2(pi*pi), sqrt2pi(std::sqrt(2*pi)) {
    a = new T[7];
    b = new T[6];
    for( int k=0; k<7; ++k ){
      a[k] = static_cast<T>(numer_a[k])/static_cast<T>(denom_a[k]);
    }
    for( int k=0; k<6; ++k ) {
      b[k] = static_cast<T>(numer_b[k])/static_cast<T>(denom_b[k]);
    }

  }

  virtual ~FermiDirac() {
    delete[] a;
    delete[] b;
  }

  // Approximate Fermi-Dirac integral F_{1/2}
  T fermi_half( const T &x ) {
    T sum = 0;
    T poly = horner(x,a,6);
    for( int k=1; k<=20; ++k ) {
      sum += std::sqrt(beta(x,k)-x);
    }

    return poly + 2*sqrt2pi*sum;
  }

  T fermi_zero( const T &x ) {
    return std::log(1+std::exp(x));
  }

  // Approximate Fermi-Dirac integral F_{-1/2}
  T fermi_minus_half( const T &x ) {
    T sum = 0;
    T poly = horner(x,b,5);
    for( int k=1; k<=20; ++k ) {
      T betak = beta(x,k);
      sum += std::sqrt(betak-x)/betak;
    } 
    return poly - sqrt2pi*sum;
  }

}; // class FermiDirac

template<typename T>
const long FermiDirac<T>::numer_a[7] = {1,-1,-13,85,3923,74141,-5990294};

template<typename T>
const long FermiDirac<T>::denom_a[7] = {770751818298,3574503105,184757992,3603084,220484,8289,7995};

template<typename T>
const long FermiDirac<T>::numer_b[6] = {-1,-1,-1,27,3923,8220};

template<typename T>
const long FermiDirac<T>::denom_b[6] = {128458636383,714900621,3553038,381503,110242,919};




#endif // FERMI_DIRAC_HPP
