#include "fft.h"

#include <iostream>
#include <complex>

using namespace std;

int main(int argc, char** argv) 
{
  const int N = 512;
  vector<complex<double> > array(N);
  
  for (int i = 0; i < N; ++i) {
    array[i] = i;
  }

  // FFT
  sp::fft1d(array.begin(), array.end());

  // do something

  // IFFT
  sp::ifft1d(array.begin(), array.end());

  return 0;
}
