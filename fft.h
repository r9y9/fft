#pragma once

#include <iterator>
#include <vector>
#include <cmath>
#include <stdexcept>

namespace sp {

namespace fft_detail {

template <class Iterator, class Complex>
inline void fft_inner(Iterator a, Iterator b, int flip)
{
  typedef typename Complex::value_type value_type;
  //const Complex I(0, 1);

  const ssize_t n = b - a;
  
  const value_type pi = atan(1.0) * 4;
  value_type theta = 2.0 * pi / n * flip;
  
  for (ssize_t m = n; m >= 2; m >>= 1) {
    int mh = m >> 1;
    for (ssize_t i = 0; i < mh; ++i) {
      Complex w(cos(i*theta), sin(i*theta));
      //Complex w = exp(i*theta*I); // slower
      for (ssize_t j = i; j < n; j += m) {
        ssize_t k = j + mh;
        Complex x = a[j] - a[k];
        a[j] += a[k];
        a[k] = w * x;
      }
    }
    theta *= 2;
  }
  for (ssize_t i = 0, j = 1; j < n - 1; j++) {
    for (ssize_t k = n >> 1; k > (i ^= k); k >>= 1);
    if (j < i) swap(a[i], a[j]);
  }

  if (flip == -1) {
    for (ssize_t i = 0; i < n; ++i) a[i] /= n;
  }
}

inline bool is_pow2(ssize_t s)
{
  if (s <= 0) return false;
  return (s & (s - 1)) == 0;
}

template <class Iterator>
inline void fft_wrap(Iterator b, Iterator e, int flip)
{
  ssize_t length = e - b;
  if (!is_pow2(length)) {
    throw std::length_error("fft: length must be power of 2");
  }
  typedef typename std::iterator_traits<Iterator>::value_type value_type;
  fft_inner<Iterator, value_type>(b, e, flip);
}

} // fft_detail

template <class Iterator>
inline void fft1d(Iterator b, Iterator e)
{
  fft_detail::fft_wrap(b, e, 1);
}

template <class Iterator>
inline void ifft1d(Iterator b, Iterator e)
{
  fft_detail::fft_wrap(b, e, -1);
}

} // sp
