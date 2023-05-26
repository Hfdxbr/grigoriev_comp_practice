#pragma once

#include <array>
#include <cmath>
#include <numeric>
#include <ostream>
#include <type_traits>

template <class T, size_t N>
class FixedVector {
 public:
  template <class... Us>
  FixedVector(Us... vals) : data_{std::forward<Us>(vals)...} {}

  FixedVector& operator+=(const FixedVector& other) {
    for (size_t i = 0; i < N; ++i)
      data_[i] += other[i];
    return *this;
  }

  T& operator[](size_t i) { return data_[i]; }

  const T& operator[](size_t i) const { return data_[i]; }

  FixedVector operator+(const FixedVector& other) const {
    FixedVector output(*this);
    output += other;
    return output;
  }

  FixedVector& operator-=(const FixedVector& other) {
    for (size_t i = 0; i < N; ++i)
      data_[i] -= other[i];
    return *this;
  }

  FixedVector operator-(const FixedVector& other) const {
    FixedVector output(*this);
    output -= other;
    return output;
  }

  FixedVector& operator*=(const FixedVector& other) {
    for (size_t i = 0; i < N; ++i)
      data_[i] *= other[i];
    return *this;
  }

  FixedVector operator*(const FixedVector& other) const {
    FixedVector output(*this);
    output *= other;
    return output;
  }

  FixedVector& operator*=(T k) {
    for (size_t i = 0; i < N; ++i)
      data_[i] *= k;
    return *this;
  }

  FixedVector operator*(T k) const {
    FixedVector output(*this);
    output *= k;
    return output;
  }

  T norm1() const {
    return std::accumulate(data_.begin(), data_.end(), 0., [](auto norm, auto v) { return norm + std::abs(v); });
  }

  T norm2() const {
    return std::sqrt(std::accumulate(data_.begin(), data_.end(), 0., [](auto norm, auto v) { return norm + v * v; }));
  }

  static T dot(const FixedVector& v1, const FixedVector& v2) {
    T sum = T(0);
    for (size_t i = 0; i < N; ++i)
      sum += v1[i] * v2[i];
    return sum;
  }

  std::array<T, N>& data() & { return data_; }
  const std::array<T, N>& data() const& { return data_; }

  template <size_t I>
  decltype(auto) get() & {
    return data_[I];
  }
  template <size_t I>
  decltype(auto) get() const& {
    return data_[I];
  }

 private:
  std::array<T, N> data_;
};

using Vector2D = FixedVector<double, 2>;
using Vector3D = FixedVector<double, 3>;
using Vector4D = FixedVector<double, 4>;
using Vector2F = FixedVector<float, 2>;
using Vector3F = FixedVector<float, 3>;
using Vector4F = FixedVector<float, 4>;

template <class T, size_t N>
std::ostream& operator<<(std::ostream& os, const FixedVector<T, N>& vec) {
  for (size_t i = 0; i < N - 1; ++i)
    os << vec[i] << ',';
  return os << vec[N - 1];
}

namespace std {
template <class T, size_t N>
struct tuple_size<::FixedVector<T, N>> : tuple_size<array<T, N>> {};

template <size_t I, class T, size_t N>
struct tuple_element<I, ::FixedVector<T, N>> : tuple_element<I, array<T, N>> {};
}  // namespace std
