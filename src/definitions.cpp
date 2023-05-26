#include "definitions.h"

#include <iostream>
#include <numeric>

double f1(const Point& p, double t, Params params) {
  const auto& [x1, x2, p1, p2] = p;
  return x2;
}

double f2(const Point& p, double t, Params params) {
  const auto& [x1, x2, p1, p2] = p;
  return p2;
}

double f3(const Point& p, double t, Params params) {
  const auto& [x1, x2, p1, p2] = p;
  double denom = 2. + std::cos(params.alpha * x1);
  return -24. * params.alpha * x2 * std::sin(params.alpha * x1) / denom / denom;
}

double f4(const Point& p, double t, Params params) {
  const auto& [x1, x2, p1, p2] = p;
  double denom = 2. + std::cos(params.alpha * x1);
  return -24. / denom - p1;
}

double L(const Point& p, double t, Params params) {
  const auto& [x1, x2, p1, p2] = p;
  return p2 * p2 - 48. * x2 / (2. + std::cos(params.alpha * x1));
}

double phi1(const Solution& sol, Params params) {
  const auto& [x1, x2, p1, p2] = sol.back().second;
  return x1 - 0.;
}

double phi2(const Solution& sol, Params params) {
  const auto& [x1, x2, p1, p2] = sol.back().second;
  return p2 - 0.;
}

double J(const Solution& sol, Params params) {
  const auto& [x10, x20, p10, p20] = sol.front().second;
  auto sum = [L_prev = (p20 * p20), t_prev = 0., params](double s, const TimePoint& tp) mutable {
    const auto& [t_curr, p] = tp;
    double L_curr = L(p, t_curr, params);
    s += (L_curr + L_prev) * (t_curr - t_prev) * 0.5;
    t_prev = t_curr;
    L_prev = L_curr;
    return s;
  };

  return std::accumulate(++sol.begin(), sol.end(), 0., sum);
}