#include "definitions.h"

#include <iostream>
#include <numeric>
#include <random>

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<double> dis(-10.0, 10.0);

void Params::Shooting::randomize() {
  for (int i = 0; i < p.size(); ++i)
    if (!std::isnan(p[i])) p[i] = dis(gen);
}

double f1(const Point& p, double t, const Params& params) {
  const auto& [x1, x2, p1, p2] = p;
  return x2;
}

double f2(const Point& p, double t, const Params& params) {
  const auto& [x1, x2, p1, p2] = p;
  return p2 * (1. + params.alpha * std::pow(t, 4));
}

double f3(const Point& p, double t, const Params& params) {  //
  return params.sp.l1;
}

double f4(const Point& p, double t, const Params& params) {
  const auto& [x1, x2, p1, p2] = p;
  return -p1;
}

double L(const Point& p, double t, const Params& params) {
  const auto& [x1, x2, p1, p2] = p;
  return p2 * p2 * (1. + params.alpha * std::pow(t, 4));
}

double phi1(const Solution& sol, const Params& params) {
  const auto& [x1, x2, p1, p2] = sol.back().second;
  return x2 - 0.;
}

double phi2(const Solution& sol, const Params& params) {
  const auto& [x1, x2, p1, p2] = sol.back().second;
  return p1 - 0.;
}

double phi3(const Solution& sol, const Params& params) {
  double sum = 0;
  for (int i = 1; i < sol.size() - 1; ++i) {
    const auto& [x1, x2, p1, p2] = sol[i].second;
    sum += 0.5 * x1 * (sol[i + 1].first - sol[i - 1].first);
  }
  {
    const auto& [x1, x2, p1, p2] = sol[0].second;
    sum += 0.5 * x1 * (sol[1].first - sol[0].first);
  }
  {
    int n = sol.size() - 1;
    const auto& [x1, x2, p1, p2] = sol[n].second;
    sum += 0.5 * x1 * (sol[n].first - sol[n - 1].first);
  }
  return sum - 1.;
}

double J(const Solution& sol, const Params& params) {
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