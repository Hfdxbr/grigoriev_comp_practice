#include "runge_kutta.h"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <random>

#include "definitions.h"

constexpr double max_allowed_step() { return (T - T0) / 100.; }

double norm(const Vector4D& v) { return v.norm1(); }

Coeffs calculate_k(Point p, double t, const Params& params, const Coeffs& k) {
  p += k;
  return {f1(p, t, params), f2(p, t, params), f3(p, t, params), f4(p, t, params)};
}

constexpr double f1_6 = 1.0 / 6.0;
constexpr double f2_3 = 2.0 / 3.0;
constexpr double f1_24 = 1.0 / 24.0;
constexpr double f1_27 = 1.0 / 27.0;
constexpr double f5_48 = 5.0 / 48.0;
constexpr double f27_56 = 27.0 / 56.0;
constexpr double fm1_336 = -1.0 / 336.0;
constexpr double f125_336 = 125.0 / 336.0;
constexpr double f1_625 = 1.0 / 625.0;
std::pair<Point, ErrVec> gammas_with_error(const Point& p, double t, const Params& params, double h) {
  Coeffs k1 = calculate_k(p, t, params) * h;
  Coeffs k2 = calculate_k(p, t + h * 0.5, params, k1 * 0.5) * h;
  Coeffs k3 = calculate_k(p, t + h * 0.5, params, (k1 + k2) * 0.25) * h;
  Coeffs k4 = calculate_k(p, t + h, params, k3 * 2. - k2) * h;
  Coeffs k5 = calculate_k(p, t + h * f2_3, params, (k1 * 7. + k2 * 10. + k4) * f1_27) * h;
  Coeffs k6 =
    calculate_k(p, t + h * 0.2, params, (k1 * 28. - k2 * 125. + k3 * 546. + k4 * 54. - k5 * 378.) * f1_625) * h;

  Point new_p = p + k1 * f1_24 + k4 * f5_48 + k5 * f27_56 + k6 * f125_336;
  ErrVec err = (k1 * 42. + k3 * 224. + k4 * 21. - k5 * 162. - k6 * 125.) * fm1_336;
  return std::make_pair(new_p, err);
}

double update_step(const ErrVec& err, double h) {
  return std::clamp(std::pow(eps_step / (eps + norm(err)), f1_6), 0.1, 10.) * 0.95 * h;
}

Solution solve(double x10, double p20, const Params& params, double& error) {
  double t = 0;
  double h = eps * 100.;
  Solution sol = {TimePoint{0, {x10, 0, 0, p20}}};
  while (t < T) {
    auto [p, err] = gammas_with_error(sol.back().second, t, params, h);
    if (double error = norm(err); error > eps_step || error < 0.1 * eps_step) {
      h = update_step(err, h);
      continue;
    }

    bool need_recalc = false;
    if (h > max_allowed_step()) {
      h = max_allowed_step();
      need_recalc = true;
    }

    if (t + h >= T) {
      h = T - t;
      need_recalc = true;
    }

    if (need_recalc)
      std::tie(p, err) = gammas_with_error(sol.back().second, t, params, h);

    error += norm(err);
    t += h;
    sol.emplace_back(t, p);
  }
  return sol;
}

Solution shooting(double& x10, double& p20, const Params& params, double& error) {
  double dx10 = 0.001;
  double dp20 = 0.001;

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> dis(-10.0, 10.0);

  while (true) {
    Solution sol = solve(x10, p20, params, error);
    Vector2D phi{phi1(sol, params), phi2(sol, params)};

    if (phi.norm1() < eps_boundary) {
      return sol;
    }

    double W11, W12, W21, W22;
    {
      double error = 0;
      Solution sol_ = solve(x10 + dx10, p20, params, error);
      Vector2D dphi = Vector2D{phi1(sol_, params), phi2(sol_, params)} - phi;
      W11 = dphi[0] / dx10;
      W21 = dphi[1] / dx10;
    }
    {
      double error = 0;
      Solution sol_ = solve(x10, p20 + dp20, params, error);
      Vector2D dphi = Vector2D{phi1(sol_, params), phi2(sol_, params)} - phi;
      W12 = dphi[0] / dx10;
      W22 = dphi[1] / dx10;
    }

    double detW = W11 * W22 - W12 * W21;
    if (detW < eps) {
      x10 = dis(gen);
      p20 = dis(gen);
      continue;
    }

    x10 -= Vector2D::dot(Vector2D{W22, -W12}, phi) / detW;
    p20 -= Vector2D::dot(Vector2D{-W21, W11}, phi) / detW;
  }
}