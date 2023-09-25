#pragma once
#include <cmath>
#include <string>

constexpr double tolerance = 1e-10;

constexpr double x10 = 0;
constexpr double x20 = 0;
constexpr double x1T = 0;
constexpr double p2T = 0;

double T;

void set_final_time(const std::string& arg) { T = std::stod(arg); }

enum State { LowerBound = -1, Intermediate, UpperBound };
State u_state;

double u(double x1, double x2, double p1, double p2) {
  if (p2 < -1) {
    u_state = LowerBound;
    return -1;
  }
  if (1 < p2) {
    u_state = UpperBound;
    return 1;
  }
  u_state = Intermediate;
  return p2;
}

double f1(double x1, double x2, double p1, double p2) { return x2; }

double f2(double x1, double x2, double p1, double p2) { return u(x1, x2, p1, p2); }

double f3(double x1, double x2, double p1, double p2) { return -x1; }

double f4(double x1, double x2, double p1, double p2) { return -p1 - x2; }

double phi1(double x1, double x2, double p1, double p2) { return x1 - x1T; }

double phi2(double x1, double x2, double p1, double p2) { return p2 - p2T; }

double L(double x1, double x2, double p1, double p2) {
  return std::pow(u(x1, x2, p1, p2), 2) - std::pow(x2, 2) - std::pow(x1, 2);
}