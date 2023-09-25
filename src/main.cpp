#include <array>
#include <fstream>
#include <iostream>

#include "params.h"

using Point = std::array<double, 4>;

Point operator*(double k, const Point& p) {
  auto out = p;
  for (auto& e : out) e *= k;
  return out;
}

Point& operator+=(Point& p1, const Point& p2) {
  for (int i = 0; i < p1.size(); ++i) p1[i] += p2[i];
  return p1;
}

Point operator+(const Point& p1, const Point& p2) {
  auto out = p1;
  out += p2;
  return out;
}

double dist(const Point& p1, const Point& p2) {
  double distance = 0;
  for (int i = 0; i < p1.size(); ++i) distance += std::pow(p1[i] - p2[i], 2);

  return std::sqrt(distance);
}

const double initial_step = std::pow(tolerance, 0.25);

void runge_kutta(Point& p, double h) {
  double hh = 0.5 * h, sh = h / 6;
  Point k1, k2, k3, k4;
  {
    auto [x1, x2, p1, p2] = p;
    k1 = {f1(x1, x2, p1, p2), f2(x1, x2, p1, p2), f3(x1, x2, p1, p2), f4(x1, x2, p1, p2)};
  }
  {
    auto [x1, x2, p1, p2] = p + hh * k1;
    k2 = {f1(x1, x2, p1, p2), f2(x1, x2, p1, p2), f3(x1, x2, p1, p2), f4(x1, x2, p1, p2)};
  }
  {
    auto [x1, x2, p1, p2] = p + hh * k2;
    k3 = {f1(x1, x2, p1, p2), f2(x1, x2, p1, p2), f3(x1, x2, p1, p2), f4(x1, x2, p1, p2)};
  }
  {
    auto [x1, x2, p1, p2] = p + h * k3;
    k4 = {f1(x1, x2, p1, p2), f2(x1, x2, p1, p2), f3(x1, x2, p1, p2), f4(x1, x2, p1, p2)};
  }
  p += sh * (k1 + 2 * (k2 + k3) + k4);
}

double detect_u_breaking_step(const Point& start_point, double h_current, State state, Point& final_point) {
  // методом половинного деления ищем точку разрыва в управлении
  double lo = 0, hi = h_current;
  Point p;
  while (hi - lo > tolerance) {
    p = start_point;
    double mid = 0.5 * (lo + hi);
    runge_kutta(p, mid);
    if (u_state != state) {
      hi = mid;
    } else {
      lo = mid;
    }
  }

  final_point = p;

  return 0.5 * (lo + hi);
}

double solve(double a, double b, Point* last, double* err1, double* err2, std::ostream* ofs, std::ostream* ocs) {
  Point p = {x10, x20, a, b};
  u(x10, x20, a, b);  // выставляем начальный u_state
  double t = 0, err = 0, I = 0;
  unsigned counter = 0;
  if (ofs) *ofs << "t,x1,x2,p1,p2\n" << t << "," << p[0] << "," << p[1] << "," << p[2] << "," << p[3] << "\n";
  double h = initial_step;
  while (t < T) {
    Point p1 = p, p2 = p;
    State cur_state = u_state;

    runge_kutta(p1, h);
    State new_state = u_state;
    if (new_state != cur_state)  //  ищем шаг до разрыва
      h = detect_u_breaking_step(p, h, cur_state, p1);

    double hh = 0.5 * h;
    runge_kutta(p2, hh);
    runge_kutta(p2, hh);

    double error = dist(p1, p2);
    if (error > tolerance && h > tolerance) {  // в случае большой ошибки уполовиниваем шаг
      h = hh;
      u_state = cur_state;
      continue;
    }

    err += error;

    p = p1;
    u_state = new_state;  // выставляем в качестве текущего состояние после точки разрыва (если он был)
    t += h;
    I += h * L(p[0], p[1], p[2], p[3]);
    if (ofs) *ofs << t << "," << p[0] << "," << p[1] << "," << p[2] << "," << p[3] << "\n";
    ++counter;
    h = std::min(initial_step, T - t);
  }
  double x1Err = phi1(p[0], p[1], p[2], p[3]), p2Err = phi2(p[0], p[1], p[2], p[3]);
  if (ocs) *ocs << x1Err << "," << p2Err << "," << I << "," << err << std::endl;

  if (err1) *err1 = x1Err;
  if (err2) *err2 = p2Err;
  if (last) *last = p;

  return std::sqrt(std::pow(x1Err, 2) + std::pow(p2Err, 2));  // суммарная ошибка
}

void shooting(double& x10, double& p20, std::ostream* ocs) {
  double err1, err2, a = x10, b = p20;
  unsigned counter = 0;

  double err = solve(a, b, nullptr, &err1, &err2, nullptr, ocs);
  for (; err > tolerance; err = solve(a, b, nullptr, &err1, &err2, nullptr, ocs)) {
    double da = 0.0001 * a;
    double db = 0.0001 * b;

    double err1a, err2a, err1b, err2b;
    solve(a + da, b, nullptr, &err1a, &err2a, nullptr, nullptr);
    solve(a, b + db, nullptr, &err1b, &err2b, nullptr, nullptr);

    double W11, W12, W21, W22, det;
    W11 = (err1a - err1) / da;
    W12 = (err1b - err1) / db;
    W21 = (err2a - err2) / da;
    W22 = (err2b - err2) / db;
    det = W11 * W22 - W12 * W21;

    if (fabs(det) < 1e-15) {
      exit(1);
    } else {
      a = a - (W22 * err1 - W12 * err2) / det;
      b = b - (W11 * err2 - W21 * err1) / det;
    }
    ++counter;
  }
  x10 = a;
  p20 = b;
}

int main(int argc, char* argv[]) {
  set_final_time(argv[1]);
  double a = 0.1, b = 0.15;
  shooting(a, b, nullptr);
  auto T_str = std::to_string(T);
  T_str.erase(T_str.find('.') + 2, T_str.size());
  std::ofstream ofs("data/T" + T_str + ".csv");
  std::ofstream stats("data/stats.csv", std::ios_base::app);
  stats << T << "," << a << "," << b << ",";
  solve(a, b, nullptr, nullptr, nullptr, &ofs, &stats);
}
