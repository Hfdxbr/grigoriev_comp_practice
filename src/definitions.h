#pragma once

#include <vector>

#include "math.h"

using Point = Vector4D;
using TimePoint = std::pair<double, Point>;
using Solution = std::vector<TimePoint>;

/// Структура для параметров задачи
struct Params {
  struct Shooting {
    void randomize();
    Shooting() = default;
    Shooting(const Shooting& other) { p = other.p; }
    FixedVector<double, 5> p{NAN, NAN, NAN, NAN, NAN};
    // Начальные параметры
    double& x1 = p[0];
    double& x2 = p[1];
    double& p1 = p[2];
    double& p2 = p[3];
    // Дополнительные параметры
    double& l1 = p[4];
  };
  double alpha;
  const Point& ip{NAN, NAN, NAN, NAN};
  Shooting sp;
};

constexpr double T0 = 0.;  /// Начальный момент времени
constexpr double T = 1.;   /// Конечный момент времени

/// Вычисление правых частей
double f1(const Point& p, double t, const Params& params);  /// в выражении производной x1
double f2(const Point& p, double t, const Params& params);  /// в выражении производной x2
double f3(const Point& p, double t, const Params& params);  /// в выражении производной p1
double f4(const Point& p, double t, const Params& params);  /// в выражении производной p2

/// Вычисление подинтегрального выражения функционала J в момент t
double L(const Point& p, double t, const Params& params);

/// Вычисление ошибки в граничное условии на правой границе
double phi1(const Solution& sol, const Params& params);  /// для x2
double phi2(const Solution& sol, const Params& params);  /// для p1
double phi3(const Solution& sol, const Params& params);  /// для I{x_1}dt

/// Вычисление значения минимизируемого функционала
double J(const Solution& sol, const Params& params);