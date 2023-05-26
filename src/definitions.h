#pragma once

#include <vector>

#include "math.h"

using Point = Vector4D;
using TimePoint = std::pair<double, Point>;
using Solution = std::vector<TimePoint>;

/**
 * Структура для постоянных параметров задачи
 */
struct Params {
  double alpha;
};

constexpr double T0 = 0.;  /// Начальный момент времени
constexpr double T = 1.;   /// Конечный момент времени

/// Вычисление правых частей
double f1(const Point& p, double t, Params params);  /// в выражении производной x1
double f2(const Point& p, double t, Params params);  /// в выражении производной x2
double f3(const Point& p, double t, Params params);  /// в выражении производной p1
double f4(const Point& p, double t, Params params);  /// в выражении производной p2

/// Вычисление подинтегрального выражения функционала J в момент t
double L(const Point& p, double t, Params params);

/// Вычисление ошибки в граничное условии на правой границе
double phi1(const Solution& sol, Params params);  /// для x1
double phi2(const Solution& sol, Params params);  /// для p2

/// Вычисление значения минимизируемого функционала
double J(const Solution& sol, Params params);