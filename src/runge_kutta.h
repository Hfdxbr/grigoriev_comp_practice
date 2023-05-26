#pragma once

#include "definitions.h"
#include "math.h"

using Coeffs = Vector4D;
using ErrVec = Vector4D;

constexpr double eps = 1e-15;  /// Минимальная возможная ошибка
constexpr double eps_step = 1e-13;  /// Допустимая ошибка при вычислении величины шага
constexpr double eps_boundary = 1e-11;  /// Допустимая ошибка при вычислении граничных условий
constexpr double max_allowed_step();  /// Ограничение на размер шага вычисления

/// Вычисление нормы вектора
double norm(const Vector4D& v);

/// Вычисление коэффициентов для метода Р-К
Coeffs calculate_k(Point p, double t, const Params& params, const Coeffs& k = {});

/// Вычисление значения точки на следующем шаге и оценки ошибки
/// Использует метод Р-К 5го порядка
std::pair<Point, ErrVec> gammas_with_error(const Point& p, double t, const Params& params, double h);

/// Вычисление нового значения шага
double update_step(const ErrVec& err, double h);

/// Поиск решения при заданных начальных значениях x1 и p2
Solution solve(double x10, double p20, const Params& params, double& error);

/// Итеративный поиск решения методом стрельбы
Solution shooting(double& x10, double& p20, const Params& params, double& error);