#include <algorithm>
#include <fstream>
#include <iomanip>
#include <sstream>

#include "definitions.h"
#include "runge_kutta.h"

void print_points(const Solution& sol, const Params& params, const std::string& filename) {
  std::ofstream ofs(filename);
  ofs << "t,x1,x2,p1,p2,L" << std::endl;
  std::for_each(sol.begin(), sol.end(), [&](const TimePoint& tp) {
    const auto& [t, p] = tp;
    ofs << t << ',' << p << ',' << L(p, t, params) << std::endl;
  });
};

template <class V, class M = decltype(std::defaultfloat)>
std::string to_string(const V& v, std::streamsize precision = 6, M modifier = std::defaultfloat) {
  std::stringstream ss;
  ss << std::setprecision(precision) << modifier << v;
  return ss.str();
}

void execute(const std::vector<double>& alphas, Params& params) {
  std::ofstream ofs_stats("data/stats.csv");
  ofs_stats << std::setprecision(8) << std::fixed;
  ofs_stats << "\\(\\alpha\\),\\(\\lambda_1\\),\\(x_1(1)\\),\\(p_1(0)\\),\\(p_2(0)\\),\\(J\\),\\(error\\)" << std::endl;

  std::for_each(alphas.begin(), alphas.end(), [&](double alpha) {
    double error = 0;
    params.alpha = alpha;
    Solution sol = shooting(params, error);
    auto alpha_s = to_string(alpha, 2, std::fixed);
    while (alpha_s.back() == '0' && alpha_s.size() > 3) alpha_s.pop_back();
    std::string plot_file = "data/points_alpha" + alpha_s + ".csv";
    print_points(sol, params, plot_file);
    auto& [x1, x2, _1, _2] = sol.back().second;
    auto& [_3, _4, p10, p20, l1] = params.sp.p;
    ofs_stats << alpha_s << ',' << Vector4D{l1, x1, p10, p20} << ',' << J(sol, params) << ','
              << to_string(error, 3, std::scientific) << std::endl;
  });
}

int main() {
  Params params{.ip = Point{0., 1., NAN, NAN}};
  params.sp.p1 = 1, params.sp.p2 = 1, params.sp.l1 = 1;
  std::vector<double> alphas{0.0, 0.01, 0.5, 1.5, 10.5};
  execute(alphas, params);

  return 0;
}