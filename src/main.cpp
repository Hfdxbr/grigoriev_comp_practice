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

void execute(double x10, double p20, const std::vector<Params>& params) {
  std::ofstream ofs_stats("data/stats.csv");
  ofs_stats << std::setprecision(8) << std::fixed;
  ofs_stats << "\\(\\alpha\\),\\(x_1(0)\\),\\(p_2(0)\\),\\(x_2(1)\\),\\(p_1(1)\\),\\(J\\),\\(error\\)" << std::endl;

  std::for_each(params.begin(), params.end(), [&](const Params& params) {
    double error = 0;
    Solution sol = shooting(x10, p20, params, error);
    auto alpha_s = to_string(params.alpha, 1, std::fixed);
    std::string plot_file = "data/points_alpha" + alpha_s + ".csv";
    print_points(sol, params, plot_file);
    ofs_stats << alpha_s << ',' << sol.back().second << ',' << J(sol, params) << ','
              << to_string(error, 3, std::scientific) << std::endl;
  });
}

int main() {
  std::vector<Params> params{{0.0}, {0.1}, {1.0}, {5.1}};
  double x10 = 1., p20 = -1.;
  execute(x10, p20, params);

  return 0;
}