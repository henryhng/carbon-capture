// interaction calibration and validation
#pragma once
#include "core/06-types.hpp"
#include <functional>

struct Fit {
    std::vector<double> values;
    double rmse;
    bool converged, at_bound;
};
Fit least_squares(const std::vector<double>& initial, const std::vector<std::array<double, 2>>& bounds,
                  const std::function<std::vector<double>(const std::vector<double>&)>& residual,
                  size_t iterations, double tolerance);
double molecular_energy(const Sheet& host, const std::vector<Molecule>& molecules, const Json& settings,
                        const fs::path& output, double seconds);

enum class Observable { adsorption, pair, mediated };
struct Point {
    Kind kind;
    std::string label;
    Sheet host;
    std::vector<Molecule> molecules;
    double target;
    bool holdout;
    Observable observable = Observable::adsorption;
};

struct Model {
    Json parameters = Json::object();
    std::map<std::string, Fit> fits;
    std::vector<double> errors;
};

void calibrate(const std::vector<Point>& points, const Json& options, const Json& settings,
               const fs::path& output, const std::function<double()>& remaining, Model& model);
