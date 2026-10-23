// molecular interactions and dilute adsorption
#include "core/06-types.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

constexpr double pi = 3.14159265358979323846;
constexpr double kb = 1.380649e-23 / 1.602176634e-19;

// trappe carbon dioxide, potoff and siepmann, doi:10.1002/aic.690470719
constexpr Vec3 sigma_co2 = {2.80, 3.05, 3.05};
constexpr Vec3 epsilon_co2 = {27, 79, 79};

// periodic pristine surface, shifted lennard-jones cutoff
Surface::Surface(const Sheet& atoms, double cutoff, double epsilon_k, double sigma_a, double gas_scale)
    : cell(atoms.cell), cutoff(cutoff) {
    if (std::any_of(atoms.numbers.begin(), atoms.numbers.end(), [](int z) { return z != 6; }))
        throw std::invalid_argument("pristine surface required");
    if (cutoff <= 0 || epsilon_k < 0 || gas_scale < 0) throw std::invalid_argument("interaction parameters");
    for (int k = 0; k < 3; ++k) {
        sigma[k] = (sigma_a + sigma_co2[k]) / 2;
        epsilon[k] = std::sqrt(epsilon_k * epsilon_co2[k] * gas_scale) * kb;
        const double sixth = std::pow(sigma[k] / cutoff, 6);
        shift[k] = sixth * sixth - sixth;
    }
    const double area = cell[0][0] * cell[1][1];
    const double reach = cutoff + 1.16 + 2 * cell[0][0];
    const int repeats = std::ceil(reach * cell[0][0] / area);
    for (int i = -repeats; i <= repeats; ++i) {
        for (int j = -repeats; j <= repeats; ++j) {
            for (auto position : atoms.positions) {
                for (int k = 0; k < 3; ++k)
                    position[k] += i * cell[0][k] + j * cell[1][k];
                carbon.push_back(position);
            }
        }
    }
}

double Surface::energy(Molecule molecule) const {
    const double fy = molecule[0][1] / cell[1][1];
    const double fx = (molecule[0][0] - fy * cell[1][0]) / cell[0][0];
    for (auto& position : molecule) {
        for (int k = 0; k < 2; ++k)
            position[k] -= std::floor(fx) * cell[0][k] + std::floor(fy) * cell[1][k];
    }

    double total = 0;
    for (int i = 0; i < 3; ++i) {
        if (epsilon[i] == 0) continue;
        for (auto position : carbon) {
            double r2 = 0;
            for (int k = 0; k < 3; ++k)
                r2 += std::pow(molecule[i][k] - position[k], 2);
            if (r2 >= cutoff * cutoff) continue;
            if (r2 == 0) return std::numeric_limits<double>::infinity();
            const double sixth = std::pow(sigma[i] * sigma[i] / r2, 3);
            total += 4 * epsilon[i] * (sixth * sixth - sixth - shift[i]);
        }
    }
    return total;
}
// points: fractional x, y, z, polar cosine, azimuth in the unit cube
double adsorption_length(const std::vector<std::array<double, 5>>& points, double height, double temperature,
                         double cutoff, double epsilon, double sigma, double gas_scale) {
    if (points.empty() || temperature <= 0) throw std::invalid_argument("quadrature inputs");
    const auto atoms = sheet(1, Kind::pristine, {}, height);
    const Surface surface(atoms, cutoff, epsilon, sigma, gas_scale);
    double integral = 0;
    for (const auto& p : points) {
        const Vec3 center = {p[0] * atoms.cell[0][0] + p[1] * atoms.cell[1][0], p[1] * atoms.cell[1][1],
                             p[2] * height};
        const double cosine = 2 * p[3] - 1, sine = std::sqrt(std::max(0.0, 1 - cosine * cosine));
        const Vec3 axis = {sine * std::cos(2 * pi * p[4]), sine * std::sin(2 * pi * p[4]), cosine};
        integral += std::exp(-surface.energy(co2(center, axis)) / (kb * temperature));
    }
    return height * integral / points.size();
}
