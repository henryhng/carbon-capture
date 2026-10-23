// adsorption sampling, capacity and statistical uncertainty
#include "runtime/03-json.hpp"

#include <boost/pending/disjoint_sets.hpp>
#include <fftw3.h>
#include <gsl/gsl_qrng.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <numeric>
#include <random>
#include <set>
#include <stdexcept>

void sample_all(Calculations& jobs, const Json& parameters, const Json& model,
                const std::map<std::string, Sheet>& hosts, std::vector<Run>& runs) {
    const auto& materials = parameters.at("surfaces");
    const bool fitted = !hosts.empty();
    if (!fitted) for (const auto& material : materials)
        if (material.value("kind", Kind::pristine) != Kind::pristine ||
            material.value("fraction", 0.) != 0.)
            throw std::runtime_error("nitrogen interactions unvalidated");
    const auto controls = parameters.at("sampling");
    const auto seeds = controls.at("seeds").get<std::vector<int>>();
    if (seeds.size() < 2) throw std::runtime_error("independent sampling seeds");
    for (const auto& material : materials) {
        std::map<double, std::set<double>> states;
        for (const auto& state : material.at("gas"))
            states[state.at("temperature_K").get<double>()].insert(state.at("pressure_Pa").get<double>());
        if (states.empty() || std::any_of(states.begin(), states.end(),
                                          [](const auto& state) { return state.second.size() < 2; }))
            throw std::runtime_error("working capacity pressures");
        const auto kind = material.value("kind", Kind::pristine);
        const int size = material.at("size");
        Sheet host;
        if (fitted) {
            host = repeat(hosts.at(material.at("name").get<std::string>()), material.at("repeats"));
        } else if (material.contains("width_A")) {
            host = slit(size, material.at("width_A"), material.at("layers"));
        } else {
            host = sheet(size, kind, {}, material.at("height_A"));
        }
        for (int seed : seeds) {
            for (size_t state = 0; state < material.at("gas").size(); ++state) {
                auto settings = controls;
                settings.erase("seeds");
                settings.update(material.at("gas")[state]);
                settings["C"] = model.at("C");
                settings["gas"] = model.at("gas");
                if (kind == Kind::pyrrolic || kind == Kind::pyridinic) settings["N"] = model.at(kind_names.at(kind));
                settings["seed"] = seed + 1000003 * state;
                settings["burn"] = material.value("burn", controls.at("burn").get<int>());
                settings["pore"] = material.contains("width_A");
                const Json request = {{"host", host}, {"settings", settings}};
                if (std::any_of(runs.begin(), runs.end(), [&](const Run& row) { return row.input == request; }))
                    continue;
                auto values = sample(host, settings, jobs.output, jobs.remaining());
                runs.push_back({material.at("name"), size, seed, settings.at("temperature_K"),
                                settings.at("pressure_Pa"), request, std::move(values)});
            }
        }
    }
}

// sobol points with a reproducible random shift
std::vector<std::array<double, 5>> sobol(size_t count, unsigned seed) {
    auto* generator = gsl_qrng_alloc(gsl_qrng_sobol, 5);
    std::mt19937 random(seed);
    std::uniform_real_distribution<double> uniform;
    std::array<double, 5> shift;
    for (auto& value : shift)
        value = uniform(random);
    std::vector<std::array<double, 5>> points(count);
    for (auto& point : points) {
        gsl_qrng_get(generator, point.data());
        for (size_t k = 0; k < point.size(); ++k)
            point[k] = std::fmod(point[k] + shift[k], 1.0);
    }
    gsl_qrng_free(generator);
    return points;
}

namespace {
std::vector<double> pool(const std::vector<Run>& rows, Stats Sample::*field) {
    if (std::any_of(rows.begin(), rows.end(), [&](const Run& row) { return (row.values.*field).ess == 0; }))
        return {};
    double mean = 0, square = 0, variance = 0;
    for (const auto& row : rows) {
        const auto& value = row.values.*field;
        mean += value.mean;
        square += value.sem * value.sem;
    }
    mean /= rows.size();
    for (const auto& row : rows)
        variance += std::pow((row.values.*field).mean - mean, 2);
    const double between = rows.size() > 1 ? std::sqrt(variance / (rows.size() * (rows.size() - 1))) : 0;
    return {mean, std::max(between, std::sqrt(square) / rows.size())};
}
}

Estimates aggregate(const std::vector<Run>& rows) {
    using Condition = std::tuple<std::string, int, double, double>;
    auto condition = [](const Run& row) {
        return Condition{row.host, row.size, row.values.xN, row.temperature};
    };
    std::map<std::pair<Condition, int>, std::vector<Run>> seeds;
    std::map<std::tuple<Condition, double, double>, std::vector<Run>> groups;
    for (const auto& row : rows) {
        seeds[{condition(row), row.seed}].push_back(row);
        groups[{condition(row), row.pressure, -1.}].push_back(row);
    }
    for (auto& [key, values] : seeds) {
        std::sort(values.begin(), values.end(), [](const Run& a, const Run& b) {
            return a.pressure < b.pressure;
        });
        const auto& released = values.front();
        for (size_t i = 1; i < values.size(); ++i) {
            auto row = values[i];
            row.values.uptake.mean -= released.values.uptake.mean;
            row.values.uptake.sem = std::hypot(row.values.uptake.sem, released.values.uptake.sem);
            groups[{condition(row), row.pressure, released.pressure}].push_back(row);
        }
    }
    Estimates values;
    for (const auto& [key, rows] : groups) {
        const auto& [state, pressure, release] = key;
        const auto& [host, size, fraction, temperature] = state;
        Estimate value{host, size, fraction, temperature, pressure, release,
                       pool(rows, &Sample::uptake), pool(rows, &Sample::cluster), pool(rows, &Sample::order)};
        (release < 0 ? values.uptake : values.working).push_back(value);
    }
    return values;
}

Distribution distribution(const std::vector<Molecule>& molecules, const std::array<Vec3, 3>& cell,
                          double plane, double layer, double cutoff) {
    Distribution result{molecules.size(), 0};
    std::vector<Vec3> centers;
    std::array<size_t, 2> faces = {0, 0};
    for (const auto& molecule : molecules) {
        Vec3 center = molecule[0];
        center[2] -= plane;
        center[2] -= std::round(center[2] / cell[2][2]) * cell[2][2];
        if (std::abs(center[2]) >= layer) continue;
        centers.push_back(center);
        ++faces[center[2] >= 0];
        Vec3 axis;
        for (int k = 0; k < 3; ++k)
            axis[k] = molecule[1][k] - molecule[2][k];
        axis = image(axis, cell);
        const double square = std::inner_product(axis.begin(), axis.end(), axis.begin(), 0.0);
        result.order_sum += (3 * axis[2] * axis[2] / square - 1) / 2;
    }
    result.adsorbed = centers.size();
    boost::disjoint_sets_with_storage<> groups(centers.size());
    for (size_t i = 0; i < centers.size(); ++i) {
        for (size_t j = 0; j < i; ++j) {
            Vec3 delta;
            for (int k = 0; k < 3; ++k)
                delta[k] = centers[i][k] - centers[j][k];
            delta = image(delta, cell);
            const double square = std::inner_product(delta.begin(), delta.end(), delta.begin(), 0.0);
            if (square < cutoff * cutoff) groups.union_set(i, j);
            const size_t bin = std::sqrt(delta[0] * delta[0] + delta[1] * delta[1]) / .2;
            if ((centers[i][2] >= 0) == (centers[j][2] >= 0) && bin < result.pair_counts.size())
                ++result.pair_counts[bin];
        }
    }
    std::vector<size_t> sizes(centers.size(), 0);
    for (size_t i = 0; i < centers.size(); ++i)
        ++sizes[groups.find_set(i)];
    for (auto size : sizes)
        result.cluster_sum += size * size;
    const double area = std::abs(cell[0][0] * cell[1][1] - cell[0][1] * cell[1][0]);
    for (auto count : faces)
        result.pair_norm += double(count) * (double(count) - 1) / (2 * area);
    return result;
}

// initial monotone paired autocorrelation sum, geyer (1992)
Stats stats(const std::vector<double>& values) {
    const size_t n = values.size();
    if (n == 0) throw std::invalid_argument("empty series");
    const double mean = std::accumulate(values.begin(), values.end(), 0.0) / n;
    size_t length = 1;
    while (length < 2 * n)
        length *= 2;
    std::vector<double> covariance(length, 0);
    double square = 0;
    for (size_t i = 0; i < n; ++i) {
        covariance[i] = values[i] - mean;
        square += covariance[i] * covariance[i];
    }
    if (square == 0) return {mean, 0, 1, double(n)};

    std::vector<std::complex<double>> spectrum(length / 2 + 1);
    auto* frequency = reinterpret_cast<fftw_complex*>(spectrum.data());
    auto forward = fftw_plan_dft_r2c_1d(length, covariance.data(), frequency, FFTW_ESTIMATE);
    auto inverse = fftw_plan_dft_c2r_1d(length, frequency, covariance.data(), FFTW_ESTIMATE);
    fftw_execute(forward);
    for (auto& value : spectrum)
        value = std::norm(value);
    fftw_execute(inverse);
    fftw_destroy_plan(forward);
    fftw_destroy_plan(inverse);

    double g = 1, previous = std::numeric_limits<double>::infinity();
    for (size_t i = 1; i + 1 < 2 * (n / 2) - 1; i += 2) {
        const double pair = (covariance[i] + covariance[i + 1]) / covariance[0];
        if (pair <= 0) break;
        previous = std::min(previous, pair);
        g += 2 * previous;
    }
    return {mean, std::sqrt(square / (n - 1) * g / n), g, n / g};
}
