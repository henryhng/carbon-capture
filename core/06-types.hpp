// structures, numerical results and interfaces
#pragma once

#include <array>
#include <cstddef>
#include <string>
#include <vector>
#include <filesystem>
#include <map>
#include "json.hpp"

namespace fs = std::filesystem;
using Json = nlohmann::json;

using Vec3 = std::array<double, 3>;
using Molecule = std::array<Vec3, 3>;

enum class Kind { pristine, vacancy, pyridinic, pyrrolic, gas };
struct Sheet {
    std::vector<int> numbers;
    std::vector<Vec3> positions;
    std::array<Vec3, 3> cell;
};

Vec3 image(Vec3 delta, const std::array<Vec3, 3>& cell);
std::vector<double> charges(const Sheet& atoms, double q);
Sheet sheet(int n, Kind kind, const std::vector<std::array<int, 2>>& sites, double height,
            double a = 2.46);
Sheet slit(int n, double width, int layers);
Sheet repeat(const Sheet& atoms, int count);
Molecule co2(Vec3 center, Vec3 axis);

struct Surface {
    std::vector<Vec3> carbon;
    std::array<Vec3, 3> cell;
    Vec3 sigma, epsilon, shift;
    double cutoff;

    Surface(const Sheet& atoms, double cutoff = 12, double epsilon_k = 28, double sigma_a = 3.4,
            double gas_scale = 1);
    double energy(Molecule molecule) const;
};

double adsorption_length(const std::vector<std::array<double, 5>>& points, double height, double temperature,
                         double cutoff = 12, double epsilon = 28, double sigma = 3.4, double gas_scale = 1);

struct Stats {
    double mean, sem, inefficiency, ess;
};

struct Distribution {
    std::size_t count, adsorbed;
    double cluster_sum = 0, order_sum = 0, pair_norm = 0;
    std::vector<double> pair_counts = std::vector<double>(50, 0);
};

Distribution distribution(const std::vector<Molecule>& molecules, const std::array<Vec3, 3>& cell,
                          double plane, double layer, double cutoff);
Stats stats(const std::vector<double>& values);
struct Sample {
    Stats count{}, adsorbed{}, uptake{}, cluster{}, order{};
    std::vector<double> radial, excess;
    size_t frames = 0;
    double xN = 0, drift = 0, drift_sem = 0;
    bool converged = false;
};

struct Run {
    std::string host;
    int size = 0, seed = 0;
    double temperature = 0, pressure = 0;
    Json input;
    Sample values;
};

struct Estimate {
    std::string host;
    int size;
    double xN, temperature, pressure, release;
    std::vector<double> uptake, cluster, order;
};
struct Estimates {
    std::vector<Estimate> uptake, working;
};
Estimates aggregate(const std::vector<Run>& rows);

std::vector<std::array<double, 5>> sobol(size_t count, unsigned seed);
Json read_json(const fs::path& path);
void write_json(const fs::path& path, const Json& value);
std::string digest(const std::string& text);
struct Scratch {
    fs::path path;
    explicit Scratch(const fs::path& parent);
    ~Scratch();
    Scratch(const Scratch&) = delete;
    Scratch& operator=(const Scratch&) = delete;
};
void run(const std::vector<std::string>& args, const fs::path& directory, const std::string& input,
         const std::string& output, double seconds, double memory_mb = 6144);

struct QeResult {
    double energy, force;
    Sheet atoms;
    std::vector<Vec3> forces;
    double magnetization = 0;
};

std::string qe_input(const Sheet& atoms, const Json& settings, const fs::path& pseudo);
QeResult read_qe(const fs::path& output, const fs::path& input);
QeResult quantum(const Sheet& atoms, const Json& settings, const fs::path& output, double seconds);
Sample sample(const Sheet& host, const Json& settings, const fs::path& directory, double seconds);
