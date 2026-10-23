// calculation state and pipeline interfaces
#pragma once
#include "model/02-model.hpp"
#include <chrono>

struct Calculations {
    fs::path output;
    Json settings, cache;
    double seconds;
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();

    Calculations(const fs::path& output, const Json& settings, double seconds);
    double remaining() const;
    QeResult solve(const std::string& name, const Sheet& atoms, const Json& changes = Json::object());
};

Sheet fragment(const Sheet& atoms, size_t begin, size_t end);
Sheet attach(Sheet host, const std::vector<Molecule>& molecules);
bool topology(const Sheet& atoms, Kind kind, double bond);
struct VacancySeed {
    Sheet atoms;
    std::array<size_t, 2> pair;
    int direction;
};
std::vector<VacancySeed> vacancy_seeds(const Sheet& atoms, Vec3 site, double displacement_A);
bool vacancy_reconstructed(const Sheet& atoms, Vec3 site);
std::vector<Point> configurations(const Sheet& host, Kind kind, const Json& settings);
struct Binding {
    std::string host;
    double height_A;
    Vec3 axis;
    Json quantum;
    double binding_eV;
    QeResult bound;
};
void binding(Calculations& jobs, const Json& atoms, const Json& settings, std::vector<Binding>& values);
void sample_all(Calculations& jobs, const Json& parameters, const Json& model,
                const std::map<std::string, Sheet>& hosts, std::vector<Run>& runs);
struct Adsorption {
    Json scope;
    std::map<std::string, Sheet> hosts;
    std::map<std::string, double> binding;
    std::map<std::string, Json> gates;
    std::map<std::string, Json> cell_checks;
    std::map<std::string, std::map<std::string, double>> methods;
    std::vector<Point> points;
};
Adsorption adsorption(Calculations& jobs, const Json& settings);
