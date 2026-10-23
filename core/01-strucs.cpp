// graphene, defect arrangements and rigid molecules
#include "core/07-pipeline.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <set>
#include <stdexcept>

Vec3 image(Vec3 delta, const std::array<Vec3, 3>& cell) {
    const double area = cell[0][0] * cell[1][1] - cell[0][1] * cell[1][0];
    const double u = std::round((delta[0] * cell[1][1] - delta[1] * cell[1][0]) / area);
    const double v = std::round((cell[0][0] * delta[1] - cell[0][1] * delta[0]) / area);
    delta[2] -= std::round(delta[2] / cell[2][2]) * cell[2][2];
    Vec3 best = delta;
    double minimum = std::numeric_limits<double>::infinity();
    for (int i = -1; i <= 1; ++i) {
        for (int j = -1; j <= 1; ++j) {
            Vec3 candidate = delta;
            for (int k = 0; k < 2; ++k)
                candidate[k] -= (u + i) * cell[0][k] + (v + j) * cell[1][k];
            const double square =
                std::inner_product(candidate.begin(), candidate.end(), candidate.begin(), 0.0);
            if (square < minimum) {
                minimum = square;
                best = candidate;
            }
        }
    }
    return best;
}

// neutral local charge model: nitrogen and its two carbon neighbors
std::vector<double> charges(const Sheet& atoms, double q) {
    std::vector<double> result(atoms.numbers.size(), 0);
    for (size_t i = 0; i < atoms.numbers.size(); ++i) {
        if (atoms.numbers[i] != 7) continue;
        std::vector<std::pair<double, size_t>> neighbors;
        for (size_t j = 0; j < atoms.numbers.size(); ++j) {
            if (atoms.numbers[j] != 6) continue;
            Vec3 delta;
            for (int k = 0; k < 3; ++k)
                delta[k] = atoms.positions[i][k] - atoms.positions[j][k];
            delta = image(delta, atoms.cell);
            neighbors.emplace_back(std::inner_product(delta.begin(), delta.end(), delta.begin(), 0.0), j);
        }
        if (neighbors.size() < 2) throw std::invalid_argument("nitrogen neighbors");
        std::partial_sort(neighbors.begin(), neighbors.begin() + 2, neighbors.end());
        result[i] += q;
        for (int j = 0; j < 2; ++j)
            result[neighbors[j].second] -= q / 2;
    }
    return result;
}

Sheet sheet(int n, Kind kind, const std::vector<std::array<int, 2>>& sites, double height,
            double a) {
    if (n < 1 || height <= 0) throw std::invalid_argument("cell size");
    if (kind == Kind::gas)
        throw std::invalid_argument("defect kind");
    if (kind == Kind::pristine && !sites.empty()) throw std::invalid_argument("pristine defects");

    Sheet atoms;
    atoms.cell = {{{n * a, 0, 0}, {-n * a / 2, n * a * std::sqrt(3) / 2, 0}, {0, 0, height}}};
    std::vector<Vec3> carbon;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            const double x = a * (i - j / 2.0), y = j * a * std::sqrt(3) / 2;
            carbon.push_back({x, y, height / 2});
            carbon.push_back({x + a / 2, y + a / (2 * std::sqrt(3)), height / 2});
        }
    }

    std::set<int> removed;
    std::vector<Vec3> nitrogen;
    for (auto site : sites) {
        const int i = (site[0] % n + n) % n, j = (site[1] % n + n) % n;
        const int first = 2 * (i * n + j);
        removed.insert(first);
        if (kind != Kind::vacancy) removed.insert(first + 1);
        Vec3 position = carbon[first];
        if (kind == Kind::pyrrolic) {
            const double shift = .5;
            position[0] += a / 4 - shift / 2;
            position[1] += a / (4 * std::sqrt(3)) + shift * std::sqrt(3) / 2;
        }
        if (kind != Kind::vacancy) nitrogen.push_back(position);
    }
    for (size_t i = 0; i < carbon.size(); ++i) {
        if (!removed.count(i)) {
            atoms.positions.push_back(carbon[i]);
            atoms.numbers.push_back(6);
        }
    }
    for (auto position : nitrogen) {
        atoms.positions.push_back(position);
        atoms.numbers.push_back(7);
    }
    return atoms;
}

Sheet repeat(const Sheet& atoms, int count) {
    if (count < 1) throw std::invalid_argument("cell repeats");
    Sheet result{{}, {}, atoms.cell};
    for (int k = 0; k < 2; ++k)
        for (auto& value : result.cell[k]) value *= count;
    for (int i = 0; i < count; ++i)
        for (int j = 0; j < count; ++j)
            for (size_t a = 0; a < atoms.numbers.size(); ++a) {
                auto position = atoms.positions[a];
                for (int k = 0; k < 3; ++k)
                    position[k] += i * atoms.cell[0][k] + j * atoms.cell[1][k];
                result.numbers.push_back(atoms.numbers[a]);
                result.positions.push_back(position);
            }
    return result;
}

Sheet slit(int n, double width, int layers) {
    if (width <= 0 || layers < 1) throw std::invalid_argument("pore dimensions");
    const auto plane = sheet(n, Kind::pristine, {}, width);
    Sheet atoms;
    atoms.cell = plane.cell;
    atoms.cell[2][2] *= layers;
    for (int i = 0; i < layers; ++i) {
        for (auto position : plane.positions) {
            position[2] += i * width;
            atoms.positions.push_back(position);
            atoms.numbers.push_back(6);
        }
    }
    return atoms;
}

Molecule co2(Vec3 center, Vec3 axis) {
    const double norm = std::sqrt(std::inner_product(axis.begin(), axis.end(), axis.begin(), 0.0));
    if (norm == 0) throw std::invalid_argument("zero axis");
    Molecule positions = {center, center, center};
    for (int k = 0; k < 3; ++k) {
        positions[1][k] += 1.16 * axis[k] / norm;
        positions[2][k] -= 1.16 * axis[k] / norm;
    }
    return positions;
}

Sheet fragment(const Sheet& atoms, size_t begin, size_t end) {
    if (begin >= end || end > atoms.numbers.size()) throw std::runtime_error("fragment range");
    return {{atoms.numbers.begin() + begin, atoms.numbers.begin() + end},
            {atoms.positions.begin() + begin, atoms.positions.begin() + end},
            atoms.cell};
}

Sheet attach(Sheet host, const std::vector<Molecule>& molecules) {
    for (const auto& molecule : molecules) {
        host.numbers.insert(host.numbers.end(), {6, 8, 8});
        host.positions.insert(host.positions.end(), molecule.begin(), molecule.end());
    }
    return host;
}

bool topology(const Sheet& atoms, Kind kind, double bond) {
    if (kind == Kind::pristine) return std::count(atoms.numbers.begin(), atoms.numbers.end(), 7) == 0;
    std::vector<std::vector<size_t>> graph(atoms.numbers.size());
    for (size_t i = 0; i < graph.size(); ++i)
        for (size_t j = 0; j < i; ++j) {
            Vec3 delta{};
            for (int k = 0; k < 3; ++k)
                delta[k] = atoms.positions[i][k] - atoms.positions[j][k];
            delta = image(delta, atoms.cell);
            const double square = std::inner_product(delta.begin(), delta.end(), delta.begin(), 0.);
            if (square < .64) return false;
            if (square < bond * bond) {
                graph[i].push_back(j);
                graph[j].push_back(i);
            }
        }
    bool found = false;
    for (size_t nitrogen = 0; nitrogen < graph.size(); ++nitrogen) {
        if (atoms.numbers[nitrogen] != 7) continue;
        found = true;
        if (graph[nitrogen].size() != 2) return false;
        std::vector<int> distance(graph.size(), -1);
        distance[nitrogen] = 0;
        const auto start = graph[nitrogen][0], end = graph[nitrogen][1];
        distance[start] = 0;
        std::vector<size_t> queue{start};
        for (size_t i = 0; i < queue.size(); ++i) {
            const auto current = queue[i];
            for (auto neighbor : graph[current])
                if (distance[neighbor] < 0) {
                    distance[neighbor] = distance[current] + 1;
                    queue.push_back(neighbor);
                }
        }
        if (distance[end] + 2 != (kind == Kind::pyrrolic ? 5 : 6)) return false;
    }
    return found;
}

namespace {
std::array<size_t, 3> vacancy_neighbors(const Sheet& atoms, Vec3 site) {
    std::vector<std::pair<double, size_t>> nearest;
    for (size_t i = 0; i < atoms.numbers.size(); ++i) {
        if (atoms.numbers[i] != 6) continue;
        Vec3 delta{};
        for (int k = 0; k < 3; ++k) delta[k] = atoms.positions[i][k] - site[k];
        delta = image(delta, atoms.cell);
        nearest.emplace_back(std::inner_product(delta.begin(), delta.end(), delta.begin(), 0.0), i);
    }
    if (nearest.size() < 3) throw std::runtime_error("vacancy neighbors");
    std::partial_sort(nearest.begin(), nearest.begin() + 3, nearest.end());
    return {nearest[0].second, nearest[1].second, nearest[2].second};
}
} // namespace

std::vector<VacancySeed> vacancy_seeds(const Sheet& atoms, Vec3 site, double displacement_A) {
    const auto neighbors = vacancy_neighbors(atoms, site);
    const std::array<std::array<size_t, 2>, 3> pairs = {
        std::array<size_t, 2>{neighbors[0], neighbors[1]},
        std::array<size_t, 2>{neighbors[1], neighbors[2]},
        std::array<size_t, 2>{neighbors[2], neighbors[0]}};
    std::vector<VacancySeed> seeds;
    for (size_t i = 0; i < pairs.size(); ++i) {
        auto seed = atoms;
        Vec3 delta{};
        for (int k = 0; k < 3; ++k)
            delta[k] = seed.positions[pairs[i][1]][k] - seed.positions[pairs[i][0]][k];
        delta = image(delta, seed.cell);
        const double length = std::sqrt(std::inner_product(delta.begin(), delta.end(), delta.begin(), 0.0));
        for (int k = 0; k < 3; ++k) {
            const double shift = displacement_A * delta[k] / length;
            seed.positions[pairs[i][0]][k] += shift;
            seed.positions[pairs[i][1]][k] -= shift;
        }
        seeds.push_back({std::move(seed), pairs[i], static_cast<int>(i)});
    }
    return seeds;
}

bool vacancy_reconstructed(const Sheet& atoms, Vec3 site) {
    constexpr double bond_A = 2.05;
    const auto neighbors = vacancy_neighbors(atoms, site);
    for (size_t i = 0; i < neighbors.size(); ++i) {
        for (size_t j = i + 1; j < neighbors.size(); ++j) {
            Vec3 delta{};
            for (int k = 0; k < 3; ++k)
                delta[k] = atoms.positions[neighbors[j]][k] - atoms.positions[neighbors[i]][k];
            delta = image(delta, atoms.cell);
            const double length = std::sqrt(std::inner_product(delta.begin(), delta.end(), delta.begin(), 0.0));
            if (length < bond_A) return true;
        }
    }
    return false;
}

std::vector<Point> configurations(const Sheet& host, Kind kind, const Json& settings) {
    std::vector<Point> points;
    const auto anchor = settings.contains("anchor_A") ? settings.at("anchor_A").get<Vec3>() :
        host.positions.at(kind == Kind::pristine ? host.positions.size() / 2 : host.positions.size() - 1);
    for (const auto& item : settings.at("single")) {
        std::vector<Molecule> molecules;
        for (const auto& pose : item.contains("molecules") ? item.at("molecules") : Json::array({item})) {
            auto center = anchor;
            const auto offset = pose.at("offset_A").get<Vec3>();
            for (int k = 0; k < 3; ++k) center[k] += offset[k];
            if (pose.value("site", std::string()) == "top") {
                Vec3 nearest_delta{};
                double nearest_square = std::numeric_limits<double>::infinity();
                for (const auto& position : host.positions) {
                    const auto delta = image({center[0] - position[0], center[1] - position[1], 0}, host.cell);
                    const double square = std::inner_product(delta.begin(), delta.end(), delta.begin(), 0.0);
                    if (square < nearest_square) {
                        nearest_square = square;
                        nearest_delta = delta;
                    }
                }
                center[0] -= nearest_delta[0];
                center[1] -= nearest_delta[1];
            }
            molecules.push_back(co2(center, pose.at("axis").get<Vec3>()));
        }
        points.push_back({kind, item.at("name"), host, molecules, 0, item.at("holdout")});
    }
    return points;
}
