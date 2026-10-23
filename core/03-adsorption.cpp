// adsorption energies and interaction targets
#include "runtime/03-json.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <stdexcept>

namespace {
struct CellCandidate {
    QeResult host, bound;
    Json quantum;
    Kind kind;
    std::array<int, 2> site;
    int size;
    double lattice, height;
    int direction;
    std::string spin_seed;
};
} // namespace

Adsorption adsorption(Calculations& jobs, const Json& settings) {
    Adsorption data;
    std::map<std::string, CellCandidate> candidates;
    data.scope = {{"geometry", "one CO2 per periodic cell"},
                  {"uptake_or_selectivity_or_capacity", false},
                  {"gcmc_or_fitted_force_field", false},
                  {"binding_gate_is_workflow_only", true}};
    const Json fixed = {{"relax", false}};
    const Json relax = {{"relax", true}};
    const double tolerance = settings.at("convergence_eV");
    const int count = settings.at("relax_count");
    bool convergence_checked = false;
    if (count < 0 || settings.at("relax_names").empty() ||
        size_t(count) > settings.at("relax_names").size())
        throw std::runtime_error("relaxation count");
    if (count == 0 && settings.contains("fit")) throw std::runtime_error("screening fit");
    QeResult gas_relaxed{};

    for (const auto& material : settings.at("hosts")) {
        const std::string name = material.at("name");
        const auto kind = material.at("kind").get<Kind>();
        const auto site = material.at("site").get<std::array<int, 2>>();
        const double a = material.value("lattice_A", 2.46), height = material.at("height_A");
        auto host_changes = material.value("quantum", Json::object());
        const auto ideal = sheet(material.at("size"), kind,
                                 kind == Kind::pristine ? std::vector<std::array<int, 2>>{} :
                                                     std::vector<std::array<int, 2>>{site}, height, a);
        const Vec3 vacancy_site = {a * (site[0] - site[1] / 2.),
                                   a * site[1] * std::sqrt(3.) / 2., height / 2.};
        std::vector<VacancySeed> geometries;
        if (kind == Kind::vacancy && material.value("reconstruction_directions", 0) > 0) {
            geometries = vacancy_seeds(ideal, vacancy_site,
                                       material.value("reconstruction_seed_A", .1));
            geometries.resize(std::min<size_t>(geometries.size(),
                                                material.value("reconstruction_directions", 3)));
        } else {
            geometries.push_back({ideal, {0, 0}, -1});
        }
        const Json spin_seeds = material.value("spin_seeds", Json::array({Json{{"name", "default"}}}));
        QeResult host{};
        double host_minimum = std::numeric_limits<double>::infinity();
        int host_direction = -1;
        std::string host_spin_seed;
        Json host_trials = Json::array();
        for (const auto& geometry : geometries) {
            for (const auto& spin_seed : spin_seeds) {
                auto trial_options = host_changes;
                trial_options.update(spin_seed.value("quantum", Json::object()));
                trial_options.erase("spin_sites");
                if (spin_seed.value("name", std::string()) == "afm") {
                    const double moment = spin_seed.value("moment", .3);
                    trial_options["spin_sites"] = {
                        {{"atom", geometry.pair[0]}, {"magnetization", moment}},
                        {{"atom", geometry.pair[1]}, {"magnetization", -moment}}};
                }
                trial_options["relax"] = true;
                std::string label = name + "/host";
                if (geometry.direction >= 0)
                    label += "/recon-" + std::to_string(geometry.direction);
                if (spin_seeds.size() > 1) label += "/" + spin_seed.at("name").get<std::string>();
                const auto trial = jobs.solve(label, geometry.atoms, trial_options);
                const bool reconstructed = kind == Kind::vacancy &&
                                           vacancy_reconstructed(trial.atoms, vacancy_site);
                host_trials.push_back({{"direction", geometry.direction},
                                       {"spin_seed", spin_seed.value("name", std::string("default"))},
                                       {"energy_eV", trial.energy},
                                       {"magnetization_muB", trial.magnetization},
                                       {"reconstructed_5_9", reconstructed}});
                if (trial.energy < host_minimum) {
                    host = trial;
                    host_minimum = trial.energy;
                    host_direction = geometry.direction;
                    host_spin_seed = spin_seed.value("name", std::string("default"));
                    host_changes = trial_options;
                    host_changes.erase("relax");
                }
            }
        }
        auto surface_fixed = host_changes;
        surface_fixed.update(fixed);
        const Vec3 anchor = {a * (site[0] - site[1] / 2.) + a / 4.,
                             a * site[1] * std::sqrt(3.) / 2. + a / (4. * std::sqrt(3.)), height / 2.};
        const Sheet empty{{}, {}, host.atoms.cell};
        const auto gas_atoms = attach(empty, {co2(anchor, {1, 0, 0})});
        if (count > 0)
            gas_relaxed = jobs.solve(name + "/relaxed gas", gas_atoms,
                                     Json{{"relax", true}, {"mesh", {1, 1, 1}}});
        data.hosts[name] = host.atoms;
        if ((kind == Kind::pyridinic || kind == Kind::pyrrolic) && !topology(host.atoms, kind, 1.8))
            throw std::runtime_error("surface topology");

        auto poses = settings.at("configurations");
        poses["anchor_A"] = anchor;
        auto local = configurations(host.atoms, kind, poses);
        for (auto& point : local) point.label = name + "/" + point.label;
        std::vector<Point> starts;
        for (const std::string label : settings.at("relax_names")) {
            const auto point = std::find_if(local.begin(), local.end(), [&](const Point& p) {
                return p.label == name + "/" + label;
            });
            if (point == local.end() || point->molecules.size() != 1)
                throw std::runtime_error("relaxation pose");
            const auto bound = jobs.solve(point->label, attach(host.atoms, point->molecules), surface_fixed);
            const auto isolated = jobs.solve(point->label + "/gas", attach(empty, point->molecules),
                                              {{"mesh", {1, 1, 1}}});
            starts.push_back(*point);
            starts.back().target = bound.energy - host.energy - isolated.energy;
            if (!std::isfinite(starts.back().target)) throw std::runtime_error("screening energy");
        }
        std::stable_sort(starts.begin(), starts.end(), [](const Point& a, const Point& b) {
            return a.target < b.target;
        });
        if (count == 0) {
            data.points.insert(data.points.end(), starts.begin(), starts.end());
            continue;
        }

        std::vector<Point> selected;
        std::set<std::string> selected_labels;
        for (const auto& group : settings.value("relax_groups", Json::array())) {
            if (selected.size() >= size_t(count)) break;
            auto group_best = starts.end();
            for (const auto& label : group) {
                const auto point = std::find_if(starts.begin(), starts.end(), [&](const Point& item) {
                    return item.label == name + "/" + label.get<std::string>();
                });
                if (point != starts.end() &&
                    (group_best == starts.end() || point->target < group_best->target))
                    group_best = point;
            }
            if (group_best != starts.end() && selected_labels.insert(group_best->label).second)
                selected.push_back(*group_best);
        }
        for (const auto& point : starts) {
            if (selected.size() >= size_t(count)) break;
            if (selected_labels.insert(point.label).second) selected.push_back(point);
        }
        QeResult best{};
        double minimum = std::numeric_limits<double>::infinity();
        Json relaxed_poses = Json::array();
        for (const auto& point : selected) {
            auto bound_options = host_changes;
            bound_options.update(relax);
            bound_options["fixed_atoms"] = host.atoms.numbers.size();
            const auto bound = jobs.solve(point.label + "/relax", attach(host.atoms, point.molecules),
                                          bound_options);
            if ((kind == Kind::pyridinic || kind == Kind::pyrrolic) &&
                !topology(fragment(bound.atoms, 0, host.atoms.numbers.size()), kind, 1.8))
                throw std::runtime_error("surface topology");
            const double energy = bound.energy - host.energy - gas_relaxed.energy;
            relaxed_poses.push_back(point.label.substr(name.size() + 1));
            if (energy < minimum) { minimum = energy; best = bound; }
        }
        if (!std::isfinite(minimum)) throw std::runtime_error("missing relaxed adsorption");
        data.binding[name] = minimum;
        candidates[name] = {host, best, host_changes, kind, site, material.at("size"), a, height,
                            host_direction, host_spin_seed};
        double surface_z = -std::numeric_limits<double>::infinity();
        for (const auto& position : host.atoms.positions) surface_z = std::max(surface_z, position[2]);
        const double distance = best.atoms.positions.at(host.atoms.numbers.size())[2] - surface_z;
        const bool retained = distance <= settings.value("desorption_A", 4.0);
        const bool promising = retained && minimum < -settings.value("binding_gate_eV", .05);
        data.gates[name] = {{"adsorbate_distance_A", distance}, {"retained", retained},
                            {"binding_gate_eV", -settings.value("binding_gate_eV", .05)},
                            {"proceed_to_convergence", promising},
                            {"host_spin_seed", host_spin_seed},
                            {"host_trials", host_trials},
                            {"relaxed_poses", relaxed_poses},
                            {"cell_size", material.at("size")},
                            {"periodic_defect_repeats", kind != Kind::pristine}};
        if (kind == Kind::vacancy)
            data.gates[name]["vacancy_reconstruction"] =
                vacancy_reconstructed(host.atoms, vacancy_site);
        if (material.value("compare_spin", false)) {
            auto nonspin = host_changes;
            nonspin["spin_polarized"] = false;
            nonspin.erase("starting_magnetization");
            nonspin.erase("magnetization");
            nonspin.erase("spin_sites");
            nonspin["relax"] = false;
            const auto host_nonspin = jobs.solve(name + "/spin-check/host", host.atoms, nonspin);
            const auto bound_nonspin = jobs.solve(name + "/spin-check/bound", best.atoms, nonspin);
            const double nonspin_binding = bound_nonspin.energy - host_nonspin.energy - gas_relaxed.energy;
            const double spin_shift = minimum - nonspin_binding;
            data.gates[name]["spin_polarized_binding_eV"] = minimum;
            data.gates[name]["spin_unpolarized_binding_eV"] = nonspin_binding;
            data.gates[name]["adsorption_spin_shift_eV"] = spin_shift;
            data.gates[name]["spin_changes_binding"] = std::abs(spin_shift) > .1;
            data.gates[name]["host_spin_energy_shift_eV"] = host.energy - host_nonspin.energy;
            data.gates[name]["host_moment_muB"] = host.magnetization;
            data.gates[name]["bound_moment_muB"] = best.magnetization;
            if (kind == Kind::vacancy) {
                const auto seeds = vacancy_seeds(ideal, vacancy_site,
                                                 material.value("reconstruction_seed_A", .1));
                const auto pair = seeds.at(std::max(host_direction, 0)).pair;
                double moment = .3;
                for (const auto& seed : spin_seeds)
                    if (seed.value("name", std::string()) == "afm") moment = seed.value("moment", .3);
                const bool selected_afm = host_spin_seed == "afm";
                auto alternate = host_changes;
                std::string alternate_name = "afm";
                if (selected_afm) {
                    alternate.erase("spin_sites");
                    alternate_name = "ferro";
                } else {
                    alternate["spin_sites"] = {{{"atom", pair[0]}, {"magnetization", moment}},
                                                {{"atom", pair[1]}, {"magnetization", -moment}}};
                }
                const auto host_alt = jobs.solve(name + "/spin-check/" + alternate_name + "-host",
                                                 host.atoms, alternate);
                const auto bound_alt = jobs.solve(name + "/spin-check/" + alternate_name + "-bound",
                                                  best.atoms, alternate);
                const double alternate_binding = bound_alt.energy - host_alt.energy - gas_relaxed.energy;
                data.gates[name][alternate_name + "_binding_eV"] = alternate_binding;
                data.gates[name][alternate_name + "_host_energy_shift_eV"] = host_alt.energy - host.energy;
                data.gates[name][alternate_name + "_host_moment_muB"] = host_alt.magnetization;
                data.gates[name][alternate_name + "_adsorption_shift_eV"] = alternate_binding - minimum;
            }
        }

        // wac18 eq. 1: matched references, unchanged geometries between methods
        for (const auto& method : settings.at("methods")) {
            if (method.value("gate", false) && !promising) continue;
            if (method.contains("hosts") &&
                std::find(method.at("hosts").begin(), method.at("hosts").end(), name) == method.at("hosts").end())
                continue;
            const std::string method_name = method.at("name");
            auto options = method.at("quantum");
            options["relax"] = false;
            auto gas_options = options;
            gas_options["mesh"] = {1, 1, 1};
            auto surface_options = options;
            surface_options.update(host_changes);
            const auto h = jobs.solve(name + "/" + method_name + "/host", host.atoms, surface_options);
            const auto g = jobs.solve(name + "/" + method_name + "/gas", gas_relaxed.atoms, gas_options);
            const auto b = jobs.solve(name + "/" + method_name + "/bound", best.atoms, surface_options);
            const double energy = b.energy - h.energy - g.energy;
            data.methods[name][method_name] = energy;
            if (method.value("convergence", false)) {
                convergence_checked = true;
                if (!(std::abs(energy - minimum) <= tolerance))
                    throw std::runtime_error("adsorption convergence");
            }
        }

        if (!settings.contains("fit")) continue;

        if (material.value("holdout", false) || kind == Kind::vacancy)
            for (auto& point : local) point.holdout = true;
        // zhang 2023, doi:10.1088/1674-1056/ac8ce2; co2-specific fit targets
        for (auto& point : local) {
            const auto bound = jobs.solve(point.label, attach(host.atoms, point.molecules), surface_fixed);
            double gas_energy;
            if (point.molecules.size() == 2) {
                gas_energy = jobs.solve(point.label + "/gas", attach(empty, point.molecules), fixed).energy;
            } else
                gas_energy = jobs.solve(point.label + "/gas", attach(empty, point.molecules),
                                        {{"mesh", {1, 1, 1}}}).energy;
            point.target = bound.energy - host.energy - gas_energy;
            data.points.push_back(point);
            if (point.molecules.size() == 2) {
                if (!point.holdout) throw std::runtime_error("pair must be held out");
                std::array<double, 2> single{}, isolated{};
                for (size_t i = 0; i < 2; ++i) {
                    const auto label = point.label + "/single-" + std::to_string(i);
                    single[i] = jobs.solve(label, attach(host.atoms, {point.molecules[i]}), surface_fixed).energy;
                    isolated[i] = jobs.solve(label + "/gas", attach(Sheet{{}, {}, host.atoms.cell},
                                                        {point.molecules[i]}), fixed).energy;
                }
                // ding 2023, section 3.2: substrate plus two adsorbates
                const double pair = bound.energy + host.energy - single[0] - single[1];
                const double gas_pair = gas_energy - isolated[0] - isolated[1];
                for (const auto observable : {Observable::pair, Observable::mediated}) {
                    auto target = point;
                    target.observable = observable;
                    target.label += "/" + observable_names.at(observable);
                    target.target = observable == Observable::pair ? pair : pair - gas_pair;
                    data.points.push_back(target);
                }
            }
        }
    }
    const bool ungated_convergence = std::any_of(
        settings.at("methods").begin(), settings.at("methods").end(),
        [](const Json& method) { return method.value("convergence", false) && !method.value("gate", false); });
    if (count > 0 && ungated_convergence && !convergence_checked)
        throw std::runtime_error("missing convergence checks");
    if (settings.contains("cell_check_size") && !candidates.empty()) {
        const int size = settings.at("cell_check_size");
        std::string pristine, defect;
        for (const auto& [name, candidate] : candidates) {
            if (candidate.kind == Kind::pristine) pristine = name;
            else if (defect.empty() || data.binding.at(name) < data.binding.at(defect)) defect = name;
        }
        std::vector<std::string> names;
        if (!pristine.empty()) names.push_back(pristine);
        if (!defect.empty()) names.push_back(defect);
        for (const auto& name : names) {
            const auto& candidate = candidates.at(name);
            const auto defects = candidate.kind == Kind::pristine ? std::vector<std::array<int, 2>>{} :
                                                                    std::vector<std::array<int, 2>>{candidate.site};
            auto large_ideal = sheet(size, candidate.kind, defects, candidate.height, candidate.lattice);
            const Vec3 site = {candidate.lattice * (candidate.site[0] - candidate.site[1] / 2.),
                               candidate.lattice * candidate.site[1] * std::sqrt(3.) / 2.,
                               candidate.height / 2.};
            auto options = candidate.quantum;
            options.erase("spin_sites");
            if (candidate.kind == Kind::vacancy) {
                double displacement = .1;
                for (const auto& material : settings.at("hosts"))
                    if (material.at("name") == name)
                        displacement = material.value("reconstruction_seed_A", displacement);
                const auto seeds = vacancy_seeds(large_ideal, site,
                                                 displacement);
                const auto direction = static_cast<size_t>(std::max(candidate.direction, 0));
                large_ideal = seeds.at(direction).atoms;
                if (candidate.spin_seed == "afm") {
                    double moment = .3;
                    for (const auto& material : settings.at("hosts"))
                        if (material.at("name") == name)
                            for (const auto& seed : material.value("spin_seeds", Json::array()))
                                if (seed.value("name", std::string()) == "afm")
                                    moment = seed.value("moment", .3);
                    const auto pair = seeds.at(direction).pair;
                    options["spin_sites"] = {{{"atom", pair[0]}, {"magnetization", moment}},
                                              {{"atom", pair[1]}, {"magnetization", -moment}}};
                }
            }
            auto mesh = jobs.settings.value("mesh", std::array<int, 3>{3, 3, 1});
            mesh[0] = std::max(1, static_cast<int>(std::lround(mesh[0] * candidate.size / double(size))));
            mesh[1] = std::max(1, static_cast<int>(std::lround(mesh[1] * candidate.size / double(size))));
            options["mesh"] = mesh;
            options["relax"] = true;
            const auto large_host = jobs.solve(name + "/cell-" + std::to_string(size) + "/host",
                                               large_ideal, options);
            double old_top = -std::numeric_limits<double>::infinity();
            double new_top = -std::numeric_limits<double>::infinity();
            for (const auto& position : candidate.host.atoms.positions) old_top = std::max(old_top, position[2]);
            for (const auto& position : large_host.atoms.positions) new_top = std::max(new_top, position[2]);
            Molecule molecule{};
            for (size_t i = 0; i < molecule.size(); ++i) {
                molecule[i] = candidate.bound.atoms.positions[candidate.host.atoms.numbers.size() + i];
                molecule[i][2] += new_top - old_top;
            }
            auto bound_options = options;
            bound_options["fixed_atoms"] = large_host.atoms.numbers.size();
            const auto large_bound = jobs.solve(name + "/cell-" + std::to_string(size) + "/bound",
                                                attach(large_host.atoms, {molecule}), bound_options);
            auto gas_options = options;
            gas_options.erase("starting_magnetization");
            gas_options.erase("magnetization");
            gas_options.erase("spin_sites");
            gas_options["spin_polarized"] = false;
            gas_options["mesh"] = {1, 1, 1};
            gas_options["relax"] = true;
            const auto gas = jobs.solve(name + "/cell-" + std::to_string(size) + "/gas",
                                        attach(Sheet{{}, {}, large_host.atoms.cell}, {molecule}), gas_options);
            const double energy = large_bound.energy - large_host.energy - gas.energy;
            data.cell_checks[name] = {{"size", size}, {"binding_eV", energy},
                                      {"change_from_original_eV", energy - data.binding.at(name)},
                                      {"defect_repeats", candidate.kind != Kind::pristine}};
        }
    }
    if (settings.contains("fit") && settings.contains("gas_pairs")) {
        auto poses = settings.at("gas_pairs");
        const double length = poses.at("cell_A");
        Sheet empty{{}, {}, {{{length, 0, 0}, {0, length, 0}, {0, 0, length}}}};
        poses["anchor_A"] = {length / 2, length / 2, length / 2};
        auto options = poses.at("quantum");
        options["relax"] = false;
        for (auto point : configurations(empty, Kind::gas, poses)) {
            point.target = jobs.solve("gas/" + point.label, attach(empty, point.molecules), options).energy;
            for (size_t i = 0; i < point.molecules.size(); ++i)
                point.target -= jobs.solve("gas/" + point.label + "/" + std::to_string(i),
                                            attach(empty, {point.molecules[i]}), options).energy;
            data.points.push_back(point);
        }
    }
    return data;
}

void binding(Calculations& jobs, const Json& atoms, const Json& configurations, std::vector<Binding>& values) {
    values.clear();
    for (const auto& settings : configurations) {
        const auto changes = settings.value("quantum", Json::object());
        auto options = jobs.settings;
        options.update(changes);
        const auto prefix = changes.empty() ? std::string() : digest(changes.dump()) + "/";
        const auto axis = settings.at("axis").get<Vec3>();
        const auto orientation = settings.at("axis").dump();
        for (const std::string name : settings.at("hosts")) {
            const auto kind = Json(name).get<Kind>();
            const auto host = jobs.solve(prefix + name, atoms.at(name).get<Sheet>(), changes);
            Vec3 center{};
            for (auto vector : host.atoms.cell)
                for (int k = 0; k < 3; ++k)
                    center[k] += vector[k] / 2;
            const auto gas = jobs.solve(prefix + name + "/gas/" + orientation,
                                        attach(Sheet{{}, {}, host.atoms.cell}, {co2(center, axis)}), changes);
            const auto anchor = host.atoms.positions[kind == Kind::pristine ? host.atoms.positions.size() / 2
                                                                        : host.atoms.positions.size() - 1];
            for (double height : settings.at("heights_A")) {
                center = anchor;
                center[2] += height;
                const auto bound = jobs.solve(prefix + name + "/" + std::to_string(height) + "/" + orientation,
                                              attach(host.atoms, {co2(center, axis)}), changes);
                values.push_back({name, height, axis, options, bound.energy - host.energy - gas.energy, bound});
            }
        }
    }
}
