// quantum espresso inputs, results and bounded native processes
#include "runtime/03-json.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>

namespace {
// codata 2006, matching the saved espresso energies
constexpr double rydberg = 13.605691932782346;
constexpr double bohr = .5291772083535413;
const std::map<int, std::pair<std::string, double>> elements = {
    {1, {"H", 1.008}}, {6, {"C", 12.011}}, {7, {"N", 14.007}}, {8, {"O", 15.999}}};

std::string read_text(const fs::path& path) {
    std::ifstream stream(path);
    if (!stream) throw std::runtime_error("missing file: " + path.string());
    return {std::istreambuf_iterator<char>(stream), {}};
}

void positions(std::istream& stream, Sheet& atoms) {
    for (size_t i = 0; i < atoms.numbers.size(); ++i) {
        std::string line, symbol;
        if (!std::getline(stream, line)) throw std::runtime_error("truncated geometry");
        std::istringstream row(line);
        Vec3 point;
        if (!(row >> symbol >> point[0] >> point[1] >> point[2]))
            throw std::runtime_error("atom coordinates");
        if (symbol.size() > 1) symbol.resize(1);
        auto element = std::find_if(elements.begin(), elements.end(),
                                    [&](const auto& item) { return item.second.first == symbol; });
        if (element == elements.end()) throw std::runtime_error("atom species");
        atoms.numbers[i] = element->first;
        atoms.positions[i] = point;
    }
}
} // namespace

std::string qe_input(const Sheet& atoms, const Json& settings, const fs::path& pseudo) {
    const bool relax = settings.value("relax", false);
    const double cutoff = settings.value("cutoff_Ry", 80.);
    const auto mesh = settings.value("mesh", std::array<int, 3>{3, 3, 1});
    const bool polarized = settings.value("spin_polarized", settings.contains("magnetization") ||
                                                             settings.contains("starting_magnetization") ||
                                                             settings.contains("spin_sites"));
    const auto functional = settings.value("input_dft", std::string("PBE"));
    const auto dispersion = settings.value("vdw_corr", std::string("grimme-d3"));
    const size_t fixed_atoms = settings.value("fixed_atoms", size_t{0});
    if (atoms.numbers.empty() || atoms.positions.size() != atoms.numbers.size() || cutoff <= 0 ||
        fixed_atoms > atoms.numbers.size() ||
        *std::min_element(mesh.begin(), mesh.end()) < 1)
        throw std::runtime_error("quantum inputs");
    std::vector<std::string> labels;
    std::map<std::string, int> species;
    std::map<std::string, double> site_moments;
    for (size_t i = 0; i < atoms.numbers.size(); ++i) {
        const auto symbol = elements.at(atoms.numbers[i]).first;
        std::string label = symbol;
        if (settings.contains("spin_sites")) {
            for (const auto& site : settings.at("spin_sites")) {
                if (site.at("atom").get<size_t>() == i) {
                    label += std::to_string(i + 1);
                    site_moments[label] = site.at("magnetization").get<double>();
                    break;
                }
            }
        }
        labels.push_back(label);
        species[label] = atoms.numbers[i];
    }
    std::ostringstream text;
    text << std::setprecision(16) << "&CONTROL\n calculation = '"
         << (relax ? "relax" : "scf") << "'\n"
         << " prefix = 'pwscf'\n outdir = '.'\n pseudo_dir = '" << pseudo.string() << "'\n"
         << " restart_mode = '" << (settings.value("restart", false) ? "restart" : "from_scratch") << "'\n"
         << " max_seconds = " << settings.value("max_seconds", 1e7) << '\n'
         << " tprnfor = .true.\n nstep = " << settings.value("steps", 200)
         << "\n etot_conv_thr = 1e-7\n forc_conv_thr = "
         << settings.value("force_eV_A", .01) / (std::sqrt(3.) * rydberg / bohr) << "\n disk_io = '"
         << (settings.value("density", false) ? "low" : "nowf") << "'\n/\n"
         << "&SYSTEM\n ibrav = 0\n nat = " << atoms.numbers.size() << "\n ntyp = " << species.size()
         << "\n ecutwfc = " << cutoff << "\n ecutrho = " << 4 * cutoff
         << "\n input_dft = '" << functional << "'\n vdw_corr = '" << dispersion << "'\n";
    if (dispersion == "grimme-d3")
        text << " dftd3_version = " << settings.value("dftd3_version", 3)
             << "\n dftd3_threebody = .false.\n";
    text << " occupations = 'smearing'\n smearing = 'gaussian'\n"
         << " degauss = " << settings.value("smearing_Ry", .003)
         << "\n nosym = .true.\n nspin = " << (polarized ? 2 : 1) << '\n';
    if (settings.contains("magnetization"))
        text << " tot_magnetization = " << settings.at("magnetization") << '\n';
    if (settings.contains("starting_magnetization")) {
        const auto moments = settings.at("starting_magnetization");
        size_t index = 0;
        for (const auto& [label, z] : species) {
            ++index;
            const auto symbol = elements.at(z).first;
            if (site_moments.count(label))
                text << " starting_magnetization(" << index << ") = " << site_moments.at(label) << '\n';
            else if (moments.contains(symbol))
                text << " starting_magnetization(" << index << ") = " << moments.at(symbol) << '\n';
        }
    } else if (!site_moments.empty()) {
        size_t index = 0;
        for (const auto& item : species) {
            ++index;
            if (site_moments.count(item.first))
                text << " starting_magnetization(" << index << ") = " << site_moments.at(item.first) << '\n';
        }
    }
    if (settings.contains("grid")) {
        const auto grid = settings.at("grid").get<std::array<int, 3>>();
        for (int k = 0; k < 3; ++k)
            text << " nr" << k + 1 << " = " << grid[k] << '\n';
    }
    text << "/\n&ELECTRONS\n conv_thr = 1e-8\n electron_maxstep = 200\n mixing_beta = 0.1\n"
         << " mixing_mode = 'plain'\n mixing_ndim = " << settings.value("mixing_ndim", 8)
         << "\n diagonalization = '" << settings.value("diagonalization", std::string("david"))
         << "'\n diago_david_ndim = 2\n/\n"
         << "&IONS\n ion_dynamics = 'bfgs'\n/\nATOMIC_SPECIES\n";
    for (const auto& [label, z] : species) {
        const auto& [symbol, mass] = elements.at(z);
        text << label << ' ' << mass << ' ' << symbol << ".pbe-dojo.UPF\n";
    }
    text << "K_POINTS ";
    if (mesh == std::array<int, 3>{1, 1, 1})
        text << "gamma\n";
    else
        text << "automatic\n" << mesh[0] << ' ' << mesh[1] << ' ' << mesh[2] << " 0 0 0\n";
    text << "CELL_PARAMETERS angstrom\n";
    for (auto vector : atoms.cell)
        text << vector[0] << ' ' << vector[1] << ' ' << vector[2] << '\n';
    text << "ATOMIC_POSITIONS angstrom\n";
    for (size_t i = 0; i < atoms.numbers.size(); ++i) {
        const auto point = atoms.positions[i];
        text << labels[i] << ' ' << point[0] << ' ' << point[1] << ' ' << point[2];
        if (i < fixed_atoms) text << " 0 0 0";
        text << '\n';
    }
    return text.str();
}

QeResult read_qe(const fs::path& output, const fs::path& input) {
    const auto data = read_text(output);
    if (data.find("JOB DONE.") == std::string::npos) throw std::runtime_error("unfinished quantum output");
    if (data.find("convergence NOT achieved") != std::string::npos)
        throw std::runtime_error("electronic convergence");
    if (data.find("Maximum CPU time exceeded") != std::string::npos ||
        data.find("maximum cpu time exceeded") != std::string::npos)
        throw std::runtime_error("compute budget");
    const auto source = read_text(input);
    std::smatch match;
    if (!std::regex_search(source, match, std::regex(R"(\bnat\s*=\s*(\d+))")))
        throw std::runtime_error("atom count");
    const int count = std::stoi(match[1]);
    QeResult result{};
    result.atoms.numbers.resize(count);
    result.atoms.positions.resize(count);
    std::istringstream source_lines(source);
    std::string line;
    while (std::getline(source_lines, line)) {
        if (line.find("CELL_PARAMETERS") == 0) {
            if (line.find("angstrom") == std::string::npos) throw std::runtime_error("cell units");
            for (auto& vector : result.atoms.cell) {
                std::getline(source_lines, line);
                std::istringstream row(line);
                if (!(row >> vector[0] >> vector[1] >> vector[2])) throw std::runtime_error("cell vectors");
            }
        } else if (line.find("ATOMIC_POSITIONS") == 0) {
            if (line.find("angstrom") == std::string::npos) throw std::runtime_error("position units");
            positions(source_lines, result.atoms);
        }
    }
    bool energy_found = false;
    std::istringstream output_lines(data);
    while (std::getline(output_lines, line)) {
        if (line.find('!') != std::string::npos && line.find("total energy") != std::string::npos) {
            result.energy = std::stod(line.substr(line.find('=') + 1)) * rydberg;
            result.forces.clear();
            energy_found = true;
        } else if (line.find("total magnetization") != std::string::npos) {
            result.magnetization = std::stod(line.substr(line.find('=') + 1));
        } else if (line.find("atom") != std::string::npos && line.find("force =") != std::string::npos) {
            std::istringstream row(line.substr(line.find("force =") + 7));
            Vec3 force;
            if (!(row >> force[0] >> force[1] >> force[2])) throw std::runtime_error("force values");
            for (auto& value : force)
                value *= rydberg / bohr;
            result.forces.push_back(force);
        } else if (line.find("ATOMIC_POSITIONS") == 0) {
            if (line.find("angstrom") == std::string::npos) throw std::runtime_error("position units");
            positions(output_lines, result.atoms);
        }
    }
    if (!energy_found || result.forces.size() != result.atoms.numbers.size())
        throw std::runtime_error("incomplete quantum result");
    for (auto force : result.forces)
        result.force = std::max(
            result.force, std::sqrt(std::inner_product(force.begin(), force.end(), force.begin(), 0.0)));
    return result;
}

QeResult quantum_in(const Sheet& atoms, const Json& settings, const fs::path& directory, double seconds) {
    const auto pseudo = fs::absolute(settings.at("pseudo").get<std::string>());
    const auto source = qe_input(atoms, settings, pseudo);
    const auto input = fs::absolute(settings.value("inputs", std::string("inputs/qe"))) /
                       digest(source) / "qe.in";
    fs::create_directories(input.parent_path());
    std::ofstream stream(input);
    stream << source;
    stream.close();
    if (!stream) throw std::runtime_error("write quantum input");
    const int ranks = settings.value("ranks", 1);
    if (ranks < 1 || ranks > 12) throw std::runtime_error("mpi ranks");
    std::vector<std::string> command = {settings.value("executable", std::string("pw.x"))};
    if (ranks > 1) command.insert(command.begin(), {"mpirun", "-np", std::to_string(ranks)});
    run(command, directory, input.string(), "qe.out", seconds);
    auto result = read_qe(directory / "qe.out", input);
    const size_t fixed_atoms = settings.value("fixed_atoms", size_t{0});
    if (settings.value("relax", false) && fixed_atoms > 0) {
        result.force = 0;
        for (size_t i = fixed_atoms; i < result.forces.size(); ++i) {
            const auto force = result.forces[i];
            result.force = std::max(
                result.force, std::sqrt(std::inner_product(force.begin(), force.end(), force.begin(), 0.0)));
        }
    }
    if (settings.value("relax", false) &&
        read_text(directory / "qe.out").find("bfgs converged") == std::string::npos)
        throw std::runtime_error("ionic convergence");
    if (settings.value("relax", false) && result.force >= settings.value("force_eV_A", .01))
        throw std::runtime_error("geometry convergence");
    return result;
}

QeResult quantum(const Sheet& atoms, const Json& settings, const fs::path& output, double seconds) {
    if (seconds < 1) throw std::runtime_error("compute budget");
    const auto directory = output / "checkpoint";
    Json potentials = Json::object();
    for (int z : std::set<int>(atoms.numbers.begin(), atoms.numbers.end())) {
        const auto symbol = elements.at(z).first;
        const auto file = fs::path(settings.at("pseudo").get<std::string>()) / (symbol + ".pbe-dojo.UPF");
        potentials[symbol] = digest(read_text(file));
    }
    auto options = settings;
    options.erase("inputs");
    const Json request = {{"atoms", atoms}, {"settings", options}, {"pseudopotentials", potentials}};
    if (fs::exists(directory / "input.json") && read_json(directory / "input.json") != request)
        throw std::runtime_error("checkpoint input changed");
    options = settings;
    options["restart"] = fs::exists(directory / "pwscf.save/data-file-schema.xml");
    if (options.at("restart").get<bool>()) std::cerr << "Resuming QE checkpoint\n";
    options["density"] = true;
    if (std::isfinite(seconds)) options["max_seconds"] = std::max(.1, seconds - 1);
    write_json(directory / "input.json", request);
    const auto result = quantum_in(atoms, options, directory, seconds);
    fs::remove_all(directory);
    return result;
}
