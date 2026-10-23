// bounded least squares using gsl trust regions
#include "runtime/03-json.hpp"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multifit_nlinear.h>
#include <algorithm>
#include <cmath>
#include <exception>
#include <stdexcept>

namespace {
struct Residual {
    const std::vector<std::array<double, 2>>& bounds;
    const std::function<std::vector<double>(const std::vector<double>&)>& function;
    std::exception_ptr error;
};

std::vector<double> unpack(const gsl_vector* x, const std::vector<std::array<double, 2>>& bounds) {
    std::vector<double> values(x->size);
    for (size_t i = 0; i < values.size(); ++i) {
        const double t = gsl_vector_get(x, i);
        const double fraction = t >= 0 ? 1 / (1 + std::exp(-t)) : std::exp(t) / (1 + std::exp(t));
        values[i] = bounds[i][0] + (bounds[i][1] - bounds[i][0]) * fraction;
    }
    return values;
}

int evaluate(const gsl_vector* x, void* data, gsl_vector* f) {
    auto& residual = *static_cast<Residual*>(data);
    try {
        const auto values = residual.function(unpack(x, residual.bounds));
        if (values.size() != f->size) throw std::runtime_error("residual size");
        for (size_t i = 0; i < values.size(); ++i) {
            if (!std::isfinite(values[i])) throw std::runtime_error("nonfinite residual");
            gsl_vector_set(f, i, values[i]);
        }
        return GSL_SUCCESS;
    } catch (...) {
        residual.error = std::current_exception();
        return GSL_EBADFUNC;
    }
}
} // namespace

Fit least_squares(const std::vector<double>& initial, const std::vector<std::array<double, 2>>& bounds,
                  const std::function<std::vector<double>(const std::vector<double>&)>& residual,
                  size_t iterations, double tolerance) {
    const auto values = residual(initial);
    if (initial.empty() || initial.size() != bounds.size() || values.size() < initial.size())
        throw std::runtime_error("fit dimensions");
    std::vector<double> transformed(initial.size());
    for (size_t i = 0; i < initial.size(); ++i) {
        if (!(bounds[i][0] < initial[i] && initial[i] < bounds[i][1])) throw std::runtime_error("fit bounds");
        const double fraction = (initial[i] - bounds[i][0]) / (bounds[i][1] - bounds[i][0]);
        transformed[i] = std::log(fraction / (1 - fraction));
    }
    gsl_set_error_handler_off();
    Residual data{bounds, residual, {}};
    gsl_multifit_nlinear_fdf fdf{};
    fdf.f = evaluate;
    fdf.n = values.size();
    fdf.p = initial.size();
    fdf.params = &data;
    auto settings = gsl_multifit_nlinear_default_parameters();
    settings.h_df = 1e-3;
    auto* work =
        gsl_multifit_nlinear_alloc(gsl_multifit_nlinear_trust, &settings, values.size(), initial.size());
    if (!work) throw std::runtime_error("fit allocation");
    auto start = gsl_vector_view_array(transformed.data(), transformed.size());
    int status = gsl_multifit_nlinear_init(&start.vector, &fdf, work), info = 0;
    if (status == GSL_SUCCESS)
        status = gsl_multifit_nlinear_driver(iterations, tolerance, tolerance, tolerance, nullptr, nullptr,
                                             &info, work);
    Fit fit{unpack(gsl_multifit_nlinear_position(work), bounds), 0, status == GSL_SUCCESS, false};
    auto* errors = gsl_multifit_nlinear_residual(work);
    for (size_t i = 0; i < errors->size; ++i)
        fit.rmse += std::pow(gsl_vector_get(errors, i), 2);
    fit.rmse = std::sqrt(fit.rmse / errors->size);
    for (size_t i = 0; i < fit.values.size(); ++i) {
        const double fraction = (fit.values[i] - bounds[i][0]) / (bounds[i][1] - bounds[i][0]);
        fit.at_bound |= fraction < 1e-3 || fraction > 1 - 1e-3;
    }
    gsl_multifit_nlinear_free(work);
    if (data.error) std::rethrow_exception(data.error);
    return fit;
}

void calibrate(const std::vector<Point>& points, const Json& options, const Json& settings,
               const fs::path& output, const std::function<double()>& remaining, Model& model) {
    for (const auto& a : points)
        for (const auto& b : points)
            if (a.holdout != b.holdout && Json(a.host) == Json(b.host) && a.molecules == b.molecules)
                throw std::runtime_error("training holdout overlap");
    if (model.parameters.empty()) model.parameters = options.at("initial");
    auto& parameters = model.parameters;
    auto predict = [&](const Point& point, const Json& trial) {
        auto controls = settings;
        controls["C"] = trial.at("C");
        controls["gas"] = trial.at("gas");
        if (std::count(point.host.numbers.begin(), point.host.numbers.end(), 7))
            controls["N"] = trial.at(kind_names.at(point.kind));
        const auto& host = point.host;
        Sheet empty{{}, {}, host.cell};
        double energy = molecular_energy(host, point.molecules, controls, output, remaining());
        if (point.observable == Observable::pair || point.observable == Observable::mediated) {
            if (point.molecules.size() != 2) throw std::runtime_error("pair configuration");
            energy += molecular_energy(host, {}, controls, output, remaining());
            for (const auto& molecule : point.molecules)
                energy -= molecular_energy(host, {molecule}, controls, output, remaining());
            if (point.observable == Observable::mediated) {
                energy -= molecular_energy(empty, point.molecules, controls, output, remaining());
                for (const auto& molecule : point.molecules)
                    energy += molecular_energy(empty, {molecule}, controls, output, remaining());
            }
        } else if (point.kind == Kind::gas) {
            for (auto molecule : point.molecules)
                energy -= molecular_energy(empty, {molecule}, controls, output, remaining());
        } else {
            energy -= molecular_energy(empty, point.molecules, controls, output, remaining());
        }
        return energy;
    };
    for (const auto& term : options.at("terms")) {
        const std::string key = term.at("key");
        const auto kind = term.at("kind").get<Kind>();
        if (std::none_of(points.begin(), points.end(), [&](const Point& p) {
                return p.kind == kind && p.holdout;
            })) throw std::runtime_error("missing holdouts");
        if (!model.fits.count(key)) {
            auto residual = [&](const std::vector<double>& values) {
                auto trial = parameters;
                trial[key] = values;
                std::vector<double> errors;
                for (const auto& point : points)
                    if (point.kind == kind && !point.holdout)
                        errors.push_back(predict(point, trial) - point.target);
                return errors;
            };
            model.fits[key] = least_squares(parameters.at(key).get<std::vector<double>>(),
                term.at("bounds").get<std::vector<std::array<double, 2>>>(), residual,
                options.at("iterations"), options.at("tolerance"));
        }
        const auto& fit = model.fits.at(key);
        parameters[key] = fit.values;
        if (!fit.converged || fit.at_bound || !(fit.rmse <= options.at("max_error_eV").get<double>()))
            throw std::runtime_error("interaction fit");
    }
    for (size_t i = 0; i < points.size(); ++i) {
        if (i == model.errors.size())
            model.errors.push_back(predict(points[i], parameters) - points[i].target);
        if (!(std::abs(model.errors[i]) <= options.at("max_error_eV").get<double>()))
            throw std::runtime_error("interaction error");
    }
}
