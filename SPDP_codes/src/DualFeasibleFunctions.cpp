#include "DualFeasibleFunctions.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace spdp {
namespace {

constexpr double kDomainTolerance = 1e-9;
constexpr double kCoefficientEqualityTolerance = 1e-12;

void validate_fs_lambda(double lambda) {
    if (!std::isfinite(lambda) || lambda <= 0.0 || lambda >= 0.5) {
        throw std::runtime_error(
            "Fekete--Schepers lambda must be finite and satisfy 0 < lambda < 0.5; "
            "lambda=0 is represented by the identity DFF."
        );
    }
}

double normalized_unit_interval(double value) {
    if (!std::isfinite(value) || value < -kDomainTolerance ||
        value > 1.0 + kDomainTolerance) {
        throw std::runtime_error(
            "A DFF edge size must be finite and lie in [0,1]."
        );
    }
    return std::clamp(value, 0.0, 1.0);
}

}  // namespace

const char* dff_family_name(DffFamily family) {
    switch (family) {
        case DffFamily::Identity:
            return "identity";
        case DffFamily::FeketeSchepers:
            return "fekete-schepers";
    }
    throw std::runtime_error("Unsupported DFF family.");
}

std::string dff_spec_name(const DffSpec& spec) {
    if (spec.family == DffFamily::Identity) {
        return dff_family_name(spec.family);
    }
    validate_fs_lambda(spec.parameter);
    std::ostringstream out;
    out << dff_family_name(spec.family) << "-lambda-"
        << std::setprecision(std::numeric_limits<double>::max_digits10)
        << spec.parameter;
    return out.str();
}

std::vector<double> canonicalize_fs_lambdas(
    const std::vector<double>& lambdas
) {
    std::vector<double> result = lambdas;
    for (double lambda : result) {
        validate_fs_lambda(lambda);
    }
    std::sort(result.begin(), result.end());
    result.erase(
        std::unique(
            result.begin(),
            result.end(),
            [](double lhs, double rhs) {
                return std::abs(lhs - rhs) <= kCoefficientEqualityTolerance;
            }
        ),
        result.end()
    );
    return result;
}

double evaluate_dff(const DffSpec& spec, double normalized_size) {
    const double x = normalized_unit_interval(normalized_size);
    switch (spec.family) {
        case DffFamily::Identity:
            return x;
        case DffFamily::FeketeSchepers:
            validate_fs_lambda(spec.parameter);
            if (x <= spec.parameter) {
                return 0.0;
            }
            if (x < 1.0 - spec.parameter) {
                return x;
            }
            return 1.0;
    }
    throw std::runtime_error("Unsupported DFF family.");
}

std::vector<double> build_dff_edge_coefficients(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const DffSpec& spec
) {
    if (!std::isfinite(data.time_limit) || data.time_limit <= 0.0) {
        throw std::runtime_error(
            "A positive finite route time limit is required for DFF coefficients."
        );
    }
    std::vector<double> coefficients;
    coefficients.reserve(graph.number_of_edges());
    for (const EdgeRecord& edge : graph.edges()) {
        coefficients.push_back(
            evaluate_dff(spec, edge.data.time / data.time_limit)
        );
    }
    return coefficients;
}

bool same_dff_coefficients(
    const std::vector<double>& lhs,
    const std::vector<double>& rhs
) {
    if (lhs.size() != rhs.size()) {
        return false;
    }
    for (std::size_t index = 0; index < lhs.size(); ++index) {
        if (std::abs(lhs[index] - rhs[index]) >
            kCoefficientEqualityTolerance) {
            return false;
        }
    }
    return true;
}

}  // namespace spdp
