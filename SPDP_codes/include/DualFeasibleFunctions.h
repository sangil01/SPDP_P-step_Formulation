#ifndef SPDP_DUAL_FEASIBLE_FUNCTIONS_H
#define SPDP_DUAL_FEASIBLE_FUNCTIONS_H

#include <string>
#include <vector>

#include "GenMultiGraph.h"
#include "ReadData.h"

namespace spdp {

enum class DffFamily {
    Identity,
    FeketeSchepers,
};

struct DffSpec {
    DffFamily family = DffFamily::Identity;
    double parameter = 0.0;
};

const char* dff_family_name(DffFamily family);
std::string dff_spec_name(const DffSpec& spec);

// Validates, sorts, and removes duplicate Fekete--Schepers parameters.
// Lambda=0 is intentionally excluded because it is the identity DFF.
std::vector<double> canonicalize_fs_lambdas(
    const std::vector<double>& lambdas
);

double evaluate_dff(const DffSpec& spec, double normalized_size);

std::vector<double> build_dff_edge_coefficients(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const DffSpec& spec
);

bool same_dff_coefficients(
    const std::vector<double>& lhs,
    const std::vector<double>& rhs
);

}  // namespace spdp

#endif
