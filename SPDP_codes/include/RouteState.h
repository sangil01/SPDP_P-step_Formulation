#ifndef SPDP_ROUTE_STATE_H
#define SPDP_ROUTE_STATE_H

#include <vector>

#include "GenMultiGraph.h"

namespace spdp {

struct OnboardSkip {
    int request_index = -1;
    int container_type = -1;
    int treatment_location = -1;
    bool is_full = true;
};

class RouteState {
public:
    RouteState() = default;
    explicit RouteState(std::vector<OnboardSkip> onboard);

    bool try_pickup(OnboardSkip skip);
    bool try_empty_request(int request_index, int treatment_location);
    int empty_at_treatment(int treatment_location);
    bool try_delivery(int container_type);

    State canonical_state() const;
    bool is_empty() const;
    const std::vector<OnboardSkip>& onboard() const;

private:
    std::vector<OnboardSkip> onboard_;
};

}  // namespace spdp

#endif
