#include "RouteState.h"

#include <algorithm>
#include <stdexcept>
#include <utility>

namespace spdp {

RouteState::RouteState(std::vector<OnboardSkip> onboard)
    : onboard_(std::move(onboard)) {
    if (onboard_.size() > 2U) {
        throw std::runtime_error("Route state exceeds vehicle capacity.");
    }
}

bool RouteState::try_pickup(OnboardSkip skip) {
    if (!skip.is_full || onboard_.size() >= 2U) {
        return false;
    }
    onboard_.push_back(skip);
    return true;
}

int RouteState::empty_at_treatment(int treatment_location) {
    int emptied_count = 0;
    for (OnboardSkip& skip : onboard_) {
        if (skip.is_full && skip.treatment_location == treatment_location) {
            skip.is_full = false;
            ++emptied_count;
        }
    }
    return emptied_count;
}

bool RouteState::try_delivery(int container_type) {
    auto matching = onboard_.end();
    for (auto it = onboard_.begin(); it != onboard_.end(); ++it) {
        if (!it->is_full && it->container_type == container_type) {
            matching = it;
        }
    }
    if (matching == onboard_.end()) {
        return false;
    }
    onboard_.erase(matching);
    return true;
}

State RouteState::canonical_state() const {
    if (onboard_.size() > 2U) {
        throw std::runtime_error("Route state exceeds vehicle capacity.");
    }

    std::vector<StateToken> tokens;
    tokens.reserve(2);
    for (const OnboardSkip& skip : onboard_) {
        if (skip.is_full) {
            tokens.push_back(StateToken{
                'F', skip.container_type, skip.treatment_location
            });
        } else {
            tokens.push_back(StateToken{'E', skip.container_type, -1});
        }
    }
    while (tokens.size() < 2U) {
        tokens.push_back(StateToken{'N', -1, -1});
    }
    State state{tokens[0], tokens[1]};
    if (state[1] < state[0]) {
        std::swap(state[0], state[1]);
    }
    return state;
}

bool RouteState::is_empty() const {
    return onboard_.empty();
}

const std::vector<OnboardSkip>& RouteState::onboard() const {
    return onboard_;
}

}  // namespace spdp
