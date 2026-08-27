#include "TestSupport.h"

#include "RouteState.h"

SPDP_TEST(route_state) {
    spdp::RouteState state;

    SPDP_CHECK(state.try_pickup(spdp::OnboardSkip{0, 2, 7, true}));
    SPDP_CHECK(!state.try_delivery(2));
    SPDP_CHECK_EQ(state.empty_at_treatment(7), 1);
    SPDP_CHECK(state.try_delivery(2));
    SPDP_CHECK(state.is_empty());

    SPDP_CHECK(state.try_pickup(spdp::OnboardSkip{0, 1, 3, true}));
    SPDP_CHECK(state.try_pickup(spdp::OnboardSkip{1, 1, 4, true}));
    SPDP_CHECK(!state.try_pickup(spdp::OnboardSkip{2, 1, 5, true}));
    SPDP_CHECK_EQ(state.onboard().size(), 2U);

    spdp::RouteState restored({
        spdp::OnboardSkip{3, 4, 8, true},
        spdp::OnboardSkip{4, 4, 9, false},
    });
    SPDP_CHECK_EQ(restored.onboard().size(), 2U);
    SPDP_CHECK_EQ(restored.empty_at_treatment(8), 1);
    SPDP_CHECK(restored.try_delivery(4));
}
