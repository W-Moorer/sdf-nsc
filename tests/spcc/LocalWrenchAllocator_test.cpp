#include "platform/backend/spcc/LocalWrenchAllocator.h"

#include <gtest/gtest.h>

namespace {

using platform::backend::spcc::LocalWrenchAllocator;
using platform::backend::spcc::ReferenceWrench;
using platform::backend::spcc::SupportWrenchPoint;
using platform::backend::spcc::WrenchAllocationResult;

TEST(LocalWrenchAllocatorTest, MatchesForceAndMomentWithNonnegativeLoads) {
    std::vector<SupportWrenchPoint> supports(2);
    supports[0].x_W = chrono::ChVector3d(-0.5, 0.0, 0.0);
    supports[0].n_W = chrono::ChVector3d(0.0, 1.0, 0.0);
    supports[0].initial_load = 1.0;
    supports[1].x_W = chrono::ChVector3d(0.5, 0.0, 0.0);
    supports[1].n_W = chrono::ChVector3d(0.0, 1.0, 0.0);
    supports[1].initial_load = 1.0;

    ReferenceWrench reference;
    reference.origin_W = chrono::ChVector3d(0.0, 0.0, 0.0);
    reference.force_W = chrono::ChVector3d(0.0, 4.0, 0.0);
    reference.moment_W = chrono::ChVector3d(0.0, 0.0, -1.0);
    reference.total_load = 4.0;

    WrenchAllocationResult result;
    LocalWrenchAllocator::Allocate(supports, reference, result);

    ASSERT_TRUE(result.feasible);
    ASSERT_EQ(result.loads.size(), 2u);
    EXPECT_NEAR(result.loads[0], 3.0, 1.0e-6);
    EXPECT_NEAR(result.loads[1], 1.0, 1.0e-6);
    EXPECT_LT(result.force_residual, 1.0e-8);
    EXPECT_LT(result.moment_residual, 1.0e-8);
}

}  // namespace
