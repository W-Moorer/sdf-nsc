#pragma once
#include_next <tbb/task.h>
#include <tbb/task_group.h>

namespace tbb {
namespace task {

struct openvdb_task_self_compat_proxy {
    void cancel_group_execution() const {
        if (auto* ctx = current_context()) {
            ctx->cancel_group_execution();
        }
    }
};

inline openvdb_task_self_compat_proxy self() {
    return openvdb_task_self_compat_proxy{};
}

}  // namespace task
}  // namespace tbb
