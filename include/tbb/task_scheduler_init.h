#pragma once

#include <tbb/info.h>

namespace tbb {

// Compatibility shim for older APIs that still include task_scheduler_init.h.
class task_scheduler_init {
  public:
    static const int automatic = -1;

    explicit task_scheduler_init(int /*threads*/ = automatic) {}
    ~task_scheduler_init() = default;

    void initialize(int /*threads*/ = automatic) {}
    void terminate() {}
    bool is_active() const { return true; }

    static int default_num_threads() {
        const auto n = tbb::info::default_concurrency();
        return n > 0 ? static_cast<int>(n) : 1;
    }
};

}  // namespace tbb
