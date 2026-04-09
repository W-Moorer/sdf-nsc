#pragma once
#include <tbb/task_arena.h>

namespace tbb {

class task_scheduler_init {
  public:
    static const int automatic = -1;

    explicit task_scheduler_init(int = automatic) {}
    ~task_scheduler_init() = default;

    void initialize(int = automatic) {}
    void terminate() {}
    bool is_active() const { return true; }

    static int default_num_threads() {
        return task_arena::automatic;
    }
};

}  // namespace tbb
