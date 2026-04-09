#pragma once

#include <string>

namespace platform {
namespace common {

enum class ContactAlgorithm {
    NativeMesh,
    SdfFirstOrder,
};

inline const char* ContactAlgorithmToCliName(ContactAlgorithm algorithm) {
    switch (algorithm) {
        case ContactAlgorithm::NativeMesh:
            return "mesh";
        case ContactAlgorithm::SdfFirstOrder:
            return "sdf_1st";
        default:
            return "unknown";
    }
}

inline bool ParseContactAlgorithm(const std::string& value, ContactAlgorithm& out_algorithm) {
    if (value == "mesh" || value == "native_mesh") {
        out_algorithm = ContactAlgorithm::NativeMesh;
        return true;
    }
    if (value == "sdf_1st" || value == "sdf1" || value == "first_order_sdf") {
        out_algorithm = ContactAlgorithm::SdfFirstOrder;
        return true;
    }
    if (value == "sdf_2nd" || value == "sdf2" || value == "hessian" || value == "hessian_sdf") {
        out_algorithm = ContactAlgorithm::SdfFirstOrder;
        return true;
    }
    return false;
}

}  // namespace common
}  // namespace platform
