#pragma once

namespace platform {
namespace backend {

class IRigidSystem {
public:
    virtual ~IRigidSystem() = default;
    
    // Initialize the system
    virtual void Initialize() = 0;
    
    // Step simulation
    virtual void StepDynamics(double step_size) = 0;
    
    // System time
    virtual double GetTime() const = 0;
    
    // Num Contacts
    virtual unsigned int GetNumContacts() const = 0;
};

} // namespace backend
} // namespace platform
