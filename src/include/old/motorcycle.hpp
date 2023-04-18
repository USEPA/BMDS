#include <string>

#ifndef MOTORCYCLE_H
#define MOTORCYCLE_H

namespace vehicles {

class Motorcycle {

private:

    /// Name
    std::string _name;

public:

    /// Constructor
    __declspec(dllexport) Motorcycle(std::string name);

    /// Get name
    /// @return Name
    __declspec(dllexport) std::string get_name() const;

    /// Ride the bike
    /// @param road Name of the road
    __declspec(dllexport) void ride(std::string road) const;
};

}

#endif
