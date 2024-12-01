#pragma once

#include <Eigen/Dense>
#include <utility>

class Atom {
private:
    double _mass;               // Atomic mass
    double _charge;             // Atomic charge
    Eigen::Vector3d _coords;    // Atomic coordinates in 3D space

public:
    // Constructor
    Atom(double mass, double charge, Eigen::Vector3d  coords)
        : _mass(mass), _charge(charge), _coords(std::move(coords)) {}

    // Default Constructor
    Atom() : _mass(0.0), _charge(0.0), _coords(Eigen::Vector3d::Zero()) {}

    // Getters
    [[nodiscard]] double getMass() const { return _mass; }
    [[nodiscard]] double getCharge() const { return _charge; }
    [[nodiscard]] Eigen::Vector3d getCoords() const { return _coords; }

    // Setters
    void setMass(double mass) { _mass = mass; }
    void setCharge(double charge) { _charge = charge; }
    void setCoords(const Eigen::Vector3d& coords) { _coords = coords; }
};