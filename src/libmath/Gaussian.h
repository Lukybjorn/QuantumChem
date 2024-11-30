#pragma once

#include <Eigen/Dense>
#include <cmath>

class Gaussian {
private:
    double _zeta;  // Width parameter (exponent coefficient)
    Eigen::Vector3f _center;  // Center of the Gaussian
    Eigen::Vector3i _angular; // Angular momentum quantum numbers

public:
    // Constructor
    Gaussian(double zeta, const Eigen::Vector3f& center, const Eigen::Vector3i& angular)
        : _zeta(zeta), _center(center), _angular(angular) {}

    // Default Constructor
    Gaussian() : _zeta(0.0), _center(Eigen::Vector3f::Zero()), _angular(Eigen::Vector3i::Zero()) {}

    // Getters
    double getZeta() const { return _zeta; }
    Eigen::Vector3f getCenter() const { return _center; }
    Eigen::Vector3i getAngular() const { return _angular; }

    // Setters
    void setZeta(double zeta) { _zeta = zeta; }
    void setCenter(const Eigen::Vector3f& center) { _center = center; }
    void setAngular(const Eigen::Vector3i& angular) { _angular = angular; }

    // Method to evaluate the Gaussian function at a given point
    double evaluate(const Eigen::Vector3f& point) const {
        Eigen::Vector3f diff = point - _center;
        double r2 = diff.squaredNorm();
        return std::exp(-_zeta * r2);
    }
};

