#pragma once

#include <Eigen/Dense>
#include <cmath>
#include <utility>

class Gaussian {
private:
    double _zeta;  // Width parameter (exponent coefficient)
    Eigen::Vector3d _center;  // Center of the Gaussian
    Eigen::Vector3i _angular; // Angular momentum quantum numbers
    double _N;

public:
    // Constructor
    Gaussian(double zeta, Eigen::Vector3d  center, Eigen::Vector3i  angular)
        : _zeta(zeta), _center(std::move(center)), _angular(std::move(angular)) {
        this->_N = std::pow(zeta * 2 / std::numbers::pi, 3.0 / 4.0);;
    }

    // Default Constructor
    Gaussian() : _zeta(0.0), _center(Eigen::Vector3d::Zero()), _angular(Eigen::Vector3i::Zero()) {}

    // Getters
    [[nodiscard]] double getZeta() const { return _zeta; }
    [[nodiscard]] Eigen::Vector3d getCenter() const { return _center; }
    [[nodiscard]] Eigen::Vector3i getAngular() const { return _angular; }
    [[nodiscard]] double getN() const { return _N; }

    // Setters
    void setZeta(double zeta) { _zeta = zeta; }
    void setCenter(const Eigen::Vector3d& center) { _center = center; }
    void setAngular(const Eigen::Vector3i& angular) { _angular = angular; }
    void setN(double N) { _N = N; }

    // Method to evaluate the Gaussian function at a given point
    [[nodiscard]] double evaluate(const Eigen::Vector3d& point) const {
        // Ensure both point and center are of the same type (Eigen::Vector3d)
        Eigen::Vector3d diff = point - _center;  // No need to cast center since it's already Eigen::Vector3d
        double r2 = diff.squaredNorm();
        return std::exp(-_zeta * r2);
    }
};
