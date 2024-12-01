#include <iostream>
#include <Eigen/Dense>
#include "libmath/Gaussian.h"
#include "libmath/overlap.h"

int main()
{
    Gaussian G_a(1, Eigen::Vector3d(0, 0, 0), Eigen::Vector3i(0, 0, 0));
    Gaussian G_b(1, Eigen::Vector3d(0, 0, 0), Eigen::Vector3i(0, 0, 0));

    // Compute the overlap integral
    double ov = overlap(G_a, G_b);

    std::cout << G_a.getN() << std::endl;
    // Output the overlap value
    std::cout << ov << std::endl;

    return 0;
}
