#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include <numbers>
#include <iostream>
#include "overlap.h"

double s_iterative(int i, int a, int b, const Eigen::Vector3d& RA, const Eigen::Vector3d& RB, double alpha, double beta) {
    double gamma = alpha + beta;
    double P_i = (alpha * RA[i] + beta * RB[i]) / gamma;
    double AB_i = RA[i] - RB[i];

    std::vector<std::vector<double>> s_table(a + 1, std::vector<double>(b + 1, 0.0));

    s_table[0][0] = 1.0;

    for (int a_idx = 0; a_idx <= a; ++a_idx) {
        for (int b_idx = 0; b_idx <= b; ++b_idx) {
            if (a_idx == 0 && b_idx == 0) continue;

            if (b_idx == 0) {
                s_table[a_idx][b_idx] = -(RA[i] - P_i) * (a_idx > 0 ? s_table[a_idx - 1][b_idx] : 0)
                                        + (a_idx > 1 ? (a_idx - 1) / (2.0 * gamma) * s_table[a_idx - 2][b_idx] : 0);
            } else {
                s_table[a_idx][b_idx] = s_table[a_idx][b_idx - 1]
                                        * AB_i + (a_idx > 0 ? s_table[a_idx - 1][b_idx - 1] : 0);
            }
        }
    }

    return s_table[a][b];
}

Eigen::Vector3d s(const Gaussian& A, const Gaussian& B) {
    double alpha = A.getZeta();
    double beta = B.getZeta();
    const Eigen::Vector3d& R_A = A.getCenter();
    const Eigen::Vector3d& R_B = B.getCenter();
    const Eigen::Vector3i& L_A = A.getAngular();
    const Eigen::Vector3i& L_B = B.getAngular();
    Eigen::Vector3d result;

    for (int i = 0; i < 3; ++i) {
        result[i] = s_iterative(i, L_A[i], L_B[i], R_A, R_B, alpha, beta);
    }

    return result;
}

double overlap(const Gaussian& A, const Gaussian& B) {
    double alpha = A.getZeta();
    double beta = B.getZeta(); // Corrected from A.getZeta() to B.getZeta()

    const Eigen::Vector3d& R_A = A.getCenter();
    const Eigen::Vector3d& R_B = B.getCenter();

    if (R_A == R_B && A.getAngular() != B.getAngular()){
        return 0;
    }

    Eigen::Vector3d s_vec = s(A, B);

    double E_AB = std::exp(-(alpha * beta / (alpha + beta)) *
                          (R_A - R_B).squaredNorm());

    double OV = A.getN() * B.getN() * E_AB * std::pow(std::numbers::pi / (alpha + beta), 1.5)
                     * s_vec.prod();  // Product of components

    return OV;
}
