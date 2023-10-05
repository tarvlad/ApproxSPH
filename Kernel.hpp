#pragma once

#include "Request.hpp"
#include <cmath>
#include <string>


namespace ApproxSPH {


enum KernelType {
    WENDLAND
};


/// Kernel DO NOT control request lifetime
class Kernel {
    Request *request;
    std::string computeResult;

    double solutionBurgers(
        double x,
        double t,
        double u0,
        double k,
        double a
    ) {
        return u0 * cos(k * x) + exp(-pow(a * k, 2) * t);
    }

public:
    Kernel(Request *request) :
            request(request) {
    }

    const std::string &getComputeRes() const {
        return computeResult;
    }

    template <
        KernelType approximationKernelType,
        size_t dim,
        size_t cN,
        size_t oN,
        size_t computeSchemeNum,
        size_t firstDerivativeApproxSchemeNum,
        size_t secondDerivativeApproxSchemeNum,
        size_t otherApproxVariant
    >
    void compute();
};


}