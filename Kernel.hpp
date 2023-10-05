#pragma once

#include "Request.hpp"
#include "Output.hpp"
#include <cmath>


namespace ApproxSPH {


enum KernelType {
    WENDLAND
};


/// Kernel DO NOT control request and output lifetime
class Kernel {
    Request *request;
    Output *output;

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
    Kernel(Request *request, Output *output) :
            request(request), output(output) {
    }

    template <
        KernelType approximationKernelType,
        size_t firstDerivativeApproxSchemeNum,
        size_t secondDerivativeApproxSchemeNum
    >
    void compute();
};


}