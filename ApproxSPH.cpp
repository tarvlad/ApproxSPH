#include <iostream>
#include "Kernel.hpp"
#include "Request.hpp"

using namespace ApproxSPH;

int main() {
    Request request;
    request.leftBorder = -1.0l;
    request.rightBorder = 2.0l;
    request.nGasParticlesPerUnitSide = 100;
    request.h = 0.1l;
    request.tau = pow(10.0l, -8);
    request.timeMoment = 100.0l * request.tau;
    request.initRho = 1;
    request.u0 = 0.1l;
    request.k = 2.0l * acos(-1.0l);
    request.a = 0.5l;
    request.etaSquaredParam = 0.001l;

    Kernel kernel(&request);
    kernel.compute<WENDLAND, 1, 2, 2, 1, 1, 1, 0>();
    freopen("output.txt", "w", stdout);
    std::cout << kernel.getComputeRes() << std::endl;

    return 0;
}