#include <iostream>
#include "Kernel.hpp"
#include "Request.hpp"

using namespace ApproxSPH;

int main() {
    /*
    Request request;
    request.leftBorder = -1.0l;
    request.rightBorder = 2.0l;
    request.nGasParticlesPerUnitSide = 1000;
    request.h = 0.1l;
    request.tau = pow(10.0l, -8);
    request.timeMoment = pow(10.0l, -5);
    request.initRho = 1.0l;
    request.u0 = 0.1l;
    request.k = 2.0l * acos(-1.0l);
    request.a = 0.5l;
    request.etaSquaredParam = 0.001l;
    */

    Request request;
    request.leftBorder = -1.0l;
    request.rightBorder = 2.0l;
    request.nGasParticlesPerUnitSide = 3162;
    request.h = 0.032l;
    request.tau = pow(10.0l, -8);
    request.timeMoment = pow(10.0l, -5);
    request.initRho = 1.0l;
    request.u0 = 0.1l;
    request.k = 2.0l * acos(-1.0l);
    request.a = 0.5l;
    request.etaSquaredParam = 0.001l;

    Kernel kernel(&request);
    kernel.computeThirdScheme<KernelType::WENDLAND, 1, 4, 2>();
    freopen("output.txt", "w", stdout);
    std::cout << kernel.getComputeRes() << std::endl;
    freopen("output_graph.csv", "w", stdout);
    std::cout << kernel.getComputeResForGraph() << std::endl;
    return 0;
}