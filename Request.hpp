#pragma once

#include <cstddef>


namespace ApproxSPH {


struct Request {
    double leftBorder;
    double rightBorder;
    double h;
    double tau;
    double timeMoment;
    double initRho;
    double u0;
    double k;
    double a;
    double etaSquaredParam;

    size_t nGasParticlesPerUnitSide;
};


}