#pragma once

#include "Request.hpp"
#include <string>
#include <cstddef>
#include <cmath>
#include <utility>
#include "WendlandKernels.hpp"
#include <iostream>
#include <sstream>
#include <iomanip>


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

private:
    double q_(double x, double x0, double h) {
        return abs(x - x0) / h;
    }


    template <
        KernelType kernelType,
        size_t dim,
        size_t cN,
        size_t oN
    >
    double W(double x, double x0, double h) {
        double q = q_(x, x0, h);
        double y = 0;

        if (q > 1) {
            return y;
        }
        if constexpr (dim == 1) {
            if constexpr (cN == 2) {
                if constexpr (oN == 2) {
                    if constexpr (kernelType == WENDLAND) {
                        y = WendlandKernels::WD1C2O2(q, h);
                    }
                }
                if constexpr (oN == 4) {
                    if constexpr (kernelType == WENDLAND) {
                        y = WendlandKernels::WD1C2O4(q, h);
                    }
                }
            }
            if constexpr (cN == 4) {
                if constexpr (oN == 2) {
                    if constexpr (kernelType == WENDLAND) {
                        y = WendlandKernels::WD1C4O2(q, h);
                    }
                }
                if constexpr (oN == 4) {
                    if constexpr (kernelType == WENDLAND) {
                        y = WendlandKernels::WD1C4O4(q, h);
                    }
                }
            }
        }

        return y;
    }


    template <
        KernelType kernelType,
        size_t dim,
        size_t cN,
        size_t oN
    >
    double dW(double x, double x0, double h) {
        double q = q_(x, x0, h);
        double y = 0;

        if (q > 1) {
            return y;
        }
        if constexpr (dim == 1) {
            if constexpr (cN == 2) {
                if constexpr (oN == 2) {
                    if constexpr (kernelType == WENDLAND) {
                        y = WendlandKernels::dWD1C2O2(x, x0, q, h);
                    }
                }
            }
            if constexpr (cN == 4) {
                if constexpr (oN == 2) {
                    if constexpr (kernelType == WENDLAND) {
                        y = WendlandKernels::dWD1C4O2(x, x0, q, h);
                    }
                }
            }
            if constexpr (cN == 4) {
                if constexpr (oN == 2) {
                    if constexpr (kernelType == WENDLAND) {
                        y = WendlandKernels::dWD1C4O2(x, x0, q, h);
                    }
                }
                if constexpr (oN == 4) {
                    if constexpr (kernelType == WENDLAND) {
                        y = WendlandKernels::dWD1C4O4(x, x0, q, h);
                    }
                }
            }
        }

        return y;
    }


    template <
        KernelType kernelType,
        size_t dim,
        size_t cN,
        size_t oN
    >
    double d2W(double x, double x0, double h) {
        double q = q_(x, x0, h);
        double y = 0;

        if (q > 1) {
            return y;
        }

        if constexpr (dim == 1) {
            if constexpr (cN == 2) {
                if constexpr (oN == 2) {
                    if constexpr (kernelType == WENDLAND) {
                        y = WendlandKernels::d2WD1C2O2(q, h);
                    }
                }
            }
            if constexpr (cN == 4) {
                if constexpr (oN == 2) {
                    if constexpr (kernelType == WENDLAND) {
                        y = WendlandKernels::d2WD1C4O2(q, h);
                    }
                }
                if constexpr (oN == 4) {
                    if constexpr (kernelType == WENDLAND) {
                        y = WendlandKernels::d2WD1C4O4(q, h);
                    }
                }
            }
        }

        return y;
    }
public:

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
    void compute() {
        double n = std::ceil(request->nGasParticlesPerUnitSide *
            (request->rightBorder - request->leftBorder));

        double m = request->initRho *
            (request->rightBorder - request->leftBorder) / n;

        size_t intN = static_cast<size_t>(n);

        double step = (request->rightBorder - request->leftBorder) / (n + 1);

        double etaSquared = request->etaSquaredParam * pow(request->h, 2);

        double *rho = new double[intN];
        for (size_t i = 0; i < intN; i++) {
            rho[i] = request->initRho;
        }

        double *xiSPH = new double[intN];
        double *prevXiSPH = new double[intN];

        for (size_t i = 0; i < intN; i++) {
            xiSPH[i] = request->leftBorder + static_cast<double>(i + 1) * step;
        }

        double *prevV = new double[intN];
        double *currV = new double[intN];


        for (size_t i = 0; i < intN; i++) {
            //TODO: maybe it need to be just arbitary callable?
            currV[i] = solutionBurgers(
                xiSPH[i],
                0,
                request->u0,
                request->k,
                request->a
            );
        }

        double nTime = ceil(request->timeMoment / request->tau);
        for (size_t n = 1; n <= nTime; n++) {
            for (size_t i = 0; i < intN; i++) {
                prevXiSPH[i] = xiSPH[i];
                prevV[i] = currV[i];
            }

            for (size_t i = 0; i < intN; i++) {
                rho[i] = 0.0l;
                for (size_t j = 0; j < intN; j++) {
                    rho[i] += W<approximationKernelType, dim, cN, oN>(
                        prevXiSPH[i], prevXiSPH[j], request->h
                    );
                }
                rho[i] *= m;
            }

            if constexpr (computeSchemeNum == 1) {
                double *vFirstDerivatives = new double[intN];
                for (size_t i = 0; i < intN; i++) {
                    vFirstDerivatives[i] = 0.0l;
                    for (size_t j = 0; j < intN; j++) {
                        if constexpr (firstDerivativeApproxSchemeNum == 1) {
                            vFirstDerivatives[i] += prevV[j] *
                                dW<approximationKernelType, dim, cN, oN>(
                                    prevXiSPH[i],
                                    prevXiSPH[j],
                                    request->h
                                ) / rho[j];
                        }
                        if constexpr (firstDerivativeApproxSchemeNum == 2) {
                            vFirstDerivatives[i] += (prevV[j] - prevV[i]) *
                                dW<approximationKernelType, dim, cN, oN>(
                                    prevXiSPH[i],
                                    prevXiSPH[j],
                                    request->h
                                ) / rho[j];
                        }
                        if constexpr (firstDerivativeApproxSchemeNum == 3) {
                            vFirstDerivatives[i] += (prevV[j] + prevV[i]) *
                                dW<approximationKernelType, dim, cN, oN>(
                                    prevXiSPH[i],
                                    prevXiSPH[j],
                                    request->h
                                ) / rho[j];
                        }
                    }
                    vFirstDerivatives[i] *= m;
                }

                for (size_t i = 0; i < intN; i++) {
                    currV[i] = 0.0l;
                    for (size_t j = 0; j < intN; j++) {
                        if constexpr (secondDerivativeApproxSchemeNum == 1) {
                            currV[i] += vFirstDerivatives[j] *
                                dW<approximationKernelType, dim, cN, oN>(
                                    prevXiSPH[i],
                                    prevXiSPH[j],
                                    request->h
                                ) / rho[j];
                        }
                        if constexpr (secondDerivativeApproxSchemeNum == 2) {
                            currV[i] += (vFirstDerivatives[j] - vFirstDerivatives[i]) *
                                dW<approximationKernelType, dim, cN, oN>(
                                    prevXiSPH[i],
                                    prevXiSPH[j],
                                    request->h
                                ) / rho[j];
                        }
                        if constexpr (secondDerivativeApproxSchemeNum == 3) {
                            currV[i] += (vFirstDerivatives[j] + vFirstDerivatives[i]) *
                                dW<approximationKernelType, dim, cN, oN>(
                                    prevXiSPH[i],
                                    prevXiSPH[j],
                                    request->h
                                ) / rho[j];
                        }
                        currV[i] *= request->tau * pow(request->a, 2) * m;
                        currV[i] += prevV[i];
                    }
                }
                delete[] vFirstDerivatives;
            }
            if constexpr (computeSchemeNum == 2) {
                for (size_t i = 0; i < intN; i++) {
                    currV[i] = 0.0l;
                    for (size_t j = 0; j < intN; j++) {
                        if constexpr (otherApproxVariant == 1) {
                            currV[i] += prevV[j] *
                                d2W<approximationKernelType, dim, cN, oN>(
                                    prevXiSPH[i],
                                    prevXiSPH[j],
                                    request->h
                                ) / rho[j];
                        }
                        if constexpr (otherApproxVariant == 2) {
                            currV[i] += (prevV[j] - prevV[i]) *
                                d2W<approximationKernelType, dim, cN, oN>(
                                    prevXiSPH[i],
                                    prevXiSPH[j],
                                    request->h
                                ) / rho[j];
                        }
                        if constexpr (otherApproxVariant == 3) {
                            currV[i] += (prevV[j] + prevV[i]) *
                                d2W<approximationKernelType, dim, cN, oN>(
                                    prevXiSPH[i],
                                    prevXiSPH[j],
                                    request->h
                                ) / rho[j];
                        }
                    }
                    currV[i] *= request->tau * pow(request->a, 2) * m;
                    currV[i] += prevV[i];
                }
            }
            if constexpr (computeSchemeNum == 3) {
                for (size_t i = 0; i < intN; i++) {
                    currV[i] = 0.0l;
                    for (size_t j = 0; j < intN; j++) {
                        if (j != i) {
                            currV[i] += (prevV[i] - prevV[j]) *
                                dW<approximationKernelType, dim, cN, oN>(
                                    prevXiSPH[i],
                                    prevXiSPH[j],
                                    request->h
                                ) / ((prevXiSPH[i] - prevXiSPH[j]) * rho[j]);
                        }
                    }
                    currV[i] *= 2.0l * request->tau * pow(request->a, 2) * m;
                    currV[i] += prevV[i];
                }
            }
            if constexpr (computeSchemeNum == 4) {
                for (size_t i = 0; i < intN; i++) {
                    currV[i] = 0.0l;
                    for (size_t j = 0; j < intN; j++) {
                        currV[i] += (prevV[i] - prevV[j]) *
                            (prevXiSPH[i] - prevXiSPH[j]) *
                            dW<approximationKernelType, dim, cN, oN>(
                                prevXiSPH[i],
                                prevXiSPH[j],
                                request->h
                            ) / ((pow(prevXiSPH[i] - prevXiSPH[j], 2) + etaSquared) * rho[j]);
                    }
                    currV[i] *= 2.0l * request->tau * pow(request->a, 2) * m;
                    currV[i] += prevV[i];
                }
            }

            for (size_t i = 0; i < intN; i++) {
                xiSPH[i] = prevXiSPH[i] + request->tau * currV[i];
            }

            std::string log =
                "n = " +
                std::to_string(n) +
                " (" +
                std::to_string(nTime) +
                ") complete, time = " +
                std::to_string(request->tau * n) +
                " (" +
                std::to_string(request->timeMoment) +
                ")\n";
            std::cout << log;
        }

        double maxErrorV = 0.0l;
        for (size_t i = 0; i < intN; i++) {
            if (0.0l <= xiSPH[i] && xiSPH[i] <= 1.0l) {
                double error = abs(
                    currV[i] - solutionBurgers(
                        xiSPH[i],
                        request->timeMoment,
                        request->u0,
                        request->k,
                        request->a)
                );
                if (error > maxErrorV) {
                    maxErrorV = error;
                }
            }
        }

        double K = pow(request->h, -1);
        double phi = pow(request->nGasParticlesPerUnitSide * request->h, -1);
        computeResult += "\n";
        computeResult +=
            "h = " +
            std::to_string(request->h) +
            ", N = " +
            std::to_string(request->nGasParticlesPerUnitSide) +
            "\n";
        computeResult +=
            "K = " +
            std::to_string(K) +
            ", log(K) = " +
            std::to_string(log10(K)) +
            "\n";
        computeResult +=
            "phi = " +
            std::to_string(phi) +
            ", log(phi) = " +
            std::to_string(log10(phi)) +
            "\n";
        computeResult +=
            "Final maximum speed error (I) = " +
            std::to_string(maxErrorV) +
            ", log(I) = " +
            std::to_string(log10(maxErrorV)) +
            "\n";
        computeResult += "\n";

        constexpr static auto indent = "  ";
        std::stringstream detailedDataBuffer;
        detailedDataBuffer << std::fixed << std::setprecision(16);
        size_t counter = 1;
        detailedDataBuffer << "# " << counter++ <<
            ") xi_SPH_g, " << indent << counter++ <<
            ") an_v, " << counter++ << indent <<
            ") num_v" << "\n" << "#\n";
        double anV = 0.0l;
        for (size_t i = 0; i < intN; i++) {
            anV = solutionBurgers(
                xiSPH[i],
                request->timeMoment,
                request->u0,
                request->k,
                request->a
            );
            detailedDataBuffer <<
                xiSPH[i] <<
                " " <<
                anV <<
                " " <<
                currV[i] <<
                "\n";
        }
        computeResult += detailedDataBuffer.str();
        detailedDataBuffer.clear();

        delete[] prevV;
        delete[] currV;

        delete[] xiSPH;
        delete[] prevXiSPH;
        delete[] rho;

        computeResult += "Finish\n";
    }
};


}