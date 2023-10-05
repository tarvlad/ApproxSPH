#include "Kernel.hpp"
#include <cstddef>
#include <cmath>
#include <utility>
#include "WendlandKernels.hpp"

using namespace ApproxSPH;


double q(double x, double x0, double h) {
    return abs(x - x0) / h;
}


template <
    KernelType kernelType,
    size_t dim,
    size_t cN,
    size_t oN
>
double W(double x, double x0, double h) {
    double q = q(x, x0, h);
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
    double q = q(x, x0, h);
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
    double q = q(x, x0, h);
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


template <
    KernelType approximationKernelType,
    size_t firstDerivativeApproxSchemeNum,
    size_t secondDerivativeApproxSchemeNum
>
void Kernel::compute() {
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


    delete[] prevV;
    delete[] currV;

    delete[] xiSPH;
    delete[] prevXiSPH;
    delete[] rho;
}