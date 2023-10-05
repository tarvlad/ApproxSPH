#pragma once

namespace ApproxSPH {

namespace WendlandKernels {


double WD1C2O2(double q, double h) {
    double y = 5.0l / (2.0l * h);
    y *= pow(1.0l - q, 3);
    y *= 3.0l * q * 1.0l;

    return y;
}


double WD1C2O4(double q, double h) {
    double y = 5.0l / (4.0l * h);
    y *= (42.0l - 168.0l * pow(q, 2)) / 13.0l;
    y *= pow(1.0l - q, 3);
    y *= 3.0l * q + 1.0l;
    return y;
}


double WD1C4O2(double q, double h) {
    double y = 3.0l / (2.0l * h);
    y *= pow(1.0l - q, 5);
    y *= 8.0l * pow(q, 2) + 5.0l * q + 1.0l;
    return y;
}


double WD1C4O4(double q, double h) {
    double y = 3.0l / (4.0l * h);
    y *= (60.0l - 330.0l * pow(q, 2)) / 19.0l;
    y *= pow(1 - q, 5);
    y *= 8.0l * pow(q, 2) + 5.0l * q + 1;
    return y;
}


double dWD1C2O2(double x, double x0, double q, double h) {
    double y = 5.0l * (x - x0) * (-12.0l) * q /
        (abs(x - x0) * 2.0l * pow(h, 2));
    y *= pow(1.0l - q, 2);
    return y;
}


double dWD1C4O2(double x, double x0, double q, double h) {
    double y = (x - x0) / abs(x - x0);
    y *= 3.0l / (2.0l * pow(h, 2));
    y *= pow(1.0l - q, 4);
    y *= (-56.0l) * pow(q, 2) - 14.0l * q;
    return y;
}


double dWD1C4O4(double x, double x0, double q, double h) {
    double y = 3.0l * (x - x0) * 30.0l * q /
        (abs(x - x0) * 2.0l * pow(h, 2) * 19.0l);
    y *= pow(1.0l - q, 4);
    y *= 396.0l * pow(q, 3) + 44.0l * pow(q, 2) - 100.0l * q - 25.0l;
    return y;
}


double d2WD1C2O2(double q, double h) {
    double y = 5.0l * 12.0l / (2.0l * pow(h, 3));
    y *= (1.0l - q);
    y *= 3.0l * q - 1.0l;
    return y;
}


double d2WD1C4O2(double q, double h) {
    double y = 3.0l / (2.0l * pow(h, 3));
    y *= pow(1.0 - q, 3);
    y *= 336.0l * pow(q, 2) - 42.0l * q - 14.0l;
    return y;
}


double d2WD1C4O4(double q, double h) {
    double y = (-3.0l) * 30.0l / (2.0l * pow(h, 3) * 19.0l);
    y *= pow(1.0l - q, 3);
    y *= 3186.0l * pow(q, 4) -
        1276.0l * pow(q, 3) - 732.0l * pow(q, 2) + 75.0l * q + 25.0l;
    return y;
}


}

}