#pragma once
#include "struct.h"

class Algorithm
{
public:
    // Вычисление ускорений
    void Calculation_A(const std::vector<Body>& bodies, std::vector<double>& ax, std::vector<double>& ay, std::vector<double>& az, bool usePNCorrection) const;

    // Рунге-Кутта 2 порядка
    void RungeKutta2(std::vector<Body>& bodies, double dt, bool usePNCorrection);

    // Рунге-Кутта 4 порядка
    void RungeKutta4(std::vector<Body>& bodies, double dt, bool usePNCorrection);

    // LeapFrog
    void Leapfrog(std::vector<Body>& bodies, double dt, bool usePNCorrection);

    // Рунге-Кутта 6 порядка
    void RungeKutta6(std::vector<Body>& bodies, double dt, bool usePNCorrection);

    // Метод Эйлера
    void EulerMethod(std::vector<Body>& bodies, double dt, bool usePNCorrection);

private:
    struct RK4Derivatives {
        double dx, dy, dz;    // производные координат
        double dvx, dvy, dvz; // производные скоростей
    };
};
