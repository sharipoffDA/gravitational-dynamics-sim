#pragma once
#include "struct.h"
#include <iostream>
#include <string>
#include <iomanip>
#include <sstream>

class Visual
{
public:

    // Конструктор
    Visual(std::vector<Body>& bodies_, std::vector<double>& time_values_, std::vector<double>& energy_values_, std::vector<double>& momentum_values_);

    // Визуализация закона сохранения энергии
    void PlotEnergy() const;

    // Визуализация закона сохранения импульса
    void PlotMomentum() const;

    // Визуализация орбит
    void PlotOrbits(const std::string& filename, const std::vector<std::string>& labels, const std::vector<std::string>& colors) const;

private:
    std::vector<Body> bodies;
    std::vector<double> time_values;
    std::vector<double> energy_values;
    std::vector<double> momentum_values;
};