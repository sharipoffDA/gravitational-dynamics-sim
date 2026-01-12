#pragma once
#include"struct.h"
#include"visualisation.h"
#include"algorithm.h"

class BodySystem {

private:
    std::vector<Body> bodies;
    std::vector<double> time_values;
    std::vector<double> energy_values;
    std::vector<double> momentum_values;

    // Определение устойчивости системы
    class OrbitState {
    public:
        double time; // Время
        double distance; // Расстояние между телами
        double eccentricity; // Эксцентриситет
    };

    // Параметры орбит
    class OrbitParameters {
    public:
        double perihelion;  // Перицентр
        double aphelion;    // Апоцентр
        double eccentricity;  // Эксцентриситет
        double semimajor_axis; // Большая полуось
    };

    std::vector<OrbitParameters> orbit_params;
    // Вычисление импульса
    double Momentum() const;

    // Вычисление кинетической энергии
    double KineticEnergy() const;

    // Вычисление потенциальной энергии
    double PotentialEnergy() const;

    // Подробная информация про то, насколько отклоняется орбита на протяжении всего времени
    void SaveToFileOrbitHistory(const std::vector<std::vector<OrbitState>>& orbit_history);

    // Центрирование системы после каждого оборота
    void recenter();

    // Запись параметров орбиты в файл
    void SaveOrbitParameters(const std::string& filename) const;

    Algorithm algorithm;

public:

    BodySystem(const std::vector<Body>& initial_bodies);

    // Запись данных в файл и расчёт с помощью метода Leapfrog
    void SaveToFile_Leapfrog(const std::string& filename, double t_end, double dt, bool usePNCorrection);

    // Запись данных в файл и расчёт с помощью метода Рунге-Кутта 2 порядка
    void SaveToFile_RK2(const std::string& filename, double t_end, double dt, bool usePNCorrection);

    // Запись данных в файл и расчёт с помощью метода Рунге-Кутта 4 порядка
    void SaveToFile_RK4(const std::string& filename, double t_end, double dt, bool usePNCorrection);

    // Запись данных в файл и расчёт с помощью метода Рунге-Кутта 6 порядка
    void SaveToFile_RK6(const std::string& filename, double t_end, double dt, bool usePNCorrection);

    // Запись данных в файл и расчёт с помощью метода Эйлера
    void SaveToFile_Euler(const std::string& filename, double t_end, double dt, bool usePNCorrection);

    Visual visual;
};