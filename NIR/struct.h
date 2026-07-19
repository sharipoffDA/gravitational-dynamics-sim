#pragma once
#define _USE_MATH_DEFINES

#ifndef NOMINMAX
#define NOMINMAX
#endif

#include <algorithm>
#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<cmath>
#include <filesystem>
#include <numbers>
#include <QString>
#include <QColor>

#ifdef _WIN32
    #include <windows.h>
#else
    #include <unistd.h>
    #include <limits.h>
#endif

namespace fs = std::filesystem;

fs::path get_executable_dir();

const double c = 299792458.0;
const double G = 6.67430e-11;
const double AU = 149597870700.0; // Астрономическая единица
const double M_SUN = 1.989e30;    // Масса Солнца (кг)
const double M_EARTH = 5.972e24;

// Параметры звезды TOI-270
const double M_STAR_TOI270 = 0.386 * M_SUN;

const double DEG_TO_RAD = M_PI / 180.0;
const double RAD_TO_DEG = 180.0 / M_PI;
const double day_to_sec = 86400.0;
const double year_to_sec = 365 * day_to_sec;

struct KeplerianElements {
    std::string name;
    double a_au;        // Большая полуось (а.е.)
    double incl_deg;    // Наклон орбиты (градусы)
    double omega_deg;   // Аргумент периастра (градусы)
    double Omega_deg;   // Долгота восходящего узла (градусы)
    double mass;        // Масса (в массах Земли)
    double ecc;         // Эксцентриситет
    double period_days; // Период обращения (дни)
    double T_p_year;    // Время прохождения периастра (год)
};

class Body {
public:
    double m;          // Масса (кг)
    double x, y, z;    // Координаты (м)
    double vx, vy, vz; // Скорости (м/с)
};

