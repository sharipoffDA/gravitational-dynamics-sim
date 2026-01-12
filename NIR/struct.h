#pragma once
#define _USE_MATH_DEFINES
#include"matplotlibcpp.h"
#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<cmath>

namespace plt = matplotlibcpp;

const double c = 299792458.0;
const double G = 6.67430e-11;
const double AU = 149597870700.0; // Астрономическая единица
const double M_SUN = 1.989e30;     // Масса Солнца (кг)
const double M_EARTH = 5.972e24;

// Параметры звезды TOI-270
const double M_STAR_TOI270 = 0.386 * M_SUN;

const double DEG_TO_RAD = M_PI / 180.0;
const double RAD_TO_DEG = 180.0 / M_PI;
const double year_to_sec = 365 * 86400.0;

class Body {
public:
    double m;  // Масса (кг)
    double x, y, z;  // Координаты (м)
    double vx, vy, vz;  // Скорости (м/с)
};

