#pragma once
#include "struct.h"
#include <iostream>
#include <string>
#include <iomanip>
#include <sstream>
#include <numbers>

class QCustomPlot; // Предварительное объявление

// Структура для передачи всех настроек интерфейса
struct PlotSettings {
    QString title;
    int projection; // 0: XY, 1: XZ, 2: YZ
    int theme;      // 0: Светлая, 1: Темная (Космос)
    int fontSize;
    bool showAxisLabels;
    bool showLegend;
    bool showGrid;

    bool keepScale = false;
    
    struct BodyStyle {
        bool visible;
        QString name;
        QColor color;
    };
    std::vector<BodyStyle> bodyStyles;
};



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
    void PlotOrbits(QCustomPlot* plot, const std::string& filename, const PlotSettings& settings) const;

    // Метод настройки стиля отрисовки
    static void SetupGWPlotStyle(QCustomPlot* plot, bool darkTheme);

    // Визуализация динамики слияния компактных объектов
    void PlotGWDynamics(QCustomPlot* plot, const QVector<double>& x1, const QVector<double>& y1, const QVector<double>& x2, const QVector<double>& y2, int tailSize, const QString& name1 = "Объект 1", const QString& name2 = "Объект 2", int fontSize = 10);

    // Визуализация гравитационных волн от слияния компактных объектов

    //void PlotGWWaveform(QCustomPlot* plot, const QVector<double>& t, const QVector<double>& hp, const QVector<double>& hx, int fontSize, const QVector<double>& ligo_t = QVector<double>(), const QVector<double>& ligo_h = QVector<double>());
    void PlotGWWaveform(QCustomPlot* plot, const QVector<double>& t, const QVector<double>& hp, const QVector<double>& hx, int fontSize = 10);

private:
    //static QString lastName;
    std::vector<Body> bodies;
    std::vector<double> time_values;
    std::vector<double> energy_values;
    std::vector<double> momentum_values;
};