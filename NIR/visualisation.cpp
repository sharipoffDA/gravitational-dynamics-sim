#include "visualisation.h"
#include "qcustomplot.h"
#include <QPen>
#include <QColor>
#include <QVector>
#include <cmath>
#include <algorithm>

Visual::Visual(std::vector<Body>& bodies_, std::vector<double>& time_values_, std::vector<double>& energy_values_, std::vector<double>& momentum_values_) :
    bodies(bodies_), time_values(time_values_), energy_values(energy_values_), momentum_values(momentum_values_) {}

std::string formatTime(double time) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(1) << std::round(time * 365 * 10) / 10;
    return oss.str();
}

//Визуализация орбит
void Visual::PlotOrbits(QCustomPlot* plot, const std::string& filename, const PlotSettings& settings) const {
    if (!plot) return;

    std::ifstream inp(filename);
    if (!inp.is_open()) {
        std::cerr << "Ошибка: не удалось открыть файл " << filename << std::endl;
        return;
    }

    std::vector<double> time;
    std::vector<QVector<double>> x, y, z;
    
    std::string line;
    int bodyCount = 0;
    int lineCounter = 0;
    std::getline(inp, line);

    // Считывание формата: Время(t) X Y Z(тело1) X Y Z(тело2) ...
    const int stepSkip = 1; 

    while (std::getline(inp, line)) {
        if (line.empty()) continue;
        
        lineCounter++;
        if (lineCounter % stepSkip != 0) continue; // Пропуск промежуточных точек

        std::istringstream iss(line);
        double t;
        if (!(iss >> t)) continue;
        
        time.push_back(t);

        double cx, cy, cz;
        int curBody = 0;
        
        while (iss >> cx >> cy >> cz) {
            if (curBody >= x.size()) {
                x.emplace_back(); y.emplace_back(); z.emplace_back();
                bodyCount++;
            }
            x[curBody].push_back(cx / AU); 
            y[curBody].push_back(cy / AU);
            z[curBody].push_back(cz / AU);
            curBody++;
        }
    }
    inp.close();

    plot->clearGraphs();
    
    // Очистка старого заголовка, если он был
    if (plot->plotLayout()->rowCount() > 1) {
        plot->plotLayout()->remove(0);
        plot->plotLayout()->simplify();
    }

    // Заголовок
    QCPLayoutElement* oldTitle = plot->plotLayout()->element(0, 0);
    if (oldTitle && dynamic_cast<QCPTextElement*>(oldTitle)) {
        plot->plotLayout()->take(oldTitle); // Удаление из слоя
        delete oldTitle;                    // Освобождение памяти
        plot->plotLayout()->simplify();     // Удаление пустой строки/ячейки
    }

    // Создание нового, если нужно
    if (!settings.title.isEmpty()) {
        plot->plotLayout()->insertRow(0); // Вставляем новую строку для заголовка
        QCPTextElement *title = new QCPTextElement(plot, settings.title, QFont("sans", settings.fontSize + 2, QFont::Bold));
        
        if (settings.theme == 1) title->setTextColor(Qt::white);
        else title->setTextColor(Qt::black);
        
        plot->plotLayout()->addElement(0, 0, title);
    }

    QString labelH, labelV;
    plot->clearPlottables();

    for (int i = 0; i < bodyCount; ++i) {
        bool visible = true;
        QString name = QString("Тело %1").arg(i + 1);
        QColor color = QColor::fromHsv((i * 50) % 360, 200, 200);

        if (i < settings.bodyStyles.size()) {
            visible = settings.bodyStyles[i].visible;
            name = settings.bodyStyles[i].name;
            color = settings.bodyStyles[i].color;
        }

        if (!visible) continue;

        QCPCurve *orbitCurve = new QCPCurve(plot->xAxis, plot->yAxis);
        
        orbitCurve->setBrush(Qt::NoBrush);
        
        // Настройка стиля линии
        orbitCurve->setPen(QPen(color, 2));
        orbitCurve->setName(name);

        // Выбор данных в зависимости от проекции
        if (settings.projection == 0) { // XY
            orbitCurve->setData(x[i], y[i]);
            labelH = "X (a.e.)"; labelV = "Y (a.e.)";
        } else if (settings.projection == 1) { // XZ
            orbitCurve->setData(x[i], z[i]);
            labelH = "X (a.e.)"; labelV = "Z (a.e.)";
        } else { // YZ
            orbitCurve->setData(y[i], z[i]);
            labelH = "Y (a.e.)"; labelV = "Z (a.e.)";
        }
    }

    // Шрифты и подписи осей
    QFont axisFont = plot->xAxis->tickLabelFont();
    axisFont.setPointSize(settings.fontSize);
    plot->xAxis->setTickLabelFont(axisFont);
    plot->yAxis->setTickLabelFont(axisFont);

    if (settings.showAxisLabels) {
        QFont labelFont = plot->xAxis->labelFont();
        labelFont.setPointSize(settings.fontSize + 2);
        plot->xAxis->setLabelFont(labelFont);
        plot->yAxis->setLabelFont(labelFont);
        plot->xAxis->setLabel(labelH);
        plot->yAxis->setLabel(labelV);
    } else {
        plot->xAxis->setLabel("");
        plot->yAxis->setLabel("");
    }

    // Сетка
    plot->xAxis->grid()->setVisible(settings.showGrid);
    plot->yAxis->grid()->setVisible(settings.showGrid);

    // Легенда
    plot->legend->setVisible(settings.showLegend);

    // Тема (Светлая / Темная)
    if (settings.theme == 1) { // Темная
        plot->setBackground(QBrush(QColor(15, 15, 25))); // Темно-синий фон
        plot->xAxis->setBasePen(QPen(Qt::white));
        plot->yAxis->setBasePen(QPen(Qt::white));
        plot->xAxis->setTickPen(QPen(Qt::white));
        plot->yAxis->setTickPen(QPen(Qt::white));
        plot->xAxis->setSubTickPen(QPen(Qt::white));
        plot->yAxis->setSubTickPen(QPen(Qt::white));
        plot->xAxis->setTickLabelColor(Qt::white);
        plot->yAxis->setTickLabelColor(Qt::white);
        plot->xAxis->setLabelColor(Qt::white);
        plot->yAxis->setLabelColor(Qt::white);
        plot->xAxis->grid()->setPen(QPen(QColor(80, 80, 100), 1, Qt::DotLine));
        plot->yAxis->grid()->setPen(QPen(QColor(80, 80, 100), 1, Qt::DotLine));
        plot->legend->setBrush(QBrush(QColor(30, 30, 45, 200)));
        plot->legend->setTextColor(Qt::white);
    } else { // Светлая
        plot->setBackground(QBrush(Qt::white));
        plot->xAxis->setBasePen(QPen(Qt::black));
        plot->yAxis->setBasePen(QPen(Qt::black));
        plot->xAxis->setTickPen(QPen(Qt::black));
        plot->yAxis->setTickPen(QPen(Qt::black));
        plot->xAxis->setSubTickPen(QPen(Qt::black));
        plot->yAxis->setSubTickPen(QPen(Qt::black));
        plot->xAxis->setTickLabelColor(Qt::black);
        plot->yAxis->setTickLabelColor(Qt::black);
        plot->xAxis->setLabelColor(Qt::black);
        plot->yAxis->setLabelColor(Qt::black);
        plot->xAxis->grid()->setPen(QPen(QColor(200, 200, 200), 1, Qt::DotLine));
        plot->yAxis->grid()->setPen(QPen(QColor(200, 200, 200), 1, Qt::DotLine));
        plot->legend->setBrush(QBrush(QColor(255, 255, 255, 200)));
        plot->legend->setTextColor(Qt::black);
    }

    // Автомасштаб и перерисовка
    if (!settings.keepScale) {
        plot->rescaleAxes();
    }
    
    // Перерисовка графика
    plot->replot();
}

void Visual::PlotEnergy() const {
    // Будет добавлено в следующей версии
}

void Visual::PlotMomentum() const {
    // Будет добавлено в следующей версии
}


void Visual::SetupGWPlotStyle(QCustomPlot* plot, bool darkTheme) {
    if (!plot) return;

    plot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    plot->legend->setVisible(true);
    
    QColor textColor = darkTheme ? Qt::white : Qt::black;
    QColor bgColor = darkTheme ? QColor(25, 25, 25) : Qt::white;
    QColor gridColor = darkTheme ? QColor(70, 70, 70) : QColor(220, 220, 220);

    plot->setBackground(QBrush(bgColor));
    
    // Настройка шрифтов
    QFont axisFont("Segoe UI", 9);
    
    auto setupAxis = [&](QCPAxis* axis) {
        axis->setBasePen(QPen(textColor));
        axis->setTickPen(QPen(textColor));
        axis->setSubTickPen(QPen(textColor));
        axis->setTickLabelColor(textColor);
        axis->setLabelColor(textColor);
        axis->setLabelFont(axisFont);
        axis->grid()->setPen(QPen(gridColor, 1, Qt::DotLine));
        axis->grid()->setSubGridVisible(false);
    };

    setupAxis(plot->xAxis);
    setupAxis(plot->yAxis);

    plot->legend->setBrush(QBrush(darkTheme ? QColor(40, 40, 40, 200) : QColor(255, 255, 255, 200)));
    plot->legend->setTextColor(textColor);
    plot->legend->setBorderPen(QPen(gridColor));
}


// Отрисовка орбит слияния компактных объектов
void Visual::PlotGWDynamics(QCustomPlot* plot, const QVector<double>& x1, const QVector<double>& y1, const QVector<double>& x2, const QVector<double>& y2, int tailSize, const QString& name1, const QString& name2, int fontSize) {
    if (!plot || x1.isEmpty()) return;

    // Определение диапазона точек согласно ползунку
    int totalPoints = x1.size();
    int pointsToDisplay = (totalPoints * tailSize) / 100;
    if (pointsToDisplay < 2) pointsToDisplay = 2;

    int startIdx = totalPoints - pointsToDisplay;

    // Прореживание для ускорения работы программы
    const int MAX_POINTS_ON_SCREEN = 5000; // Лимит точек для четкого рисунка
    QVector<double> x1_final, y1_final, x2_final, y2_final;

    if (pointsToDisplay > MAX_POINTS_ON_SCREEN) {
        int step = pointsToDisplay / MAX_POINTS_ON_SCREEN;
        
        for (int i = startIdx; i < totalPoints; i += step) {
            x1_final.append(x1[i]);
            y1_final.append(y1[i]);
            x2_final.append(x2[i]);
            y2_final.append(y2[i]);
        }
        // Обязательно добавляем самую последнюю точку, чтобы спираль была доведена до конца
        x1_final.append(x1.last());
        y1_final.append(y1.last());
        x2_final.append(x2.last());
        y2_final.append(y2.last());
    } else {
        // Если точек немного, просто копируем нужный кусок
        x1_final = x1.mid(startIdx);
        y1_final = y1.mid(startIdx);
        x2_final = x2.mid(startIdx);
        y2_final = y2.mid(startIdx);
    }

    // Отрисовка
    QCPCurve *orbit1 = (plot->plottableCount() > 0) ? qobject_cast<QCPCurve*>(plot->plottable(0)) : nullptr;
    QCPCurve *orbit2 = (plot->plottableCount() > 1) ? qobject_cast<QCPCurve*>(plot->plottable(1)) : nullptr;

    if (!orbit1) {
        plot->clearPlottables();
        orbit1 = new QCPCurve(plot->xAxis, plot->yAxis);
        // Задаем сплошную линию для первого объекта
        orbit1->setPen(QPen(QColor(255, 50, 50), 2, Qt::SolidLine));
        
        orbit2 = new QCPCurve(plot->xAxis, plot->yAxis);
        // Задаем штриховую линию для второго объекта
        orbit2->setPen(QPen(QColor(50, 170, 255), 2, Qt::DashLine));

        plot->legend->setVisible(true);
    }

    plot->xAxis->setLabel("X / M");
    plot->yAxis->setLabel("Y / M");
    
    QFont font = plot->xAxis->labelFont();
    font.setPointSize(fontSize);
    plot->xAxis->setLabelFont(font);
    plot->yAxis->setLabelFont(font);
    plot->xAxis->setTickLabelFont(font);
    plot->yAxis->setTickLabelFont(font);
    if (plot->legend) plot->legend->setFont(font);

    orbit1->setName(name1); 
    orbit2->setName(name2); 

    orbit1->setData(x1_final, y1_final);
    orbit2->setData(x2_final, y2_final);

    plot->replot(QCustomPlot::rpQueuedReplot);
}

// Отрисовка гравитационных волн
void Visual::PlotGWWaveform(QCustomPlot* plot, const QVector<double>& t, const QVector<double>& hp, const QVector<double>& hx, int fontSize) {
    if (!plot || t.isEmpty()) return;

    // Прореживание для ускорения работы программы 
    // const int MAX_WAVE_POINTS = 150000;
    // QVector<double> t_final, hp_final;

    // if (t.size() > MAX_WAVE_POINTS) {
    //     int step = t.size() / MAX_WAVE_POINTS;
    //     for (int i = 0; i < t.size(); i += step) {        <- Временно отлючено ->
    //         t_final.append(t[i]);
    //         hp_final.append(hp[i]);
    //     }
    //     t_final.append(t.last());
    //     hp_final.append(hp.last());
    // } else {
    //     t_final = t;
    //     hp_final = hp;
    // }

    QVector<double> t_final = t;
    QVector<double> hp_final = hp;

    QCPGraph *graph = (plot->graphCount() > 0) ? plot->graph(0) : plot->addGraph();
    graph->setData(t_final, hp_final);
    graph->setPen(QPen(Qt::red, 1));

    // Подписи осей
    plot->xAxis->setLabel("Время (t/M)");
    plot->yAxis->setLabel("Деформация r·h/M");

    // Регулировка размера шрифта
    QFont font = plot->xAxis->labelFont();
    font.setPointSize(fontSize);
    plot->xAxis->setLabelFont(font);
    plot->yAxis->setLabelFont(font);
    plot->xAxis->setTickLabelFont(font);
    plot->yAxis->setTickLabelFont(font);
    plot->legend->setFont(font);

    plot->xAxis->setLabel("Время, t/M");
    plot->yAxis->setLabel("Амплитуда, rh₊/M");

    plot->yAxis->setRange(-0.003, 0.003);
    plot->xAxis->setRange(-1500, 250);

    QSharedPointer<QCPAxisTickerFixed> fixedTicker(new QCPAxisTickerFixed);
    fixedTicker->setTickStep(500.0);
    plot->xAxis->setTicker(fixedTicker);

    //plot->rescaleAxes();
    plot->replot(QCustomPlot::rpQueuedReplot);
}