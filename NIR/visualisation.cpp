#include "visualisation.h"

Visual::Visual(std::vector<Body>& bodies_, std::vector<double>& time_values_, std::vector<double>& energy_values_, std::vector<double>& momentum_values_) :
    bodies(bodies_), time_values(time_values_), energy_values(energy_values_), momentum_values(momentum_values_) {}

std::string formatTime(double time) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(1) << std::round(time * 365 * 10) / 10;
    return oss.str();
}

//Визуализация орбит
void Visual::PlotOrbits(const std::string& filename, const std::vector<std::string>& labels, const std::vector<std::string>& colors) const {
    std::ifstream inp(filename);
    if (!inp.is_open()) {
        std::cerr << "Ошибка: не удалось открыть файл " << filename << std::endl;
        return;
    }

    std::vector<double> time;
    std::vector<std::vector<double>> x, y, z;


    for (int i = 0; i < bodies.size(); ++i) {
        x.push_back(std::vector<double>());
        y.push_back(std::vector<double>());
        z.push_back(std::vector<double>());
    }

    // Код для сохранения изображений планет в каждый момент времени 

    //int frame_count = 0;
    //const int skip_frames = 300;
    //plt::figure_size(1300, 1300);
    //double t;
    //while (inp >> t) {
    //    time.push_back(t);
    //    for (int i = 0; i < bodies.size(); ++i) {
    //        double xi, yi, zi;
    //        inp >> xi >> yi >> zi;
    //        x[i].push_back(xi / AU);
    //        y[i].push_back(yi / AU);
    //        z[i].push_back(zi / AU);
    //    }

    //    if (frame_count++ % skip_frames != 0) continue;

    //    plt::clf();
    //    plt::rcparams({ { "font.size", "30" } });

    //    for (int i = 0; i < bodies.size(); ++i) {
    //        if (x[i].size() < 2) continue;

    //        // Линия без последнего отрезка
    //        std::vector<double> x_part(x[i].begin(), x[i].end() - 1);
    //        std::vector<double> y_part(y[i].begin(), y[i].end() - 1);

    //        plt::plot(x_part, y_part, {
    //            {"color", colors[i]},
    //            {"linewidth", "3"},
    //            {"label", labels[i]}
    //            });
    //    }

    //    // 2. Рисуем ПОСЛЕДНИЙ ОТРЕЗОК отдельно (тоньше)
    //    for (int i = 0; i < bodies.size(); ++i) {
    //        if (x[i].size() < 2) continue;

    //        std::vector<double> last_segment_x = { x[i][x[i].size() - 2], x[i].back() };
    //        std::vector<double> last_segment_y = { y[i][y[i].size() - 2], y[i].back() };

    //        plt::plot(last_segment_x, last_segment_y, {
    //            {"color", colors[i]},
    //            {"linewidth", "2.5"}  // Тонкая линия
    //            });
    //    }

    //    // 3. Рисуем точки (кружки) поверх
    //    for (int i = 1; i < bodies.size(); ++i) {
    //        if (x[i].empty()) continue;

    //        // Большой белый круг (подложка)
    //        plt::plot({ x[i].back() }, { y[i].back() }, {
    //            {"color", "white"},
    //            {"marker", "o"},
    //            {"markersize", "0.5"},
    //            {"markeredgewidth", "0"},
    //            {"linestyle", "none"}
    //            });

    //        // Цветной круг с обводкой
    //        plt::plot({ x[i].back() }, { y[i].back() }, {
    //            {"color", colors[i]},
    //            {"marker", "o"},
    //            {"markersize", "25"},
    //            {"markeredgecolor", "black"},
    //            {"markeredgewidth", "1.5"},
    //            {"linestyle", "none"}
    //            });
    //    }

    //    //for (int i = 5; i < bodies.size(); ++i) {
    //    //    if (x[i].empty()) continue;

    //    //    // Большой белый круг (подложка)
    //    //    plt::plot({ x[i].back() }, { y[i].back() }, {
    //    //        {"color", "white"},
    //    //        {"marker", "o"},
    //    //        {"markersize", "0.5"},
    //    //        {"markeredgewidth", "0"},
    //    //        {"linestyle", "none"}
    //    //        });

    //    //    // Цветной круг с обводкой
    //    //    plt::plot({ x[i].back() }, { y[i].back() }, {
    //    //        {"color", colors[i]},
    //    //        {"marker", "o"},
    //    //        {"markersize", "25"},
    //    //        {"markeredgecolor", "black"},
    //    //        {"markeredgewidth", "1.5"},
    //    //        {"linestyle", "none"}
    //    //        });
    //    //}

    //    std::string str = "planet" + std::to_string(t) + ".png";
    //    std::string title = "t = " + formatTime(t) + " (годы)";
    //    plt::title(title, {{"fontsize", "30"}});
    //    plt::xlim(-0.1, 0.1);
    //    plt::ylim(-0.1, 0.1);
    //    plt::xlabel("X (а.е.)");
    //    plt::ylabel("Y (а.е.)");

    //    plt::legend({
    //    {"loc", "upper right"},
    //    {"fontsize", "15"}
    //        });
    //    plt::grid(true);
    //    //plt::axis("equal");
    //    plt::backend("Agg");
    //    plt::save("D:\\image_toi_270\\"+str);
    //}
    //plt::legend();

    //int frame_count = 0;
    //const int skip_frames = 300;
    //plt::figure_size(1300, 1300);
    //double t;
    //while (inp >> t) {
    //    time.push_back(t);
    //    for (int i = 0; i < bodies.size(); ++i) {
    //        double xi, yi, zi;
    //        inp >> xi >> yi >> zi;
    //        x[i].push_back(xi / AU);
    //        y[i].push_back(yi / AU);
    //        z[i].push_back(zi / AU);
    //    }

    //    if (frame_count++ % skip_frames != 0) continue;

    //    plt::clf();
    //    plt::rcparams({ { "font.size", "20" } });

    //    for (int i = 0; i < bodies.size(); ++i) {
    //        if (x[i].size() < 2) continue;

    //        // Линия без последнего отрезка
    //        std::vector<double> x_part(x[i].begin(), x[i].end() - 1);
    //        std::vector<double> y_part(z[i].begin(), z[i].end() - 1);

    //        plt::plot(x_part, y_part, {
    //            {"color", colors[i]},
    //            {"linewidth", "3"},
    //            {"label", labels[i]}
    //            });
    //    }

    //    // 2. Рисуем ПОСЛЕДНИЙ ОТРЕЗОК отдельно (тоньше)
    //    for (int i = 0; i < bodies.size(); ++i) {
    //        if (x[i].size() < 2) continue;

    //        std::vector<double> last_segment_x = { x[i][x[i].size() - 2], x[i].back() };
    //        std::vector<double> last_segment_y = { z[i][z[i].size() - 2], z[i].back() };

    //        plt::plot(last_segment_x, last_segment_y, {
    //            {"color", colors[i]},
    //            {"linewidth", "2.5"}  // Тонкая линия
    //            });
    //    }

    //    // 3. Рисуем точки (кружки) поверх
    //    for (int i = 1; i < bodies.size(); ++i) {
    //        if (x[i].empty()) continue;

    //        // Большой белый круг (подложка)
    //        plt::plot({ x[i].back() }, { z[i].back() }, {
    //            {"color", "white"},
    //            {"marker", "o"},
    //            {"markersize", "0.5"},
    //            {"markeredgewidth", "0"},
    //            {"linestyle", "none"}
    //            });

    //        // Цветной круг с обводкой
    //        plt::plot({ x[i].back() }, { z[i].back() }, {
    //            {"color", colors[i]},
    //            {"marker", "o"},
    //            {"markersize", "25"},
    //            {"markeredgecolor", "black"},
    //            {"markeredgewidth", "1.5"},
    //            {"linestyle", "none"}
    //            });
    //    }

    //    std::string str = "planet" + std::to_string(t) + ".png";
    //    std::string title = "t = " + formatTime(t) + " (дни)";
    //    plt::title(title, { {"fontsize", "20"} });
    //    plt::xlim(-0.1, 0.1);
    //    plt::ylim(-0.1, 0.1);
    //    plt::xlabel("X (а.е.)");
    //    plt::ylabel("Z (а.е.)");

    //    plt::legend({
    //    {"loc", "upper right"},
    //    {"fontsize", "15"}
    //        });
    //    plt::grid(true);
    //    plt::axis("equal");
    //    plt::backend("Agg");
    //    plt::save("D:\\image_toi_270\\" + str);
    //}
    //plt::legend();

    double t;
    while (inp >> t) {
        time.push_back(t);
        for (int i = 0; i < bodies.size(); ++i) {
            double xi, yi, zi;
            inp >> xi >> yi >> zi;
            x[i].push_back(xi / AU);
            y[i].push_back(yi / AU);     // Для визуализации в астрономических единицах
            z[i].push_back(zi / AU);
            /*x[i].push_back(xi / 1000);
            y[i].push_back(yi / 1000);   // Для визуализации в километрах
            z[i].push_back(zi / 1000);*/
        }
    }

    inp.close();

    //График орбит на осях x и y
    plt::figure_size(1000, 1000);
    plt::rcparams({ { "font.size", "21" } });

    for (int i = 0; i < bodies.size(); ++i) {
        if (i == 0) plt::scatter(x[i], y[i], 100.0,{ {"label", labels[i]}, {"color", colors[i]} });
        else plt::plot(x[i], y[i], { {"label", labels[i]}, {"color", colors[i]} });
        //else plt::scatter(x[i], y[i], 50.0, { {"label", labels[i]}, {"color", colors[i]} });
    }
    //plt::title("Орбиты тел (X и Y)");
    plt::xlim(-30, 30);
    plt::xlabel("X (a.e.)");
    plt::ylabel("Y (a.e.)");
    plt::legend({ { "fontsize", "15" } });
    //plt::legend();
    plt::grid(true);
    plt::axis("equal");

    // График орбит на осях x и z
    plt::figure_size(1000, 1000);
    plt::rcparams({ { "font.size", "18" } });

    for (int i = 0; i < bodies.size(); ++i) {
        plt::plot(x[i], z[i], { {"label", labels[i]}, {"color", colors[i]} });
    }
    plt::title("Орбиты тел (X и Z)");
    plt::xlim(-30 / AU, 30 / AU);
    plt::xlabel("X (а.е.)");
    plt::ylabel("Z (а.е.)");
    plt::legend({ { "fontsize", "15" } });
    plt::grid(true);
    plt::axis("equal");
    plt::show();
}



void Visual::PlotEnergy() const {
    if (time_values.empty() || energy_values.empty()) {
        std::cerr << "Ошибка: данные для построения графика энергии отсутствуют.\n";
        return;
    }
    plt::figure_size(1000, 1000);
    plt::rcparams({ { "font.size", "25" } });

    plt::plot(time_values, energy_values, "-");
    plt::title("Закон сохранения энергии");
    plt::xlabel("Время (годы)");
    plt::ylabel("Полная энергия (Дж)");
    plt::grid(true);
    plt::tight_layout();
    plt::show();
}

void Visual::PlotMomentum() const {
    if (time_values.empty() || momentum_values.empty()) {
        std::cerr << "Ошибка: данные для построения графика импульса отсутствуют.\n";
        return;
    }

    plt::figure_size(900, 900);
    plt::rcparams({ { "font.size", "25" } });

    plt::plot(time_values, momentum_values, "-");
    plt::title("Закон сохранения импульса");
    plt::xlabel("Время (годы)");
    plt::ylabel("Импульс");
    plt::grid(true);
    plt::tight_layout();
    plt::show();
}