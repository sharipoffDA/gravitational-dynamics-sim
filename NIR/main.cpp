#include"system.h"
#include <numbers>
//#include <QMainWindow>
//#include <QLabel>
//#include <QApplication>
//#include <QMainWindow>
//#include <QPushButton>
//#include <QVBoxLayout>
//#include <QWidget>
//#include <QDebug>
#include <Python.h>
#include <cstdlib>
const double pi = std::numbers::pi;

struct KeplerianElements {
    std::string name;
    double a_au;        // Большая полуось (а.е.)
    double incl_deg;    // Наклон орбиты (градусы)
    double omega_deg;   // Аргумент периастра (градусы)
    double Omega_deg;   // Долгота восходящего узла (градусы)
    double mass;  // Масса (в массах Земли)
    double ecc;         // Эксцентриситет
    double period_days; // Период обращения (дни)
    double T_p_year;    // Время прохождения периастра (год)
};

// Для экзопланетных систем
//void compute_initial_conditions(
//    const NasaPlanetData& data,
//    double M_body,
//    double& x, double& y, double& z,
//    double& vx, double& vy, double& vz
//) {
//
//    double a = data.a_au * AU;
//    double i = data.incl_deg * pi / 180.0;
//    double omega = data.omega_deg * pi / 180.0;
//    double Omega = data.Omega_deg * pi / 180.0;
//    double e = data.ecc;
//    double mass = data.mass_earth * M_EARTH;
//
//    // Положение в перицентре (θ = 0)
//    double r_peri = a * (1 - e); // Расстояние от фокуса эллипса до перицентра
//    x = r_peri * (cos(Omega) * cos(omega) - sin(Omega) * sin(omega) * cos(i));
//    y = r_peri * (sin(Omega) * cos(omega) + cos(Omega) * sin(omega) * cos(i));
//    z = r_peri * sin(omega) * sin(i);
//
//    // Скорость в перицентре
//    double mu = G * (M_body + mass); // Гравитационный параметр
//    double v_peri = sqrt((mu / a) * (1 + e) / (1 - e)); // Скорость в перицентре
//    vx = -v_peri * (cos(Omega) * sin(omega) + sin(Omega) * cos(omega) * cos(i));
//    vy = -v_peri * (sin(Omega) * sin(omega) - cos(Omega) * cos(omega) * cos(i));
//    vz = v_peri * cos(omega) * sin(i);
//}

// Для звёзд и компактных сверхмассивных объектов
void compute_initial_conditions(
    const KeplerianElements& data,
    double M_star,
    double epoch_year,  // Начальная эпоха моделирования (например, 2002.0)
    double& x, double& y, double& z,
    double& vx, double& vy, double& vz
) {
    // Перевод в радианы
    double a = data.a_au * AU;
    //double a = data.a_au;
    double i = data.incl_deg * DEG_TO_RAD;
    double omega = data.omega_deg * DEG_TO_RAD;
    double Omega = data.Omega_deg * DEG_TO_RAD;
    double e = data.ecc;
    //double mass = data.mass_earth * M_EARTH; 
    double mass = data.mass * M_SUN; // исправление на массы солнца
    double T_p = data.T_p_year;  

    // 1. Вычислить среднюю аномалию на epoch_year
    double period_years = data.period_days / 365.25;
    double n = 2.0 * M_PI / period_years;  // среднее движение (рад/год)
    double M = n * (epoch_year - T_p);     // средняя аномалия
    M = fmod(M + 2 * M_PI, 2 * M_PI);

    // 2. Решить уравнение Кеплера для эксцентрической аномалии (E)
    double E = M;
    double tol = 1e-8;
    for (int iter = 0; iter < 100; ++iter) {
        double dE = (E - e * sin(E) - M) / (1 - e * cos(E));
        E -= dE;
        if (fabs(dE) < tol) break;
    }

    // 3. Найти истинную аномалию (ν)
    double nu = 2.0 * atan2(sqrt(1 + e) * sin(E / 2), sqrt(1 - e) * cos(E / 2));

    // 4. Вычислить радиус-вектор
    double r = a * (1 - e * cos(E));

    // 5. Координаты в орбитальной плоскости
    double x_orb = r * cos(nu);
    double y_orb = r * sin(nu);

    // 6. Преобразование в 3D пространство
    x = (cos(Omega) * cos(omega) - sin(Omega) * sin(omega) * cos(i)) * x_orb +
        (-cos(Omega) * sin(omega) - sin(Omega) * cos(omega) * cos(i)) * y_orb;
    y = (sin(Omega) * cos(omega) + cos(Omega) * sin(omega) * cos(i)) * x_orb +
        (-sin(Omega) * sin(omega) + cos(Omega) * cos(omega) * cos(i)) * y_orb;
    z = (sin(omega) * sin(i)) * x_orb + (cos(omega) * sin(i)) * y_orb;

    // 7. Скорости (производные)
    double mu = G * (M_star + mass);
    double p = a * (1 - e * e);  // фокальный параметр
    double h = sqrt(mu * p);     // удельный момент импульса

    double sqrt_mu_p = sqrt(mu / p);
    double v_x_orb = -sqrt_mu_p * sin(nu);
    double v_y_orb = sqrt_mu_p * (e + cos(nu));

    vx = (cos(Omega) * cos(omega) - sin(Omega) * sin(omega) * cos(i)) * v_x_orb +
        (-cos(Omega) * sin(omega) - sin(Omega) * cos(omega) * cos(i)) * v_y_orb;
    vy = (sin(Omega) * cos(omega) + cos(Omega) * sin(omega) * cos(i)) * v_x_orb +
        (-sin(Omega) * sin(omega) + cos(Omega) * cos(omega) * cos(i)) * v_y_orb;
    vz = sin(omega) * sin(i) * v_x_orb + cos(omega) * sin(i) * v_y_orb;
}

//class MainWindow : public QMainWindow {
//    Q_OBJECT
//public:
//    MainWindow(QWidget* parent = nullptr);
//    ~MainWindow();
//
//private slots:
//    void visualize();
//
//private:
//    QWidget* plotWidget;
//};
//
//class MainWindow : public QMainWindow {
//    Q_OBJECT
//public:
//    MainWindow(QWidget* parent = nullptr) : QMainWindow(parent) {
//        // Настройка основного окна
//        QWidget* centralWidget = new QWidget(this);
//        setCentralWidget(centralWidget);
//        QVBoxLayout* layout = new QVBoxLayout(centralWidget);
//
//        // Кнопка
//        QPushButton* button = new QPushButton("Визуализировать", this);
//        layout->addWidget(button);
//
//        // Пустой QWidget для графика
//        plotWidget = new QWidget(this);
//        layout->addWidget(plotWidget);
//
//        // Подключение кнопки
//        connect(button, &QPushButton::clicked, this, &MainWindow::visualize);
//
//        // Инициализация Python
//        Py_Initialize();
//        if (!Py_IsInitialized()) {
//            qDebug() << "Ошибка: Python не инициализирован";
//        }
//
//        // Установка QtAgg backend
//        PyRun_SimpleString("import matplotlib; matplotlib.use('Qt5Agg')");
//    }
//
//    ~MainWindow() {
//        Py_Finalize();
//    }
//
//private slots:
//    void visualize() {
//        // Данные для графика (пример)
//        std::vector<double> x = { 1, 2, 3, 4 };
//        std::vector<double> y = { 1, 4, 2, 3 };
//
//        // Создание графика с matplotlib-cpp
//        plt::figure(); // Новая фигура
//        plt::plot(x, y, "-");
//        plt::title("Пример графика");
//        plt::xlabel("X");
//        plt::ylabel("Y");
//        //plt::legend();
//
//        // Показать график (интеграция с QtAgg)
//        plt::show();
//    }
//
//private:
//    QWidget* plotWidget;
//};

int main(int argc, char* argv[]) {

    // Инициализация QApplication перед Python
    /*QApplication app(argc, argv);

    _putenv("PYTHONHOME=D:/Python");
    Py_Initialize();

    MainWindow window;
    window.resize(800, 600);
    window.show();
    return app.exec();*/
    /*QApplication app(argc, argv);

    QMainWindow window;
    QLabel label("Я использую Qt бесплатно!");
    window.setCentralWidget(&label);
    window.show();*/

    //return app.exec();
    std::vector<Body> solar_system_bodies = {
            {1.989e30, -8.572865039469250E+05 * 1000, -7.346088835335051E+05 * 1000, 2.685423265526889E+04 * 1000, 1.239859639798033E-02 * 1000, -6.348466611140617E-03 * 1000, -2.037876555553517E-04 * 1000}, // Солнце         
            {3.302e23, -5.879699189509091E+07 * 1000, -2.492820404148239E+07 * 1000, 3.364042452841429E+06 * 1000, 8.711982611106873E+00 * 1000, -4.284856986770977E+01 * 1000, -4.299279282370732E+00 * 1000}, // Меркурий        
            {4.867e24, 6.697319534635594E+07 * 1000, 8.337171945245868E+07 * 1000, -2.731933993919346E+06 * 1000, -2.735548307021769E+01 * 1000, 2.182743070706988E+01 * 1000, 1.878804135388283E+00 * 1000}, // Венера        
            {5.972e24, -2.758794880287251E+07 * 1000, 1.439239583084676E+08 * 1000, 1.921064327326417E+04 * 1000, -2.977686364628585E+01 * 1000, -5.535813340802556E+00 * 1000, -1.943387942073826E-04 * 1000}, // Земля        
            {6.417e23, -7.890038131682467E+07 * 1000, 2.274372361241295E+08 * 1000, 6.722196400986686E+06 * 1000, -2.199759485544059E+01 * 1000, -5.787405095467102E+00 * 1000, 4.184257990348734E-01 * 1000}, // Марс       
            {1.898e27, 1.571230833020991E+08 * 1000, 7.429840488421507E+08 * 1000, -6.597049828231782E+06 * 1000, -1.293244436609816E+01 * 1000, 3.325781476287804E+00 * 1000, 2.755437569190042E-01 * 1000}, // Юпитер 
            {5.683e26, 1.414498231862034E+09 * 1000, -2.647172137275474E+08 * 1000, -5.171551879510410E+07 * 1000, 1.240660798615463E+00 * 1000, 9.473546595187154E+00 * 1000, -2.135791731559418E-01 * 1000}, // Сатурн
            {8.681e25, 1.660221941282016E+09 * 1000, 2.406965914394302E+09 * 1000, -1.256903703858864E+07 * 1000, -5.655920544550284E+00 * 1000, 3.549247198028906E+00 * 1000, 8.651776992957205E-02 * 1000}, // Уран
            {1.024e26, 4.469116222588663E+09 * 1000, -9.560778256566879E+07 * 1000, -1.010264767638457E+08 * 1000, 8.064561683471368E-02 * 1000, 5.465730017544922E+00 * 1000, -1.151205185674022E-01 * 1000}, // Нептун
    };


    BodySystem solar_system(solar_system_bodies);


    double t_end = 200 * year_to_sec; 
    double dt = 6*3600;  // Шаг времени

    
    solar_system.SaveToFile_RK6("solar_system_orbits.txt", t_end, dt, false);

    
    std::vector<std::string> labels = { "Солнце", "Меркурий", "Венера", "Земля", "Марс", "Юпитер", "Сатурн", "Уран", "Нептун" };
    
    std::vector<std::string> colors = { "yellow", "gray", "orange", "blue", "red", "brown", "purple", "black", "navy" };
    
    
    solar_system.visual.PlotOrbits("solar_system_orbits.txt", labels, colors);
    solar_system.visual.PlotEnergy();
    solar_system.visual.PlotMomentum();


    /*const std::vector<NasaPlanetData> TOI_270 = {
        {"b", 0.03197, 89.39, 0.0, 0.0, 1.58, 0.034, 3.3601538},
        {"c", 0.04526, 89.36, 10.0, 0.0, 6.15, 0.027, 5.6605731},
        {"d", 0.07210, 89.73, -6.0, 0.0, 4.78, 0.032, 11.379573},
    };

    std::vector<Body> toi270_system_bodies;

    toi270_system_bodies.push_back({ M_STAR_TOI270, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 });

    double x, y, z, vx, vy, vz;

    for (const auto& planet : TOI_270) {
        compute_initial_conditions(planet, M_STAR_TOI270, x, y, z, vx, vy, vz);

        toi270_system_bodies.push_back({ planet.mass_earth * M_EARTH, x, y, z, vx, vy, vz });
    }

    PlanetarySystem toi270_system(toi270_system_bodies);

    double t_end = 35 * 24 * 3600;
    double dt = 60.0;

    toi270_system.SaveToFile_RK6("toi270_orbits.txt", t_end, dt, false);

    std::vector<std::string> labels = { "TOI-270", "b", "c", "d" };
    std::vector<std::string> colors = { "red", "orange", "yellow", "green" };

    toi270_system.visual.PlotOrbits("toi270_orbits.txt", labels, colors);*/

    //struct KeplerianElements {
    //    std::string name;
    //    double a_au;        // Большая полуось (а.е.)
    //    double incl_deg;    // Наклон орбиты (градусы)
    //    double omega_deg;   // Аргумент периастра (градусы)
    //    double Omega_deg;   // Долгота восходящего узла (градусы)
    //    double mass;  // Масса (в массах Земли)
    //    double ecc;         // Эксцентриситет
    //    double period_days; // Период обращения (дни)
    //    double T_p_year;    // Время прохождения периастра (год)
    //};

    // Эксперимент: звёзды вокруг чёрной дыры 
    /*const std::vector<KeplerianElements> STAR_SYSTEM = {
        {"S2", 1030.4, 136.78, 71.36, 234.50, 15.0, 0.884, 5897.0, 2002.32},
        {"S38", 1144.0, 166.22, 18.4, 101.8, 10.0, 0.818, 6790.0, 2003.30},
        {"S55/S0-102", 890.0, 141.7, 133.5, 129.9, 8.0, 0.74, 4380.0, 2009.31},
    };
    
    std::vector<Body> star_system_bodies;

    star_system_bodies.push_back({ 4.15e6 * M_SUN, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 });

    double x, y, z, vx, vy, vz;

    for (const auto& body : STAR_SYSTEM) {
        compute_initial_conditions(body, 4.15e6 * M_SUN, 2002.0, x, y, z, vx, vy, vz);

        star_system_bodies.push_back({ body.mass * M_SUN, x, y, z, vx, vy, vz });
    }

    BodySystem star_system(star_system_bodies);

    double t_end = 50 * year_to_sec;
    double dt = 6 * 3600.0;

    star_system.SaveToFile_RK4("star_system_orbits.txt", t_end, dt, true);

    std::vector<std::string> labels = { "Sgr A*", "S2", "S38", "S55/S0-102" };
    std::vector<std::string> colors = { "red", "orange", "yellow", "green" };

    star_system.visual.PlotOrbits("star_system_orbits.txt", labels, colors);*/


    //std::string name;
    //double a_au;        // Большая полуось (а.е.)
    //double incl_deg;    // Наклон орбиты (градусы)
    //double omega_deg;   // Аргумент периастра (градусы)
    //double Omega_deg;   // Долгота восходящего узла (градусы)
    //double mass_earth;  // Масса (в массах Земли)
    //double ecc;         // Эксцентриситет
    //double period_days; // Период обращения (дни)
    //double T_p_year;    // Время прохождения периастра (год)

    //const std::vector<NasaPlanetData> NSNS_SYSTEM = {
    //{"NS1", 5.4e4, 0.0, 0.0, 0.0, 1.499, 0.0, 2.43e-7, 0.0},  // NS1, a в единицах M_SUN
    //{"NS2", 5.4e4, 0.0, 0.0, 0.0, 1.499, 0.0, 2.43e-7, 0.0}   // NS2, идентичные параметры
    //};

    //std::vector<Body> nsns_system_bodies;

    //double M_total = 1.499 * M_SUN + 1.499 * M_SUN;  // 2.998 M_SUN
    //nsns_system_bodies.push_back({ 1.499 * M_SUN, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 });  // NS1 (временная заглушка)

    //double x1, y1, z1, vx1, vy1, vz1;
    //double x2, y2, z2, vx2, vy2, vz2;

    //// Вычисляем начальные условия для обеих звёзд
    //for (size_t i = 0; i < NSNS_SYSTEM.size(); ++i) {
    //    const auto& body = NSNS_SYSTEM[i];
    //    double x, y, z, vx, vy, vz;
    //    compute_initial_conditions(body, M_total, 0.0, x, y, z, vx, vy, vz);  // Используем M_total как референс

    //    if (i == 0) {
    //        // NS1: смещаем влево относительно центра масс
    //        x1 = -0.5 * x;
    //        y1 = y;
    //        z1 = z;
    //        vx1 = -0.5 * vx;
    //        vy1 = vy;
    //        vz1 = vz;
    //        nsns_system_bodies[0].x = x1;
    //        nsns_system_bodies[0].vx = vx1;
    //        nsns_system_bodies[0].vy = vy1;
    //        nsns_system_bodies[0].vz = vz1;
    //    }
    //    else {
    //        // NS2: смещаем вправо относительно центра масс
    //        x2 = 0.5 * x;
    //        y2 = y;
    //        z2 = z;
    //        vx2 = 0.5 * vx;
    //        vy2 = vy;
    //        vz2 = vz;
    //        nsns_system_bodies.push_back({ 1.499 * M_SUN, x2, y2, z2, vx2, vy2, vz2 });  // NS2
    //    }
    //}

    //PlanetarySystem nsns_system(nsns_system_bodies);

    //// Параметры времени
    //double t_end = 0.05;  // ~2.38 орбиты для тестирования
    //double dt = 1e-6;     // Уменьшен до ~P/21000 для стабильности с 2.5PN

    //nsns_system.SaveToFile_RK4("nsns_system_orbits.txt", t_end, dt, true);  // true для 1PN/2.5PN

    //std::vector<std::string> labels = { "NS1", "NS2" };
    //std::vector<std::string> colors = { "red", "orange" };

    //nsns_system.visual.PlotOrbits("nsns_system_orbits.txt", labels, colors);

    //const std::vector<KeplerianElements> BHNS_SYSTEM = {
    //{"NS", 5.57e5 / AU, 0.0, 0.0, 0.0, 2.20, 0.0, 3.82e-7, 0.0}  // a в метрах
    //};
    //std::vector<Body> bhns_system_bodies;
    //double M_total = 7.40 * M_SUN + 2.20 * M_SUN;
    //bhns_system_bodies.push_back({ 7.40 * M_SUN, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 });
    //double x, y, z, vx, vy, vz;
    //for (const auto& body : BHNS_SYSTEM) {
    //    compute_initial_conditions(body, 7.40 * M_SUN, 0.0, x, y, z, vx, vy, vz);
    //    // Сдвиг в центр масс
    //    double x_cm = (7.40 * M_SUN * 0.0 + 2.20 * M_SUN * x) / M_total;
    //    double x_bh = -(2.20 * M_SUN / M_total) * x;
    //    double x_ns = (7.40 * M_SUN / M_total) * x;
    //    double vx_bh = -(2.20 * M_SUN / M_total) * vx;
    //    double vy_bh = -(2.20 * M_SUN / M_total) * vy;
    //    double vz_bh = -(2.20 * M_SUN / M_total) * vz;
    //    double vx_ns = (7.40 * M_SUN / M_total) * vx;
    //    double vy_ns = (7.40 * M_SUN / M_total) * vy;
    //    double vz_ns = (7.40 * M_SUN / M_total) * vz;
    //    bhns_system_bodies[0].x = x_bh;    // Обновляем позицию BH
    //    bhns_system_bodies[0].vx = vx_bh;  // Обновляем скорость BH
    //    bhns_system_bodies[0].vy = vy_bh;
    //    bhns_system_bodies[0].vz = vz_bh;
    //    bhns_system_bodies.push_back({ 2.20 * M_SUN, x_ns, y, z, vx_ns, vy_ns, vz_ns }); // Используем скорректированные значения для NS
    //    //bhns_system_bodies.push_back({ 2.20 * M_SUN, x, y, z, vx, vy, vz }); // Используем скорректированные значения для NS
    //}
    //BodySystem bhns_system(bhns_system_bodies);
    //// Исправленные параметры времени
    //double t_end = 0.25;  // 0.05 с для 5 орбит
    //double dt = 0.02592/1e3;   // ~P/1000
    //bhns_system.SaveToFile_RK4("bhns_system_orbits.txt", t_end, dt, true);
    //std::vector<std::string> labels = { "ЧД", "НЗ" };
    //std::vector<std::string> colors = { "red", "orange" };
    //bhns_system.visual.PlotOrbits("bhns_system_orbits.txt", labels, colors);


    //const std::vector<KeplerianElements> NSNS_SYSTEM = {
    //{"NS", 51600.66 / AU, 0.0, 0.0, 0.0, 1.388, 0.0, 4.030588e-08, 0.0}
    //};
    //std::vector<Body> nsns_system_bodies;
    //double M_total = 1.982 * M_SUN + 1.388 * M_SUN;
    //nsns_system_bodies.push_back({ 1.982 * M_SUN, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 });
    //double x, y, z, vx, vy, vz;
    //for (const auto& body : NSNS_SYSTEM) {
    //    compute_initial_conditions(body, 1.982 * M_SUN, 0.0, x, y, z, vx, vy, vz);
    //    // Сдвиг в центр масс
    //    double x_cm = (1.982 * M_SUN * 0.0 + 1.388 * M_SUN * x) / M_total;
    //    double x_bh = -(1.388 * M_SUN / M_total) * x;
    //    double x_ns = (1.982 * M_SUN / M_total) * x;
    //    double vx_bh = -(1.388 * M_SUN / M_total) * vx;
    //    double vy_bh = -(1.388 * M_SUN / M_total) * vy;
    //    double vz_bh = -(1.388 * M_SUN / M_total) * vz;
    //    double vx_ns = (1.982 * M_SUN / M_total) * vx;
    //    double vy_ns = (1.982 * M_SUN / M_total) * vy;
    //    double vz_ns = (1.982 * M_SUN / M_total) * vz;
    //    nsns_system_bodies[0].x = x_bh;    // Обновляем позицию BH
    //    nsns_system_bodies[0].vx = vx_bh;  // Обновляем скорость BH
    //    nsns_system_bodies[0].vy = vy_bh;
    //    nsns_system_bodies[0].vz = vz_bh;
    //    nsns_system_bodies.push_back({ 1.388 * M_SUN, x_ns, y, z, vx_ns, vy_ns, vz_ns }); // Используем скорректированные значения для NS
    //}
    //BodySystem bhns_system(nsns_system_bodies);

    //double t_end = 0.05; 
    //double dt = 4.030588e-08;

    //bhns_system.SaveToFile_RK6("nsns_system_orbits.txt", t_end, dt, true);
    //std::vector<std::string> labels = { "НЗ", "НЗ" };
    //std::vector<std::string> colors = { "red", "orange" };
    //bhns_system.visual.PlotOrbits("nsns_system_orbits.txt", labels, colors);

    //const std::vector<KeplerianElements> BHBH_SYSTEM = {
    //    {"BH2", 7.58e-6, 0.0, 0.0, 0.0, 29.0, 8e-4, 9.5e-7, 0.0}
    //};
    //std::vector<Body> bhbh_system_bodies;
    //double M_total = 36.0 * M_SUN + 29.0 * M_SUN;
    //bhbh_system_bodies.push_back({ 36.0 * M_SUN, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 });
    //double x, y, z, vx, vy, vz;
    //for (const auto& body : BHBH_SYSTEM) {
    //    compute_initial_conditions(body, 36.0 * M_SUN, 0.0, x, y, z, vx, vy, vz);
    //    // Сдвиг в центр масс
    //    double x_cm = (36.0 * M_SUN * 0.0 + 1.388 * M_SUN * x) / M_total;
    //    double x_bh = -(29.0 * M_SUN / M_total) * x;
    //    double x_ns = (36.0 * M_SUN / M_total) * x;
    //    double vx_bh = -(29.0 * M_SUN / M_total) * vx;
    //    double vy_bh = -(29.0 * M_SUN / M_total) * vy;
    //    double vz_bh = -(29.0 * M_SUN / M_total) * vz;
    //    double vx_ns = (36.0 * M_SUN / M_total) * vx;
    //    double vy_ns = (36.0 * M_SUN / M_total) * vy;
    //    double vz_ns = (36.0 * M_SUN / M_total) * vz;
    //    bhbh_system_bodies[0].x = x_bh;    // Обновляем позицию BH
    //    bhbh_system_bodies[0].vx = vx_bh;  // Обновляем скорость BH
    //    bhbh_system_bodies[0].vy = vy_bh;
    //    bhbh_system_bodies[0].vz = vz_bh;
    //    bhbh_system_bodies.push_back({ 29.0 * M_SUN, x_ns, y, z, vx_ns, vy_ns, vz_ns }); // Используем скорректированные значения для NS
    //}
    //BodySystem bhns_system(bhbh_system_bodies);

    //double t_end = 0.5;
    //double dt = 0.07776e-3;

    //bhns_system.SaveToFile_RK6("bhbh_system_orbits.txt", t_end, dt, true);
    //std::vector<std::string> labels = { "ЧД", "ЧД" };
    //std::vector<std::string> colors = { "red", "orange" };
    //bhns_system.visual.PlotOrbits("bhbh_system_orbits.txt", labels, colors);
    return 0;
}

//#include "main.moc"