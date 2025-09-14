#include"system.h"
#include <numbers>
const double pi = std::numbers::pi;

struct NasaPlanetData {
    std::string name;
    double a_au;        // Большая полуось (а.е.)
    double incl_deg;    // Наклон орбиты (градусы)
    double omega_deg;   // Аргумент периастра (градусы)
    double Omega_deg;   // Долгота восходящего узла (градусы)
    double mass_earth;  // Масса (в массах Земли)
    double ecc;         // Эксцентриситет
    double period_days; // Период обращения (дни)
};


void compute_initial_conditions(
    const NasaPlanetData& data,
    double M_star,
    double& x, double& y, double& z,
    double& vx, double& vy, double& vz
) {
    // Перевод в единицы изммерения, используемые в программе
    double a = data.a_au * AU;
    double i = data.incl_deg * pi / 180.0;
    double omega = data.omega_deg * pi / 180.0;
    double Omega = data.Omega_deg * pi / 180.0;
    double e = data.ecc;
    double mass = data.mass_earth * M_EARTH;

    // Положение в перицентре (θ = 0)
    double r_peri = a * (1 - e); // Расстояние от фокуса эллипса до перицентра
    x = r_peri * (cos(Omega) * cos(omega) - sin(Omega) * sin(omega) * cos(i));
    y = r_peri * (sin(Omega) * cos(omega) + cos(Omega) * sin(omega) * cos(i));
    z = r_peri * sin(omega) * sin(i);

    // Скорость в перицентре
    double mu = G * (M_star + mass); // Гравитационный параметр
    double v_peri = sqrt((mu / a) * (1 + e) / (1 - e)); // Скорость в перицентре
    vx = -v_peri * (cos(Omega) * sin(omega) + sin(Omega) * cos(omega) * cos(i));
    vy = -v_peri * (sin(Omega) * sin(omega) - cos(Omega) * cos(omega) * cos(i));
    vz = v_peri * cos(omega) * sin(i);
}

int main() {

    //std::vector<Body> solar_system_bodies = {
    //        {1.989e30, -8.572865039469250E+05 * 1000, -7.346088835335051E+05 * 1000, 2.685423265526889E+04 * 1000, 1.239859639798033E-02 * 1000, -6.348466611140617E-03 * 1000, -2.037876555553517E-04 * 1000}, // Солнце       
    //        {3.302e23, -5.879699189509091E+07 * 1000, -2.492820404148239E+07 * 1000, 3.364042452841429E+06 * 1000, 8.711982611106873E+00 * 1000, -4.284856986770977E+01 * 1000, -4.299279282370732E+00 * 1000}, // Меркурий        
    //        {4.867e24, 6.697319534635594E+07 * 1000, 8.337171945245868E+07 * 1000, -2.731933993919346E+06 * 1000, -2.735548307021769E+01 * 1000, 2.182743070706988E+01 * 1000, 1.878804135388283E+00 * 1000}, // Венера        
    //        {5.972e24, -2.758794880287251E+07 * 1000, 1.439239583084676E+08 * 1000, 1.921064327326417E+04 * 1000, -2.977686364628585E+01 * 1000, -5.535813340802556E+00 * 1000, -1.943387942073826E-04 * 1000}, // Земля        
    //        {6.417e23, -7.890038131682467E+07 * 1000, 2.274372361241295E+08 * 1000, 6.722196400986686E+06 * 1000, -2.199759485544059E+01 * 1000, -5.787405095467102E+00 * 1000, 4.184257990348734E-01 * 1000}, // Марс       
    //        {1.898e27, 1.571230833020991E+08 * 1000, 7.429840488421507E+08 * 1000, -6.597049828231782E+06 * 1000, -1.293244436609816E+01 * 1000, 3.325781476287804E+00 * 1000, 2.755437569190042E-01 * 1000}, // Юпитер 
    //        {5.683e26, 1.414498231862034E+09 * 1000, -2.647172137275474E+08 * 1000, -5.171551879510410E+07 * 1000, 1.240660798615463E+00 * 1000, 9.473546595187154E+00 * 1000, -2.135791731559418E-01 * 1000}, // Сатурн
    //        {8.681e25, 1.660221941282016E+09 * 1000, 2.406965914394302E+09 * 1000, -1.256903703858864E+07 * 1000, -5.655920544550284E+00 * 1000, 3.549247198028906E+00 * 1000, 8.651776992957205E-02 * 1000}, // Уран
    //        {1.024e26, 4.469116222588663E+09 * 1000, -9.560778256566879E+07 * 1000, -1.010264767638457E+08 * 1000, 8.064561683471368E-02 * 1000, 5.465730017544922E+00 * 1000, -1.151205185674022E-01 * 1000}, // Нептун
    //};


    //PlanetarySystem solar_system(solar_system_bodies);

    //
    //double t_end = 200 * year_to_sec; 
    //double dt = 12*3600;  // Шаг времени

    //
    //solar_system.SaveToFile_RK6("solar_system_orbits.txt", t_end, dt);

    //
    //std::vector<std::string> labels = { "Солнце", "Меркурий", "Венера", "Земля", "Марс", "Юпитер", "Сатурн", "Уран", "Нептун" };
    //std::vector<std::string> colors = { "yellow", "gray", "orange", "blue", "red", "brown", "purple", "black", "navy" };
    //solar_system.visual.PlotOrbits("solar_system_orbits.txt", labels, colors);
    //solar_system.visual.PlotEnergy();
    //solar_system.visual.PlotMomentum();


    const std::vector<NasaPlanetData> TOI_270 = {
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

    toi270_system.SaveToFile_RK6("toi270_orbits.txt", t_end, dt);

    std::vector<std::string> labels = { "TOI-270", "b", "c", "d" };
    std::vector<std::string> colors = { "red", "orange", "yellow", "green" };

    toi270_system.visual.PlotOrbits("toi270_orbits.txt", labels, colors);

    /*std::vector<double> x1 = { 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150};
    std::vector<double> y1 = { 2.7, 4.55, 5.8, 6.6, 7.3, 7.9, 8.4, 8.8, 9.2, 9.6, 10.0, 10.3, 10.5, 10.9, 11.1, 11.4, 11.7, 11.9, 12.1, 12.4, 12.6, 12.8, 12.945, 13.2, 13.4, 13.5, 13.7, 13.9, 14.1, 14.2, 14.4 };

    std::vector<double> y_log(y1.size());
    std::vector<double> x_log(x1.size());
    for (size_t i = 0; i < y1.size(); ++i) {
        y_log[i] = log10(y1[i]);
        x_log[i] = log10(x1[i]);
    }*/

    /*plt::plot(x_log, y_log);
    plt::title("Определение устойчивости планетной системы");
    plt::xlabel("M (в массах Солнца)");
    plt::ylabel("r (a.e.)");
    plt::legend();
    plt::grid(true);
    plt::axis("equal");
    plt::show();*/

    //toi270_system.visual.PlotOrbits("toi270_orbits.txt", labels, colors);

    return 0;
}