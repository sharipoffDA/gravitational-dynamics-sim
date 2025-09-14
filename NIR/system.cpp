#include"system.h"

PlanetarySystem::PlanetarySystem(const std::vector<Body>& initial_bodies) : bodies(initial_bodies), visual(bodies, time_values, energy_values, momentum_values)
{
    orbit_params.resize(bodies.size(), { std::numeric_limits<double>::max(), 0.0, 0.0 });
}

// Вычисление кинетической энергии
double PlanetarySystem::KineticEnergy() const {
    double KE = 0.0;
    for (const auto& body : bodies) {
        KE += 0.5 * body.m * (body.vx * body.vx + body.vy * body.vy + body.vz * body.vz);
    }
    return KE;
}

double PlanetarySystem::Momentum() const {
    double P = 0.0;
    for (const auto& body : bodies) {
        P += body.m * sqrt((body.vx * body.vx + body.vy * body.vy + body.vz * body.vz));
    }
    return P;
}

// Вычисление потенциальной энергии
double PlanetarySystem::PotentialEnergy() const {
    double PE = 0.0;
    size_t n = bodies.size();
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            double dx = bodies[j].x - bodies[i].x;
            double dy = bodies[j].y - bodies[i].y;
            double dz = bodies[j].z - bodies[i].z;
            double r = sqrt(dx * dx + dy * dy + dz * dz);
            PE -= G * bodies[i].m * bodies[j].m / r;
        }
    }
    return PE;
}

void PlanetarySystem::SaveToFile_Leapfrog(const std::string& filename, double t_end, double dt) {
    std::ofstream out(filename);
    for (double t = 0; t < t_end; t += dt) {

        recenter();

        algorithm.Leapfrog(bodies, dt);

        out << t / year_to_sec << " ";
        for (const auto& body : bodies) {
            out << body.x << " " << body.y << " " << body.z << " ";
        }
        out << "\n";

        double total_energy = KineticEnergy() + PotentialEnergy();

        time_values.push_back(t / year_to_sec);
        energy_values.push_back(total_energy);
        momentum_values.push_back(Momentum());

        for (size_t i = 0; i < bodies.size(); ++i) {

            double dx = bodies[i].x - bodies[0].x;
            double dy = bodies[i].y - bodies[0].y;
            double dz = bodies[i].z - bodies[0].z;
            double r = std::sqrt(dx * dx + dy * dy + dz * dz);

            if (r < orbit_params[i].perihelion) {
                orbit_params[i].perihelion = r;
            }
            if (r > orbit_params[i].aphelion) {
                orbit_params[i].aphelion = r;
            }
        }

    }
    out.close();
    std::cout << "Результаты записаны в " << filename << ".\n";

    for (size_t i = 0; i < bodies.size(); ++i) {
        orbit_params[i].eccentricity = (orbit_params[i].aphelion - orbit_params[i].perihelion) / (orbit_params[i].aphelion + orbit_params[i].perihelion);
    }

    for (size_t i = 0; i < bodies.size(); ++i) {
        orbit_params[i].semimajor_axis = (orbit_params[i].aphelion + orbit_params[i].perihelion) * 0.5;
    }

    visual = Visual(bodies, time_values, energy_values, momentum_values);

    SaveOrbitParameters("orbit_parameters.txt");
}


// Запись параметров орбиты в файл
void PlanetarySystem::SaveOrbitParameters(const std::string& filename) const {
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Ошибка: не удалось открыть файл " << filename << " для записи.\n";
        return;
    }

    out << "Планета\tПеригелий (м)\tАфелий (м)\tЭксцентриситет\tБольшая полуось";
    for (size_t i = 0; i < bodies.size(); ++i) {
        out << "Тело " << i + 1 << "\t"
            << orbit_params[i].perihelion << "\t"
            << orbit_params[i].aphelion << "\t"
            << orbit_params[i].eccentricity << "\t"
            << orbit_params[i].semimajor_axis << "\n";
    }
    out.close();
    std::cout << "Параметры орбит записаны в " << filename << ".\n";
}

// Запись данных в файл и расчёт с помощью метода Эйлера
void PlanetarySystem::SaveToFile_Euler(const std::string& filename, double t_end, double dt) {
    std::ofstream out(filename);
    for (double t = 0; t < t_end; t += dt) {

        recenter();

        algorithm.EulerMethod(bodies, dt);

        out << t / year_to_sec << " ";
        for (const auto& body : bodies) {
            out << body.x << " " << body.y << " " << body.z << " ";
        }
        out << "\n";

        double total_energy = KineticEnergy() + PotentialEnergy();

        time_values.push_back(t / year_to_sec);
        energy_values.push_back(total_energy);
        momentum_values.push_back(Momentum());

        for (size_t i = 0; i < bodies.size(); ++i) {

            double dx = bodies[i].x - bodies[0].x;
            double dy = bodies[i].y - bodies[0].y;
            double dz = bodies[i].z - bodies[0].z;
            double r = std::sqrt(dx * dx + dy * dy + dz * dz);

            if (r < orbit_params[i].perihelion) {
                orbit_params[i].perihelion = r;
            }
            if (r > orbit_params[i].aphelion) {
                orbit_params[i].aphelion = r;
            }
        }
    }
    out.close();
    std::cout << "Результаты записаны в " << filename << ".\n";

    for (size_t i = 0; i < bodies.size(); ++i) {
        orbit_params[i].eccentricity = (orbit_params[i].aphelion - orbit_params[i].perihelion) / (orbit_params[i].aphelion + orbit_params[i].perihelion);
    }

    for (size_t i = 0; i < bodies.size(); ++i) {
        orbit_params[i].semimajor_axis = (orbit_params[i].aphelion + orbit_params[i].perihelion) * 0.5;
    }

    visual = Visual(bodies, time_values, energy_values, momentum_values);

    SaveOrbitParameters("orbit_parameters.txt");

}

// Запись данных в файл и расчёт с помощью метода Рунге-Кутта 2 порядка
void PlanetarySystem::SaveToFile_RK2(const std::string& filename, double t_end, double dt) {
    std::ofstream out(filename);
    for (double t = 0; t < t_end; t += dt) {

        recenter();

        algorithm.RungeKutta2(bodies, dt);

        out << t / year_to_sec << " ";
        for (const auto& body : bodies) {
            out << body.x << " " << body.y << " " << body.z << " ";
        }
        out << "\n";

        double total_energy = KineticEnergy() + PotentialEnergy();

        time_values.push_back(t / year_to_sec);
        energy_values.push_back(total_energy);
        momentum_values.push_back(Momentum());

        for (size_t i = 0; i < bodies.size(); ++i) {

            double dx = bodies[i].x - bodies[0].x;
            double dy = bodies[i].y - bodies[0].y;
            double dz = bodies[i].z - bodies[0].z;
            double r = std::sqrt(dx * dx + dy * dy + dz * dz);

            if (r < orbit_params[i].perihelion) {
                orbit_params[i].perihelion = r;
            }
            if (r > orbit_params[i].aphelion) {
                orbit_params[i].aphelion = r;
            }
        }


    }
    out.close();
    std::cout << "Результаты записаны в " << filename << ".\n";


    for (size_t i = 0; i < bodies.size(); ++i) {
        orbit_params[i].eccentricity = (orbit_params[i].aphelion - orbit_params[i].perihelion) / (orbit_params[i].aphelion + orbit_params[i].perihelion);
    }

    for (size_t i = 0; i < bodies.size(); ++i) {
        orbit_params[i].semimajor_axis = (orbit_params[i].aphelion + orbit_params[i].perihelion) * 0.5;
    }

    visual = Visual(bodies, time_values, energy_values, momentum_values);

    SaveOrbitParameters("orbit_parameters.txt");
}

// Запись данных в файл и расчёт с помощью метода Рунге-Кутта 4 порядка
void PlanetarySystem::SaveToFile_RK4(const std::string& filename, double t_end, double dt) {
    std::ofstream out(filename);
    int cnt = 0;
    for (double t = 0; t < t_end; t += dt) {

        recenter();

        algorithm.RungeKutta4(bodies, dt);

        out << t / year_to_sec << " ";
        for (const auto& body : bodies) {
            out << body.x << " " << body.y << " " << body.z << " ";
        }
        out << "\n";

        double total_energy = KineticEnergy() + PotentialEnergy();

        time_values.push_back(t / year_to_sec);
        energy_values.push_back(total_energy);
        momentum_values.push_back(Momentum());

        cnt++;

        for (size_t i = 0; i < bodies.size(); ++i) {

            double dx = bodies[i].x - bodies[0].x;
            double dy = bodies[i].y - bodies[0].y;
            double dz = bodies[i].z - bodies[0].z;
            double r = std::sqrt(dx * dx + dy * dy + dz * dz);

            if (r < orbit_params[i].perihelion) {
                orbit_params[i].perihelion = r;
            }
            if (r > orbit_params[i].aphelion) {
                orbit_params[i].aphelion = r;
            }
        }
    }
    out.close();
    std::cout << "Результаты записаны в " << filename << ".\n";

    for (size_t i = 0; i < bodies.size(); ++i) {
        orbit_params[i].eccentricity = (orbit_params[i].aphelion - orbit_params[i].perihelion) / (orbit_params[i].aphelion + orbit_params[i].perihelion);
    }

    for (size_t i = 0; i < bodies.size(); ++i) {
        orbit_params[i].semimajor_axis = (orbit_params[i].aphelion + orbit_params[i].perihelion) * 0.5;
    }

    std::cout << "\n" << cnt << std::endl;

    visual = Visual(bodies, time_values, energy_values, momentum_values);

    SaveOrbitParameters("orbit_parameters.txt");
}

void PlanetarySystem::SaveToFile_RK6(const std::string& filename, double t_end, double dt) {
    std::ofstream out(filename);

    std::vector<std::vector<OrbitState>> orbit_history(bodies.size());

    std::vector<double> initial_distances(bodies.size(), 0.0);
    std::vector<bool> is_unstable(bodies.size(), false);

    for (double t = 0; t < t_end; t += dt) {

        recenter();

        algorithm.RungeKutta6(bodies, dt);

        out << t / year_to_sec << " ";
        for (const auto& body : bodies) {
            out << body.x << " " << body.y << " " << body.z << " ";
        }
        out << "\n";

        time_values.push_back(t / year_to_sec);
        energy_values.push_back(KineticEnergy() + PotentialEnergy());
        momentum_values.push_back(Momentum());

        for (size_t i = 0; i < bodies.size(); ++i) {
            double dx = bodies[i].x - bodies[0].x;
            double dy = bodies[i].y - bodies[0].y;
            double dz = bodies[i].z - bodies[0].z;
            double r = std::sqrt(dx * dx + dy * dy + dz * dz);


            if (r < orbit_params[i].perihelion) orbit_params[i].perihelion = r;
            if (r > orbit_params[i].aphelion) orbit_params[i].aphelion = r;


            OrbitState state;
            state.time = t;
            state.distance = r;
            orbit_history[i].push_back(state);


            if (t > 0) {
                double rel_change = fabs(r - initial_distances[i]) / initial_distances[i];
                if (rel_change > 0.1) {
                    is_unstable[i] = true;
                }
            }
            else {
                initial_distances[i] = r;
            }
        }
    }
    out.close();

    // Расчет эксцентриситета и большой полуоси
    for (size_t i = 0; i < bodies.size(); ++i) {
        orbit_params[i].eccentricity = (orbit_params[i].aphelion - orbit_params[i].perihelion) /
            (orbit_params[i].aphelion + orbit_params[i].perihelion);
        orbit_params[i].semimajor_axis = (orbit_params[i].aphelion + orbit_params[i].perihelion) * 0.5;
    }

    visual = Visual(bodies, time_values, energy_values, momentum_values);

    // Сохранение параметров орбиты
    SaveOrbitParameters("orbit_parameters.txt");

    // Анализ и вывод информации об устойчивости
    std::cout << "\nАнализ устойчивости системы:\n";
    for (size_t i = 1; i < bodies.size(); ++i) { // начинаем с 1, так как 0 - это звезда
        std::cout << "Тело " << i << ": ";
        if (is_unstable[i]) {
            std::cout << "НЕУСТОЙЧИВАЯ орбита (отклонение > 10%)\n";
        }
        else {
            std::cout << "устойчивая орбита\n";
        }
        std::cout << "  Начальное расстояние: " << initial_distances[i] / AU << " а.е.\n";
        std::cout << "  Эксцентриситет: " << orbit_params[i].eccentricity << "\n";
        std::cout << "  Большая полуось: " << orbit_params[i].semimajor_axis / AU << " а.е.\n";
    }
}

void PlanetarySystem::SaveToFileOrbitHistory(const std::vector<std::vector<OrbitState>>& orbit_history) {
    std::ofstream stability_out("stability_analysis.txt");
    stability_out << "# Анализ устойчивости орбит\n";
    stability_out << "# Тело Время(лет) Расстояние(а.е.) Относительное_изменение\n";

    for (size_t i = 1; i < orbit_history.size(); ++i) {
        double initial_r = orbit_history[i][0].distance;

        for (const auto& state : orbit_history[i]) {
            double rel_change = (state.distance - initial_r) / initial_r;
            stability_out << i << " " << state.time / year_to_sec << " "
                << state.distance / AU << " " << rel_change << "\n";
        }
    }
    stability_out.close();

    std::cout << "Детальный анализ устойчивости сохранен в stability_analysis.txt\n";
}


void PlanetarySystem::recenter() {
    size_t star_idx = 0;

    double star_x = bodies[star_idx].x;
    double star_y = bodies[star_idx].y;
    double star_z = bodies[star_idx].z;
    double star_vx = bodies[star_idx].vx;
    double star_vy = bodies[star_idx].vy;
    double star_vz = bodies[star_idx].vz;

    for (size_t i = 0; i < bodies.size(); ++i) {
        if (i != star_idx) {
            bodies[i].x -= star_x;
            bodies[i].y -= star_y;
            bodies[i].z -= star_z;
            bodies[i].vx -= star_vx;
            bodies[i].vy -= star_vy;
            bodies[i].vz -= star_vz;
        }
    }

    bodies[star_idx].x = 0;
    bodies[star_idx].y = 0;
    bodies[star_idx].z = 0;
    bodies[star_idx].vx = 0;
    bodies[star_idx].vy = 0;
    bodies[star_idx].vz = 0;
}