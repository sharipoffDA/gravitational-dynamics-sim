#include"algorithm.h"

// Вычисление ускорений
void Algorithm::Calculation_A(const std::vector<Body>& bodies, std::vector<double>& ax, std::vector<double>& ay, std::vector<double>& az) const {
    size_t n = bodies.size();
    ax.assign(n, 0.0);
    ay.assign(n, 0.0);
    az.assign(n, 0.0);

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            if (j != i) {
                double dx = bodies[j].x - bodies[i].x;
                double dy = bodies[j].y - bodies[i].y;
                double dz = bodies[j].z - bodies[i].z;
                double r = std::sqrt(dx * dx + dy * dy + dz * dz);
                double force = G * bodies[j].m / (r * r * r);
                ax[i] += force * dx;
                ay[i] += force * dy;
                az[i] += force * dz;
            }
        }
    }
}

// Рунге-Кутта 6 порядка
void Algorithm::RungeKutta6(std::vector<Body>& bodies, double dt) {

    size_t n = bodies.size();
    std::vector<double> ax(n), ay(n), az(n);
    std::vector<Body> k1 = bodies, k2 = bodies, k3 = bodies, k4 = bodies, k5 = bodies, k6 = bodies;

    // Шаг 1
    Calculation_A(bodies, ax, ay, az);
    for (size_t i = 0; i < n; ++i) {
        k1[i].vx = ax[i] * dt;
        k1[i].vy = ay[i] * dt;
        k1[i].vz = az[i] * dt;
        k1[i].x = bodies[i].vx * dt;
        k1[i].y = bodies[i].vy * dt;
        k1[i].z = bodies[i].vz * dt;
    }

    // Шаг 2
    for (size_t i = 0; i < n; ++i) {
        k2[i].x = bodies[i].x + 0.25 * k1[i].x;
        k2[i].y = bodies[i].y + 0.25 * k1[i].y;
        k2[i].z = bodies[i].z + 0.25 * k1[i].z;
        k2[i].vx = bodies[i].vx + 0.25 * k1[i].vx;
        k2[i].vy = bodies[i].vy + 0.25 * k1[i].vy;
        k2[i].vz = bodies[i].vz + 0.25 * k1[i].vz;
    }
    Calculation_A(k2, ax, ay, az);
    for (size_t i = 0; i < n; ++i) {
        k2[i].vx = ax[i] * dt;
        k2[i].vy = ay[i] * dt;
        k2[i].vz = az[i] * dt;
        k2[i].x = (bodies[i].vx + 0.25 * k1[i].vx) * dt;
        k2[i].y = (bodies[i].vy + 0.25 * k1[i].vy) * dt;
        k2[i].z = (bodies[i].vz + 0.25 * k1[i].vz) * dt;
    }

    // Шаг 3
    for (size_t i = 0; i < n; ++i) {
        k3[i].x = bodies[i].x + (3.0 / 32.0) * k1[i].x + (9.0 / 32.0) * k2[i].x;
        k3[i].y = bodies[i].y + (3.0 / 32.0) * k1[i].y + (9.0 / 32.0) * k2[i].y;
        k3[i].z = bodies[i].z + (3.0 / 32.0) * k1[i].z + (9.0 / 32.0) * k2[i].z;
        k3[i].vx = bodies[i].vx + (3.0 / 32.0) * k1[i].vx + (9.0 / 32.0) * k2[i].vx;
        k3[i].vy = bodies[i].vy + (3.0 / 32.0) * k1[i].vy + (9.0 / 32.0) * k2[i].vy;
        k3[i].vz = bodies[i].vz + (3.0 / 32.0) * k1[i].vz + (9.0 / 32.0) * k2[i].vz;
    }
    Calculation_A(k3, ax, ay, az);
    for (size_t i = 0; i < n; ++i) {
        k3[i].vx = ax[i] * dt;
        k3[i].vy = ay[i] * dt;
        k3[i].vz = az[i] * dt;
        k3[i].x = (bodies[i].vx + (3.0 / 32.0) * k1[i].vx + (9.0 / 32.0) * k2[i].vx) * dt;
        k3[i].y = (bodies[i].vy + (3.0 / 32.0) * k1[i].vy + (9.0 / 32.0) * k2[i].vy) * dt;
        k3[i].z = (bodies[i].vz + (3.0 / 32.0) * k1[i].vz + (9.0 / 32.0) * k2[i].vz) * dt;
    }

    // Шаг 4
    for (size_t i = 0; i < n; ++i) {
        k4[i].x = bodies[i].x + (1932.0 / 2197.0) * k1[i].x - (7200.0 / 2197.0) * k2[i].x + (7296.0 / 2197.0) * k3[i].x;
        k4[i].y = bodies[i].y + (1932.0 / 2197.0) * k1[i].y - (7200.0 / 2197.0) * k2[i].y + (7296.0 / 2197.0) * k3[i].y;
        k4[i].z = bodies[i].z + (1932.0 / 2197.0) * k1[i].z - (7200.0 / 2197.0) * k2[i].z + (7296.0 / 2197.0) * k3[i].z;
        k4[i].vx = bodies[i].vx + (1932.0 / 2197.0) * k1[i].vx - (7200.0 / 2197.0) * k2[i].vx + (7296.0 / 2197.0) * k3[i].vx;
        k4[i].vy = bodies[i].vy + (1932.0 / 2197.0) * k1[i].vy - (7200.0 / 2197.0) * k2[i].vy + (7296.0 / 2197.0) * k3[i].vy;
        k4[i].vz = bodies[i].vz + (1932.0 / 2197.0) * k1[i].vz - (7200.0 / 2197.0) * k2[i].vz + (7296.0 / 2197.0) * k3[i].vz;
    }
    Calculation_A(k4, ax, ay, az);
    for (size_t i = 0; i < n; ++i) {
        k4[i].vx = ax[i] * dt;
        k4[i].vy = ay[i] * dt;
        k4[i].vz = az[i] * dt;
        k4[i].x = (bodies[i].vx + (1932.0 / 2197.0) * k1[i].vx - (7200.0 / 2197.0) * k2[i].vx + (7296.0 / 2197.0) * k3[i].vx) * dt;
        k4[i].y = (bodies[i].vy + (1932.0 / 2197.0) * k1[i].vy - (7200.0 / 2197.0) * k2[i].vy + (7296.0 / 2197.0) * k3[i].vy) * dt;
        k4[i].z = (bodies[i].vz + (1932.0 / 2197.0) * k1[i].vz - (7200.0 / 2197.0) * k2[i].vz + (7296.0 / 2197.0) * k3[i].vz) * dt;
    }

    // Шаг 5
    for (size_t i = 0; i < n; ++i) {
        k5[i].x = bodies[i].x + (439.0 / 216.0) * k1[i].x - 8.0 * k2[i].x + (3680.0 / 513.0) * k3[i].x - (845.0 / 4104.0) * k4[i].x;
        k5[i].y = bodies[i].y + (439.0 / 216.0) * k1[i].y - 8.0 * k2[i].y + (3680.0 / 513.0) * k3[i].y - (845.0 / 4104.0) * k4[i].y;
        k5[i].z = bodies[i].z + (439.0 / 216.0) * k1[i].z - 8.0 * k2[i].z + (3680.0 / 513.0) * k3[i].z - (845.0 / 4104.0) * k4[i].z;
        k5[i].vx = bodies[i].vx + (439.0 / 216.0) * k1[i].vx - 8.0 * k2[i].vx + (3680.0 / 513.0) * k3[i].vx - (845.0 / 4104.0) * k4[i].vx;
        k5[i].vy = bodies[i].vy + (439.0 / 216.0) * k1[i].vy - 8.0 * k2[i].vy + (3680.0 / 513.0) * k3[i].vy - (845.0 / 4104.0) * k4[i].vy;
        k5[i].vz = bodies[i].vz + (439.0 / 216.0) * k1[i].vz - 8.0 * k2[i].vz + (3680.0 / 513.0) * k3[i].vz - (845.0 / 4104.0) * k4[i].vz;
    }
    Calculation_A(k5, ax, ay, az);
    for (size_t i = 0; i < n; ++i) {
        k5[i].vx = ax[i] * dt;
        k5[i].vy = ay[i] * dt;
        k5[i].vz = az[i] * dt;
        k5[i].x = (bodies[i].vx + (439.0 / 216.0) * k1[i].vx - 8.0 * k2[i].vx + (3680.0 / 513.0) * k3[i].vx - (845.0 / 4104.0) * k4[i].vx) * dt;
        k5[i].y = (bodies[i].vy + (439.0 / 216.0) * k1[i].vy - 8.0 * k2[i].vy + (3680.0 / 513.0) * k3[i].vy - (845.0 / 4104.0) * k4[i].vy) * dt;
        k5[i].z = (bodies[i].vz + (439.0 / 216.0) * k1[i].vz - 8.0 * k2[i].vz + (3680.0 / 513.0) * k3[i].vz - (845.0 / 4104.0) * k4[i].vz) * dt;
    }

    // Шаг 6
    for (size_t i = 0; i < n; ++i) {
        k6[i].x = bodies[i].x - (8.0 / 27.0) * k1[i].x + 2.0 * k2[i].x - (3544.0 / 2565.0) * k3[i].x + (1859.0 / 4104.0) * k4[i].x - (11.0 / 40.0) * k5[i].x;
        k6[i].y = bodies[i].y - (8.0 / 27.0) * k1[i].y + 2.0 * k2[i].y - (3544.0 / 2565.0) * k3[i].y + (1859.0 / 4104.0) * k4[i].y - (11.0 / 40.0) * k5[i].y;
        k6[i].z = bodies[i].z - (8.0 / 27.0) * k1[i].z + 2.0 * k2[i].z - (3544.0 / 2565.0) * k3[i].z + (1859.0 / 4104.0) * k4[i].z - (11.0 / 40.0) * k5[i].z;
        k6[i].vx = bodies[i].vx - (8.0 / 27.0) * k1[i].vx + 2.0 * k2[i].vx - (3544.0 / 2565.0) * k3[i].vx + (1859.0 / 4104.0) * k4[i].vx - (11.0 / 40.0) * k5[i].vx;
        k6[i].vy = bodies[i].vy - (8.0 / 27.0) * k1[i].vy + 2.0 * k2[i].vy - (3544.0 / 2565.0) * k3[i].vy + (1859.0 / 4104.0) * k4[i].vy - (11.0 / 40.0) * k5[i].vy;
        k6[i].vz = bodies[i].vz - (8.0 / 27.0) * k1[i].vz + 2.0 * k2[i].vz - (3544.0 / 2565.0) * k3[i].vz + (1859.0 / 4104.0) * k4[i].vz - (11.0 / 40.0) * k5[i].vz;
    }
    Calculation_A(k6, ax, ay, az);
    for (size_t i = 0; i < n; ++i) {
        k6[i].vx = ax[i] * dt;
        k6[i].vy = ay[i] * dt;
        k6[i].vz = az[i] * dt;
        k6[i].x = (bodies[i].vx - (8.0 / 27.0) * k1[i].vx + 2.0 * k2[i].vx - (3544.0 / 2565.0) * k3[i].vx + (1859.0 / 4104.0) * k4[i].vx - (11.0 / 40.0) * k5[i].vx) * dt;
        k6[i].y = (bodies[i].vy - (8.0 / 27.0) * k1[i].vy + 2.0 * k2[i].vy - (3544.0 / 2565.0) * k3[i].vy + (1859.0 / 4104.0) * k4[i].vy - (11.0 / 40.0) * k5[i].vy) * dt;
        k6[i].z = (bodies[i].vz - (8.0 / 27.0) * k1[i].vz + 2.0 * k2[i].vz - (3544.0 / 2565.0) * k3[i].vz + (1859.0 / 4104.0) * k4[i].vz - (11.0 / 40.0) * k5[i].vz) * dt;
    }


    for (size_t i = 0; i < n; ++i) {
        bodies[i].x += (16.0 / 135.0) * k1[i].x + (6656.0 / 12825.0) * k3[i].x + (28561.0 / 56430.0) * k4[i].x - (9.0 / 50.0) * k5[i].x + (2.0 / 55.0) * k6[i].x;
        bodies[i].y += (16.0 / 135.0) * k1[i].y + (6656.0 / 12825.0) * k3[i].y + (28561.0 / 56430.0) * k4[i].y - (9.0 / 50.0) * k5[i].y + (2.0 / 55.0) * k6[i].y;
        bodies[i].z += (16.0 / 135.0) * k1[i].z + (6656.0 / 12825.0) * k3[i].z + (28561.0 / 56430.0) * k4[i].z - (9.0 / 50.0) * k5[i].z + (2.0 / 55.0) * k6[i].z;
        bodies[i].vx += (16.0 / 135.0) * k1[i].vx + (6656.0 / 12825.0) * k3[i].vx + (28561.0 / 56430.0) * k4[i].vx - (9.0 / 50.0) * k5[i].vx + (2.0 / 55.0) * k6[i].vx;
        bodies[i].vy += (16.0 / 135.0) * k1[i].vy + (6656.0 / 12825.0) * k3[i].vy + (28561.0 / 56430.0) * k4[i].vy - (9.0 / 50.0) * k5[i].vy + (2.0 / 55.0) * k6[i].vy;
        bodies[i].vz += (16.0 / 135.0) * k1[i].vz + (6656.0 / 12825.0) * k3[i].vz + (28561.0 / 56430.0) * k4[i].vz - (9.0 / 50.0) * k5[i].vz + (2.0 / 55.0) * k6[i].vz;
    }
}

// Рунге-Кутта 4 порядка
void Algorithm::RungeKutta4(std::vector<Body>& bodies, double dt) {
    size_t n = bodies.size();
    std::vector<double> ax(n), ay(n), az(n);
    std::vector<Body> k1 = bodies, k2 = bodies, k3 = bodies, k4 = bodies;

    // Шаг 1
    Calculation_A(bodies, ax, ay, az);
    for (size_t i = 0; i < n; ++i) {
        k1[i].vx = ax[i] * dt;
        k1[i].vy = ay[i] * dt;
        k1[i].vz = az[i] * dt;
        k1[i].x = bodies[i].vx * dt;
        k1[i].y = bodies[i].vy * dt;
        k1[i].z = bodies[i].vz * dt;
    }

    // Шаг 2
    for (size_t i = 0; i < n; ++i) {
        k2[i].x = bodies[i].x + 0.5 * k1[i].x;
        k2[i].y = bodies[i].y + 0.5 * k1[i].y;
        k2[i].z = bodies[i].z + 0.5 * k1[i].z;
        k2[i].vx = bodies[i].vx + 0.5 * k1[i].vx;
        k2[i].vy = bodies[i].vy + 0.5 * k1[i].vy;
        k2[i].vz = bodies[i].vz + 0.5 * k1[i].vz;
    }
    Calculation_A(k2, ax, ay, az);
    for (size_t i = 0; i < n; ++i) {
        k2[i].vx = ax[i] * dt;
        k2[i].vy = ay[i] * dt;
        k2[i].vz = az[i] * dt;
        k2[i].x = (bodies[i].vx + 0.5 * k1[i].vx) * dt;
        k2[i].y = (bodies[i].vy + 0.5 * k1[i].vy) * dt;
        k2[i].z = (bodies[i].vz + 0.5 * k1[i].vz) * dt;
    }

    // Шаг 3
    for (size_t i = 0; i < n; ++i) {
        k3[i].x = bodies[i].x + 0.5 * k2[i].x;
        k3[i].y = bodies[i].y + 0.5 * k2[i].y;
        k3[i].z = bodies[i].z + 0.5 * k2[i].z;
        k3[i].vx = bodies[i].vx + 0.5 * k2[i].vx;
        k3[i].vy = bodies[i].vy + 0.5 * k2[i].vy;
        k3[i].vz = bodies[i].vz + 0.5 * k2[i].vz;
    }
    Calculation_A(k3, ax, ay, az);
    for (size_t i = 0; i < n; ++i) {
        k3[i].vx = ax[i] * dt;
        k3[i].vy = ay[i] * dt;
        k3[i].vz = az[i] * dt;
        k3[i].x = (bodies[i].vx + 0.5 * k2[i].vx) * dt;
        k3[i].y = (bodies[i].vy + 0.5 * k2[i].vy) * dt;
        k3[i].z = (bodies[i].vz + 0.5 * k2[i].vz) * dt;
    }

    // Шаг 4
    for (size_t i = 0; i < n; ++i) {
        k4[i].x = bodies[i].x + k3[i].x;
        k4[i].y = bodies[i].y + k3[i].y;
        k4[i].z = bodies[i].z + k3[i].z;
        k4[i].vx = bodies[i].vx + k3[i].vx;
        k4[i].vy = bodies[i].vy + k3[i].vy;
        k4[i].vz = bodies[i].vz + k3[i].vz;
    }
    Calculation_A(k4, ax, ay, az);
    for (size_t i = 0; i < n; ++i) {
        k4[i].vx = ax[i] * dt;
        k4[i].vy = ay[i] * dt;
        k4[i].vz = az[i] * dt;
        k4[i].x = (bodies[i].vx + k3[i].vx) * dt;
        k4[i].y = (bodies[i].vy + k3[i].vy) * dt;
        k4[i].z = (bodies[i].vz + k3[i].vz) * dt;
    }


    for (size_t i = 0; i < n; ++i) {
        bodies[i].x += (k1[i].x + 2 * k2[i].x + 2 * k3[i].x + k4[i].x) / 6;
        bodies[i].y += (k1[i].y + 2 * k2[i].y + 2 * k3[i].y + k4[i].y) / 6;
        bodies[i].z += (k1[i].z + 2 * k2[i].z + 2 * k3[i].z + k4[i].z) / 6;
        bodies[i].vx += (k1[i].vx + 2 * k2[i].vx + 2 * k3[i].vx + k4[i].vx) / 6;
        bodies[i].vy += (k1[i].vy + 2 * k2[i].vy + 2 * k3[i].vy + k4[i].vy) / 6;
        bodies[i].vz += (k1[i].vz + 2 * k2[i].vz + 2 * k3[i].vz + k4[i].vz) / 6;
    }
}

// Рунге-Кутта 2 порядка
void Algorithm::RungeKutta2(std::vector<Body>& bodies, double dt) {
    size_t n = bodies.size();
    std::vector<double> ax(n), ay(n), az(n);
    std::vector<Body> k1 = bodies, k2 = bodies;

    // Шаг 1
    Calculation_A(bodies, ax, ay, az);
    for (size_t i = 0; i < n; ++i) {
        k1[i].vx = ax[i] * dt;
        k1[i].vy = ay[i] * dt;
        k1[i].vz = az[i] * dt;
        k1[i].x = bodies[i].vx * dt;
        k1[i].y = bodies[i].vy * dt;
        k1[i].z = bodies[i].vz * dt;
    }

    // Шаг 2
    for (size_t i = 0; i < n; ++i) {
        k2[i].x = bodies[i].x + 0.5 * k1[i].x;
        k2[i].y = bodies[i].y + 0.5 * k1[i].y;
        k2[i].z = bodies[i].z + 0.5 * k1[i].z;
        k2[i].vx = bodies[i].vx + 0.5 * k1[i].vx;
        k2[i].vy = bodies[i].vy + 0.5 * k1[i].vy;
        k2[i].vz = bodies[i].vz + 0.5 * k1[i].vz;
    }
    Calculation_A(k2, ax, ay, az);
    for (size_t i = 0; i < n; ++i) {
        k2[i].vx = ax[i] * dt;
        k2[i].vy = ay[i] * dt;
        k2[i].vz = az[i] * dt;
        k2[i].x = (bodies[i].vx + 0.5 * k1[i].vx) * dt;
        k2[i].y = (bodies[i].vy + 0.5 * k1[i].vy) * dt;
        k2[i].z = (bodies[i].vz + 0.5 * k1[i].vz) * dt;
    }


    for (size_t i = 0; i < n; ++i) {
        bodies[i].x += k2[i].x;
        bodies[i].y += k2[i].y;
        bodies[i].z += k2[i].z;
        bodies[i].vx += k2[i].vx;
        bodies[i].vy += k2[i].vy;
        bodies[i].vz += k2[i].vz;
    }
}

void Algorithm::EulerMethod(std::vector<Body>& bodies, double dt) {
    size_t n = bodies.size();
    std::vector<double> ax(n), ay(n), az(n);

    Calculation_A(bodies, ax, ay, az);

    for (size_t i = 0; i < n; ++i) {

        bodies[i].vx += ax[i] * dt;
        bodies[i].vy += ay[i] * dt;
        bodies[i].vz += az[i] * dt;

        bodies[i].x += bodies[i].vx * dt;
        bodies[i].y += bodies[i].vy * dt;
        bodies[i].z += bodies[i].vz * dt;
    }
}

void Algorithm::Leapfrog(std::vector<Body>& bodies, double dt)
{
    size_t n = bodies.size();

    std::vector<double> ax(n), ay(n), az(n);
    Calculation_A(bodies, ax, ay, az);

    Calculation_A(bodies, ax, ay, az);

    for (size_t i = 0; i < n; ++i) {
        bodies[i].vx += 0.5 * ax[i] * dt;
        bodies[i].vy += 0.5 * ay[i] * dt;
        bodies[i].vz += 0.5 * az[i] * dt;
    }

    for (size_t i = 0; i < n; ++i) {
        bodies[i].x = bodies[i].x + bodies[i].vx * dt;
        bodies[i].y = bodies[i].y + bodies[i].vy * dt;
        bodies[i].z = bodies[i].z + bodies[i].vz * dt;
    }

    Calculation_A(bodies, ax, ay, az);
    for (size_t i = 0; i < n; ++i) {
        bodies[i].vx += 0.5 * ax[i] * dt;
        bodies[i].vy += 0.5 * ay[i] * dt;
        bodies[i].vz += 0.5 * az[i] * dt;
    }
}