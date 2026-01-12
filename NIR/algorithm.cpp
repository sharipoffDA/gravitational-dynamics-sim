#include"algorithm.h"

// Также правильная реализация (как в статье)
//void Algorithm::Calculation_A(const std::vector<Body>& bodies,
//    std::vector<double>& ax, std::vector<double>& ay, std::vector<double>& az,
//    bool use_relativistic) const {
//    size_t n = bodies.size();
//    ax.assign(n, 0.0);
//    ay.assign(n, 0.0);
//    az.assign(n, 0.0);
//
//    for (size_t i = 0; i < n; ++i) {
//        for (size_t j = i + 1; j < n; ++j) {
//            double dx = bodies[j].x - bodies[i].x;
//            double dy = bodies[j].y - bodies[i].y;
//            double dz = bodies[j].z - bodies[i].z;
//            double r = std::sqrt(dx * dx + dy * dy + dz * dz);
//            double r3 = r * r * r;
//
//            double force = G * bodies[i].m * bodies[j].m / r3;
//
//            ax[i] += force * dx;
//            ay[i] += force * dy;
//            az[i] += force * dz;
//
//            ax[j] -= force * dx;
//            ay[j] -= force * dy;
//            az[j] -= force * dz;
//        }
//
//        ax[i] /= bodies[i].m;
//        ay[i] /= bodies[i].m;
//        az[i] /= bodies[i].m;
//    }
//
//    if (use_relativistic) {
//
//        for (size_t i = 0; i < n; ++i) {
//            double pn_ax = 0.0, pn_ay = 0.0, pn_az = 0.0;
//
//            for (size_t j = 0; j < n; ++j) {
//                if (j == i) continue;
//
//                // Относительные положение и скорость (r_ij, v_ij)
//                double rx = bodies[j].x - bodies[i].x;
//                double ry = bodies[j].y - bodies[i].y;
//                double rz = bodies[j].z - bodies[i].z;
//
//                double vx = bodies[j].vx - bodies[i].vx;
//                double vy = bodies[j].vy - bodies[i].vy;
//                double vz = bodies[j].vz - bodies[i].vz;
//
//                double r = std::sqrt(rx * rx + ry * ry + rz * rz);
//                double v_sq = vx * vx + vy * vy + vz * vz;
//                double dot_rv = rx * vx + ry * vy + rz * vz;  // v_ij · r_ij
//
//                double common_factor = (G * bodies[j].m) / (c * c * r * r * r);
//
//                double a_pn_x = common_factor * (
//                    (4.0 * G * bodies[j].m / r - v_sq) * rx + 4.0 * dot_rv * vx);
//
//                double a_pn_y = common_factor * (
//                    (4.0 * G * bodies[j].m / r - v_sq) * ry + 4.0 * dot_rv * vy);
//
//                double a_pn_z = common_factor * (
//                    (4.0 * G * bodies[j].m / r - v_sq) * rz + 4.0 * dot_rv * vz);
//
//                pn_ax += a_pn_x;
//                pn_ay += a_pn_y;
//                pn_az += a_pn_z;
//            }
//
//            ax[i] += pn_ax;
//            ay[i] += pn_ay;
//            az[i] += pn_az;
//        }
//    }
//}

void Algorithm::Calculation_A(const std::vector<Body>& bodies, std::vector<double>& ax,
    std::vector<double>& ay, std::vector<double>& az,
    bool use_relativistic) const {
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

                double a_newton = G * bodies[j].m / (r * r);
                ax[i] += a_newton * dx / r;
                ay[i] += a_newton * dy / r;
                az[i] += a_newton * dz / r;
            }
        }
    }

    if (use_relativistic) {
        // Копируем ньютоновские ускорения для приближения в 1PN-терминах
        std::vector<double> ax_newt = ax;
        std::vector<double> ay_newt = ay;
        std::vector<double> az_newt = az;

        // Постньютоновская поправка 1PN (EIH) для всех тел
        for (size_t i = 0; i < n; ++i) {  // тело A = i
            double pn_ax = 0.0, pn_ay = 0.0, pn_az = 0.0;

            // Третий термин: 7/2 * sum_B G m_B a_B / r_AB (векторный)
            double sum3_x = 0.0, sum3_y = 0.0, sum3_z = 0.0;
            for (size_t j = 0; j < n; ++j) {
                if (j == i) continue;
                double dx = bodies[j].x - bodies[i].x;
                double dy = bodies[j].y - bodies[i].y;
                double dz = bodies[j].z - bodies[i].z;
                double r = std::sqrt(dx * dx + dy * dy + dz * dz);

                double factor = G * bodies[j].m / r;
                sum3_x += factor * ax_newt[j];
                sum3_y += factor * ay_newt[j];
                sum3_z += factor * az_newt[j];
            }
            pn_ax += 3.5 * sum3_x;  // 7/2 = 3.5
            pn_ay += 3.5 * sum3_y;
            pn_az += 3.5 * sum3_z;

            for (size_t j = 0; j < n; ++j) {  // тело B = j
                if (j == i) continue;

                // Относительные положение
                double dx = bodies[j].x - bodies[i].x;  // x_B - x_A
                double dy = bodies[j].y - bodies[i].y;
                double dz = bodies[j].z - bodies[i].z;
                double r = std::sqrt(dx * dx + dy * dy + dz * dz);
                double r2 = r * r;

                // n_BA = (x_B - x_A) / r (от A к B)
                double nx_BA = dx / r;
                double ny_BA = dy / r;
                double nz_BA = dz / r;

                // n_AB = -n_BA (от B к A)
                double nx_AB = -nx_BA;
                double ny_AB = -ny_BA;
                double nz_AB = -nz_BA;

                // Скорости v_A и v_B
                double vx_A = bodies[i].vx, vy_A = bodies[i].vy, vz_A = bodies[i].vz;
                double vx_B = bodies[j].vx, vy_B = bodies[j].vy, vz_B = bodies[j].vz;

                double vA2 = vx_A * vx_A + vy_A * vy_A + vz_A * vz_A;
                double vB2 = vx_B * vx_B + vy_B * vy_B + vz_B * vz_B;
                double dot_vA_vB = vx_A * vx_B + vy_A * vy_B + vz_A * vz_B;

                // (n_AB · v_B)^2
                double dot_nAB_vB = nx_AB * vx_B + ny_AB * vy_B + nz_AB * vz_B;
                double term_nab_vb2 = dot_nAB_vB * dot_nAB_vB;

                // sum_C != A G m_C / r_AC
                double sum_Gm_rAC = 0.0;
                for (size_t k = 0; k < n; ++k) {
                    if (k == i) continue;
                    double dxk = bodies[k].x - bodies[i].x;
                    double dyk = bodies[k].y - bodies[i].y;
                    double dzk = bodies[k].z - bodies[i].z;
                    double rk = std::sqrt(dxk * dxk + dyk * dyk + dzk * dzk);
                    sum_Gm_rAC += G * bodies[k].m / rk;
                }

                // sum_C != B G m_C / r_BC
                double sum_Gm_rBC = 0.0;
                for (size_t k = 0; k < n; ++k) {
                    if (k == j) continue;
                    double dxk = bodies[k].x - bodies[j].x;
                    double dyk = bodies[k].y - bodies[j].y;
                    double dzk = bodies[k].z - bodies[j].z;
                    double rk = std::sqrt(dxk * dxk + dyk * dyk + dzk * dzk);
                    sum_Gm_rBC += G * bodies[k].m / rk;
                }

                // 1/2 * ((x_B - x_A) · a_B)
                double dot_xBA_aB = dx * ax_newt[j] + dy * ay_newt[j] + dz * az_newt[j];
                double half_dot_xBA_aB = 0.5 * dot_xBA_aB;

                // Скобка для первого 1PN-термина
                double bracket = vA2 + 2.0 * vB2 - 4.0 * dot_vA_vB - 1.5 * term_nab_vb2
                    - 4.0 * sum_Gm_rAC - sum_Gm_rBC + half_dot_xBA_aB;

                // G m_B n_BA / r^2 * bracket
                double gm_n_r2 = G * bodies[j].m / r2;
                double contrib1_x = gm_n_r2 * bracket * nx_BA;
                double contrib1_y = gm_n_r2 * bracket * ny_BA;
                double contrib1_z = gm_n_r2 * bracket * nz_BA;

                // Второй термин: G m_B / r^2 * [n_AB · (4 v_A - 3 v_B)] * (v_A - v_B)
                double dot_nAB_4vA_3vB = nx_AB * (4.0 * vx_A - 3.0 * vx_B) +
                    ny_AB * (4.0 * vy_A - 3.0 * vy_B) +
                    nz_AB * (4.0 * vz_A - 3.0 * vz_B);
                double factor2 = gm_n_r2 * dot_nAB_4vA_3vB;

                double dvx = vx_A - vx_B;
                double dvy = vy_A - vy_B;
                double dvz = vz_A - vz_B;

                double contrib2_x = factor2 * dvx;
                double contrib2_y = factor2 * dvy;
                double contrib2_z = factor2 * dvz;

                // Добавляем вклады
                pn_ax += contrib1_x + contrib2_x;
                pn_ay += contrib1_y + contrib2_y;
                pn_az += contrib1_z + contrib2_z;
            }

            // Делим все PN-термины на c^2
            ax[i] += pn_ax / (c * c);
            ay[i] += pn_ay / (c * c);
            az[i] += pn_az / (c * c);
        }
    }
}

// Рунге-Кутта 6 порядка
void Algorithm::RungeKutta6(std::vector<Body>& bodies, double dt, bool usePNCorrection) {

    size_t n = bodies.size();
    std::vector<double> ax(n), ay(n), az(n);
    std::vector<Body> k1 = bodies, k2 = bodies, k3 = bodies, k4 = bodies, k5 = bodies, k6 = bodies;

    // Шаг 1
    Calculation_A(bodies, ax, ay, az, usePNCorrection);
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
    Calculation_A(k2, ax, ay, az, usePNCorrection);
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
    Calculation_A(k3, ax, ay, az, usePNCorrection);
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
    Calculation_A(k4, ax, ay, az, usePNCorrection);
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
    Calculation_A(k5, ax, ay, az, usePNCorrection);
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
    Calculation_A(k6, ax, ay, az, usePNCorrection);
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
void Algorithm::RungeKutta4(std::vector<Body>& bodies, double dt, bool usePNCorrection) {
    size_t n = bodies.size();
    std::vector<double> ax(n), ay(n), az(n);
    std::vector<Body> k1 = bodies, k2 = bodies, k3 = bodies, k4 = bodies;

    // Шаг 1
    Calculation_A(bodies, ax, ay, az, usePNCorrection);
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
    Calculation_A(k2, ax, ay, az, usePNCorrection);
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
    Calculation_A(k3, ax, ay, az, usePNCorrection);
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
    Calculation_A(k4, ax, ay, az, usePNCorrection);
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
void Algorithm::RungeKutta2(std::vector<Body>& bodies, double dt, bool usePNCorrection) {
    size_t n = bodies.size();
    std::vector<double> ax(n), ay(n), az(n);
    std::vector<Body> k1 = bodies, k2 = bodies;

    // Шаг 1
    Calculation_A(bodies, ax, ay, az, usePNCorrection);
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
    Calculation_A(k2, ax, ay, az, usePNCorrection);
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

void Algorithm::EulerMethod(std::vector<Body>& bodies, double dt, bool usePNCorrection) {
    size_t n = bodies.size();
    std::vector<double> ax(n), ay(n), az(n);

    Calculation_A(bodies, ax, ay, az, usePNCorrection);

    for (size_t i = 0; i < n; ++i) {

        bodies[i].vx += ax[i] * dt;
        bodies[i].vy += ay[i] * dt;
        bodies[i].vz += az[i] * dt;

        bodies[i].x += bodies[i].vx * dt;
        bodies[i].y += bodies[i].vy * dt;
        bodies[i].z += bodies[i].vz * dt;
    }
}

void Algorithm::Leapfrog(std::vector<Body>& bodies, double dt, bool usePNCorrection)
{
    size_t n = bodies.size();

    std::vector<double> ax(n), ay(n), az(n);
    Calculation_A(bodies, ax, ay, az, usePNCorrection);

    Calculation_A(bodies, ax, ay, az, usePNCorrection);

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

    Calculation_A(bodies, ax, ay, az, usePNCorrection);
    for (size_t i = 0; i < n; ++i) {
        bodies[i].vx += 0.5 * ax[i] * dt;
        bodies[i].vy += 0.5 * ay[i] * dt;
        bodies[i].vz += 0.5 * az[i] * dt;
    }
}