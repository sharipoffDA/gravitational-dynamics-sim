#include"system.h"
#include <numbers>
#include <QApplication>
#include "MainWindow.h"
#include <cstdlib>
#include <clocale>
#include <QLocale>
const double pi = std::numbers::pi;

int main(int argc, char* argv[]) {

        std::setlocale(LC_NUMERIC, "C");
        
        QLocale::setDefault(QLocale::c());
    try {
        std::cout.flush();
        QApplication app(argc, argv);
        std::cout.flush();

    MainWindow w;
    w.show();

    return app.exec();
    }
    catch (const std::exception& e) {
        std::cerr << "Ошибка: " << e.what() << std::endl;
        return 1;
    }
    catch (...) {
        std::cerr << "Неизвестная ошибка!" << std::endl;
        return 1;
    }
}