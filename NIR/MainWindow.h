#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#ifdef ERROR
#undef ERROR
#endif

#ifdef __cplusplus
extern "C" {
#endif

#include "TEOBResumS.h"
#undef I  // Исправляет ошибки в Qt

#ifdef __cplusplus
}
#endif
// ---------------------------

#include <QMainWindow>
#include "struct.h"
#include "qcustomplot.h"
#include <QPen>
#include <QColor>
#include <QVector>
#include <QEvent>
#include <QFileDialog>
#include <QTextStream>
#include <QRegularExpression>
#include <QtConcurrent>
#include <QFutureWatcher>
#include <QProgressBar>
#include "ui_MainWindow.h"
#include "visualisation.h"
#include <QFile>

extern "C" EOBParameters *EOBPars;

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow();
    
    bool eventFilter(QObject *obj, QEvent *event) override;

    void resizeEvent(QResizeEvent *event) override;                 // Для пересчета путей при ресайзе

private slots:

    bool highlightAndCheckFields(const QList<QLineEdit*>& fields);  // Функция для проверки LineEdit

    // Кнопки "Импорт", "Добавить", "Удалить" - Декартовы координаты
    void on_pushButton_Add_XYZ_clicked();
    void on_deleteButton_XYZ_clicked();
    void on_pushButtonImport_clicked();

    // Выбор учёта релятивистских эффектов
    void on_radioButton_noEffect_clicked();
    void on_radioButton_yesEffect_clicked();

    // Кнопки "Рассчитать", "Рассчитать и визуализировать", "Визуализировать"
    void on_pushButtonCalculate_clicked();
    void on_pushButtonCalculate_Visual_clicked();
    void on_pushButtonVisual_clicked();

    // Кнопки "Папка сохранения" и "Файл визуализации"
    void on_pushButton__Folder_clicked();
    void on_pushButton_FileVisual_clicked();
    
    // Кнопки "Импорт", "Добавить", "Удалить" - Кеплеровы элементы
    void on_pushButtonImport_Kepler_clicked();
    void on_deleteButton_Kepler_clicked();
    void on_pushButton_Add_Kepler_clicked();

    // Кнопка "Применить"
    void on_pushButtonApplySettings_clicked();
    // Кнопка "Синхронизировать с таблицей данных"
    //void on_pushButton_SyncBodies_clicked();
    // Кнопка "Авто-масштаб"
    void on_pushButton_PlotRescale_clicked();
    // Кнопка "Сохранить график (PNG)"
    void on_pushButton_PlotSavePNG_clicked();

    // Кнопка "Рассчитать TEOBresumS"
    void on_gw_btn_calculate_clicked();
    void on_dekart_cell_changed(int row, int col);
    void on_gw_combo_systemType_currentIndexChanged(int index);
    void on_tableWidget_PlotSettings_cellClicked(int row, int col);
    void syncPlotSettingsTable();

    // Обработка движения ползунка
    void on_gw_slider_tail_valueChanged(int value);

    // Кнопки "Сохранить график гравитационных волн" и "Сохранить график динамики слияния"
    void on_pushButton_save_GW_clicked();
    void on_pushButton_save_Dynamic_GW_clicked();

    // Функции 
    void on_gw_chk_showWave_stateChanged(int state);
    void on_gw_combo_waveColor_currentIndexChanged(const QString &color);
    void on_gw_spin_lineWidth_valueChanged(int value);
    void on_gw_combo_plotLayout_currentIndexChanged(int index);


    void on_btnBuildStudyPlot_clicked();
private:
    QVector<double> last_t;
    QVector<double> last_hp;
    QVector<double> last_hx;

    // Указатель на класс визуализации
    Visual *vis = nullptr;

    double x, y, z, vx, vy, vz;
    double step, time_interval;
    int curIndexRelativistic, curIndexMethods;                  // Индексы 
    int curIndexMethodsInitParams;
    std::string filename, pathFolder = "", pathFileVisual = "";
    std::vector<Body> body_system_XYZ;                          // Массив тел в "Декартовых координатах"
    std::vector<KeplerianElements> body_system_Keplerian;       // Массив тел в "Кеплеровых элементах орбиты"
    std::vector<std::string> labels_XYZ, labels_Kepler;
    QVector<double> last_x1, last_y1, last_x2, last_y2;

    QString filePath_InitParam_XYZ, filePath_InitParam_Kepler, qfilename;
    
    // Массивы для хранения имен и указателей на LineEdit для проверки при добавлении тел
    QList<QLineEdit*> requiredForAdd_XYZ;
    QList<QLineEdit*> requiredForAdd_Kepler;
    QList<QLineEdit*> requiredForStart;

    QProgressBar *progressBar;                 
    void visualizeData(bool keepScale);
    void startCalculation(bool autoVisualize);
    void updatePathLabels();

    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
