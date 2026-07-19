#include "system.h"
#include "MainWindow.h"
#include <QDoubleValidator>
#include <QDebug>
#include <QMouseEvent>
#include <fstream>
#include <cstdint>

// Объявление глобальных переменных TEOBResumS для правильной инициализации
extern "C" {
    extern EOBParameters *EOBPars; 
    extern NQCdata *NQC;             
}

const double pi = std::numbers::pi;

MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent), ui(new Ui::MainWindow) {
    ui->setupUi(this);

    ui->gw_plot_dynamics->setMinimumSize(50, 50);
    ui->gw_plot_waveform->setMinimumSize(50, 50);

    requiredForAdd_XYZ = {
        ui->lineEdit_ObjectName,
        ui->lineEdit_M,
        ui->lineEdit_X,
        ui->lineEdit_Y,
        ui->lineEdit_Z,
        ui->lineEdit_VX,
        ui->lineEdit_VY,
        ui->lineEdit_VZ
    };

    requiredForAdd_Kepler = {
        ui->lineEdit_ObjectName_Kepler,
        ui->lineEdit_a,
        ui->lineEdit_i,
        ui->lineEdit_omega,
        ui->lineEdit_Omega,
        ui->lineEdit_T_p_year,
        ui->lineEdit_M_Kepler,
        ui->lineEdit_ecc,
        ui->lineEdit_period_days
    };

    requiredForStart = {
        ui->lineEdit_SimulationTime,
        ui->lineEdit_TimeStep
    };

    auto addDoubleValidator = [](QLineEdit* le) {
        auto v = new QDoubleValidator(le);
        v->setNotation(QDoubleValidator::StandardNotation);
        le->setValidator(v);
    };

    ui->tableWidget_info_dekart->setColumnCount(8);
    ui->tableWidget_info_dekart->setHorizontalHeaderLabels({"Объект", "M", "X", "Y", "Z", "Vx", "Vy", "Vz"});

    ui->tableWidget_info_dekart->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
    ui->tableWidget_info_dekart->setEditTriggers(QAbstractItemView::DoubleClicked);     
    ui->tableWidget_info_dekart->setSelectionBehavior(QAbstractItemView::SelectRows);  
    ui->tableWidget_info_dekart->setAlternatingRowColors(true);                          
    ui->tableWidget_info_dekart->verticalHeader()->setVisible(false); 
    
    ui->tableWidget_info_Kepler->setColumnCount(9);
    ui->tableWidget_info_Kepler->setHorizontalHeaderLabels({"Объект", "a", "i", "ω", "Ω", "M", "e", "T", "T_p"});

    ui->tableWidget_info_Kepler->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
    ui->tableWidget_info_Kepler->setEditTriggers(QAbstractItemView::DoubleClicked);     
    ui->tableWidget_info_Kepler->setSelectionBehavior(QAbstractItemView::SelectRows);  
    ui->tableWidget_info_Kepler->setAlternatingRowColors(true);                          
    ui->tableWidget_info_Kepler->verticalHeader()->setVisible(false); 

    ui->comboBox_PN_correction->setEnabled(false);
    ui->customPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    ui->gw_plot_dynamics->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    ui->gw_plot_waveform->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);

    for (auto le : requiredForAdd_XYZ) le->installEventFilter(this);
    for (auto le : requiredForAdd_Kepler) le->installEventFilter(this);
    for (auto le : requiredForStart) le->installEventFilter(this);
    
    for (auto le : requiredForAdd_XYZ) {
        if (le != ui->lineEdit_ObjectName) addDoubleValidator(le);
    }
    for (auto le : requiredForAdd_Kepler) {
        if (le != ui->lineEdit_ObjectName_Kepler) addDoubleValidator(le);
    }
    for (auto le : requiredForStart) {
        addDoubleValidator(le);
    }

    ui->lineEdit_ObjectName->setMaxLength(20);
    ui->lineEdit_ObjectName_Kepler->setMaxLength(20);
    ui->lineEdit_PlotTitle->setMaxLength(20);

    progressBar = new QProgressBar(this);
    progressBar->setRange(0, 0); // Режим "бесконечной загрузки"
    progressBar->setVisible(false);
    progressBar->setMaximumWidth(200);
    ui->statusbar->addPermanentWidget(progressBar);

    connect(ui->gw_combo_systemType, QOverload<int>::of(&QComboBox::currentIndexChanged), this, &MainWindow::on_gw_combo_systemType_currentIndexChanged);
    connect(ui->tableWidget_info_dekart, &QTableWidget::cellChanged, this, &MainWindow::on_dekart_cell_changed);
    connect(ui->tabWidget, &QTabWidget::currentChanged, this, &MainWindow::syncPlotSettingsTable);
    connect(ui->gw_slider_tail, &QSlider::valueChanged, this, &MainWindow::on_gw_slider_tail_valueChanged);

    ui->tableWidget_info_dekart->setStyleSheet("QTableView::item:selected { background-color: #3399ff; } QLineEdit { background: white; color: black; }");
    ui->tableWidget_info_Kepler->setStyleSheet("QTableView::item:selected { background-color: #3399ff; } QLineEdit { background: white; color: black; }");
    ui->tableWidget_PlotSettings->setStyleSheet("QTableView::item:selected { background-color: #3399ff; } QLineEdit { background: white; color: black; }");

    ui->gw_plot_dynamics->setOpenGl(true);
    ui->gw_plot_waveform->setOpenGl(true);

    // Отключаем сглаживание для всего, кроме легенды (сильно экономит ресурсы)
    ui->gw_plot_dynamics->setNotAntialiasedElements(QCP::aeAll);
    ui->gw_plot_waveform->setNotAntialiasedElements(QCP::aeAll);
}

MainWindow::~MainWindow() {
    delete ui;
}

void MainWindow::on_radioButton_noEffect_clicked() {
    if (ui->radioButton_noEffect->isChecked()) {
        ui->comboBox_PN_correction->setEnabled(false);
    }   
}

void MainWindow::on_radioButton_yesEffect_clicked(){
    if (ui->radioButton_yesEffect->isChecked()) {
        ui->comboBox_PN_correction->setEnabled(true);
    }
}

// ДОБАВИТЬ В БУДУЩЕМ
// Для звёзд и компактных сверхмассивных объектов
void compute_initial_conditions(
    const KeplerianElements& data,
    double M_star,
    double epoch_year,
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
    double mass = data.mass * M_EARTH; // !!!!!!!!ПОКА В МАССАХ ЗЕМЛИ!!!!!!!!!!
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

// Функция для получения директории исполняемого файла
fs::path get_executable_dir() {
    #ifdef _WIN32
        wchar_t buffer[MAX_PATH];
        GetModuleFileNameW(NULL, buffer, MAX_PATH);
        return fs::path(buffer).parent_path();
    #else
        char buffer[PATH_MAX];
        ssize_t count = readlink("/proc/self/exe", buffer, PATH_MAX);
        if (count != -1) {
            buffer[count] = '\0';
            return fs::path(buffer).parent_path();
        }
        return fs::current_path();
    #endif
}

bool MainWindow::eventFilter(QObject *obj, QEvent *event)
{
    QLineEdit *lineEdit = qobject_cast<QLineEdit*>(obj);
    
    if (lineEdit) {
        if (event->type() == QEvent::FocusIn || event->type() == QEvent::MouseButtonPress) {
            if (!lineEdit->styleSheet().isEmpty()) {
                lineEdit->setStyleSheet("");
            }
        }
    }
    
    return QMainWindow::eventFilter(obj, event);
}


//<======================< Функция для проверки LineEdit >==========================>

bool MainWindow::highlightAndCheckFields(const QList<QLineEdit*>& fields)
{
    bool hasError = false;

    for (QLineEdit* le : fields) {
        QString text = le->text().trimmed();

        if (text.isEmpty()) {
            le->setStyleSheet("QLineEdit { border: 2px solid red; border-radius: 4px; padding: 1px; }");
            hasError = true;
        } else {
            le->setStyleSheet("");
        }
    }

    if (hasError) {
        QMessageBox::warning(this, "Ошибка ввода", 
                             "Не все обязательные поля заполнены!\n"
                             "Пожалуйста, введите значения в подсвеченные поля.");
    }

    return hasError;
}

//<======================< !Функция для проверки LineEdit >==========================>

//<======================< Кнопка "Добавить (Декартовы координаты)" >==========================>

void MainWindow::on_pushButton_Add_XYZ_clicked()
{
    if (highlightAndCheckFields(requiredForAdd_XYZ)) {
        return;
    }

    int row = ui->tableWidget_info_dekart->rowCount();
    ui->tableWidget_info_dekart->insertRow(row);

    // Добавление в таблицу
    for (int i = 0; i < requiredForAdd_XYZ.size(); ++i)
    {
        ui->tableWidget_info_dekart->setItem(row, i, new QTableWidgetItem(requiredForAdd_XYZ[i]->text()));
    }

    for (int col = 0; col < requiredForAdd_XYZ.size(); ++col) {
        ui->tableWidget_info_dekart->item(row, col)->setTextAlignment(Qt::AlignCenter);
    }

    labels_XYZ.push_back(ui->lineEdit_ObjectName->text().toStdString());

    body_system_XYZ.push_back(
        {ui->lineEdit_M->text().toDouble(), 
        ui->lineEdit_X->text().toDouble(), ui->lineEdit_Y->text().toDouble(), ui->lineEdit_Z->text().toDouble(),
        ui->lineEdit_VX->text().toDouble(), ui->lineEdit_VY->text().toDouble(), ui->lineEdit_VZ->text().toDouble()
        }
    );

    ui->tableWidget_info_dekart->scrollToItem(ui->tableWidget_info_dekart->item(row, 0));

    for (QLineEdit *le : requiredForAdd_XYZ) {
        le->blockSignals(true);
        le->clear();
        le->setStyleSheet("");
        le->blockSignals(false);
    }

    syncPlotSettingsTable();
}
//<======================< !Кнопка "Добавить (Декартовы координаты)" >==========================>

//<======================< Кнопка "Удалить (Декартовы координаты)" >==========================>

void MainWindow::on_deleteButton_XYZ_clicked()
{
    QItemSelectionModel *select = ui->tableWidget_info_dekart->selectionModel();
    QModelIndexList selectedRows = select->selectedRows();

    if (selectedRows.isEmpty()) return;

    std::sort(selectedRows.begin(), selectedRows.end(), [](const QModelIndex &a, const QModelIndex &b) {
        return a.row() > b.row();
    });

    QAbstractItemModel *model = ui->tableWidget_info_dekart->model();
    
    for (const QModelIndex &index : selectedRows) {
        int row = index.row();

        if (row < body_system_XYZ.size()) {
            labels_XYZ.erase(labels_XYZ.begin() + row);
            body_system_XYZ.erase(body_system_XYZ.begin() + row);
        }

        model->removeRow(row);
    }

    syncPlotSettingsTable();
}

//<======================< !Кнопка "Удалить (Декартовы координаты)" >==========================>

//<======================< Кнопка "Импорт (Декартовы координаты)" >==========================>
void MainWindow::on_pushButtonImport_clicked()
{
    filePath_InitParam_XYZ = QFileDialog::getOpenFileName(this, "Выберите файл для импорта", QDir::homePath(), "Данные (*.txt *.csv);;Все файлы (*.*)");
    if (filePath_InitParam_XYZ.isEmpty()) return; 

    QFile file(filePath_InitParam_XYZ);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        statusBar()->showMessage("Ошибка: не удалось открыть файл для импорта!");
        return;
    }

    QTextStream in(&file);
    QList<QStringList> validLines;
    int lineNum = 0;

    // 1. Предварительная проверка файла
    while (!in.atEnd()) {
        QString line = in.readLine().trimmed();
        lineNum++;
        if (line.isEmpty()) continue; 

        QStringList parts = line.split(QRegularExpression("\\s+"), Qt::SkipEmptyParts);
        
        if (parts.size() != 8) {
            QMessageBox::warning(this, "Ошибка формата", QString("Неверный формат в строке %1. Ожидается: Имя и 7 чисел (M, X, Y, Z, Vx, Vy, Vz).").arg(lineNum));
            return; // Прерываем загрузку
        }

        bool ok = true;
        for (int i = 1; i <= 7; ++i) {
            parts[i].toDouble(&ok);
            if (!ok) {
                QMessageBox::warning(this, "Ошибка данных", QString("Ошибка чтения числа в строке %1, колонка %2").arg(lineNum).arg(i + 1));
                return; // Прерываем загрузку
            }
        }
        validLines.append(parts);
    }
    file.close();

    // 2. Очищаем старые данные и загружаем проверенные
    ui->tableWidget_info_dekart->setRowCount(0);
    body_system_XYZ.clear();
    labels_XYZ.clear();

    for (const QStringList& parts : validLines) {
        labels_XYZ.push_back(parts[0].toStdString());
        body_system_XYZ.push_back({parts[1].toDouble(), parts[2].toDouble(), parts[3].toDouble(), parts[4].toDouble(), parts[5].toDouble(), parts[6].toDouble(), parts[7].toDouble()});

        int row = ui->tableWidget_info_dekart->rowCount();
        ui->tableWidget_info_dekart->insertRow(row);
        for (int col = 0; col < 8; ++col) {
            QTableWidgetItem *item = new QTableWidgetItem(parts[col]);
            item->setTextAlignment(Qt::AlignCenter);
            ui->tableWidget_info_dekart->setItem(row, col, item);
        }
    }
    statusBar()->showMessage(QString("Успешно загружен файл: %1").arg(filePath_InitParam_XYZ));
}

//<======================< Кнопка "Импорт (Кеплеровы элементы)" >==========================>
void MainWindow::on_pushButtonImport_Kepler_clicked()
{
    filePath_InitParam_Kepler = QFileDialog::getOpenFileName(this, "Выберите файл для импорта", QDir::homePath(), "Данные (*.txt *.csv);;Все файлы (*.*)");
    if (filePath_InitParam_Kepler.isEmpty()) return; 

    QFile file(filePath_InitParam_Kepler);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        statusBar()->showMessage("Ошибка: не удалось открыть файл для импорта!");
        return;
    }

    QTextStream in(&file);
    QList<QStringList> validLines;
    int lineNum = 0;

    // 1. Предварительная проверка
    while (!in.atEnd()) {
        QString line = in.readLine().trimmed();
        lineNum++;
        if (line.isEmpty()) continue; 

        QStringList parts = line.split(QRegularExpression("\\s+"), Qt::SkipEmptyParts);
        
        if (parts.size() != 9) {
            QMessageBox::warning(this, "Ошибка формата", QString("Неверный формат в строке %1. Ожидается: Имя и 8 чисел (a, i, ω, Ω, M, e, T, T_p).").arg(lineNum));
            return;
        }

        bool ok = true;
        for (int i = 1; i <= 8; ++i) {
            parts[i].toDouble(&ok);
            if (!ok) {
                QMessageBox::warning(this, "Ошибка данных", QString("Ошибка чтения числа в строке %1, колонка %2").arg(lineNum).arg(i + 1));
                return;
            }
        }
        validLines.append(parts);
    }
    file.close();

    // 2. Загрузка
    ui->tableWidget_info_Kepler->setRowCount(0);
    body_system_Keplerian.clear();
    labels_Kepler.clear();

    for (const QStringList& parts : validLines) {
        labels_Kepler.push_back(parts[0].toStdString());
        body_system_Keplerian.push_back({parts[0].toStdString(), parts[1].toDouble(), parts[2].toDouble(), parts[3].toDouble(), parts[4].toDouble(), parts[5].toDouble(), parts[6].toDouble(), parts[7].toDouble(), parts[8].toDouble()});

        int row = ui->tableWidget_info_Kepler->rowCount();
        ui->tableWidget_info_Kepler->insertRow(row);
        for (int col = 0; col < 9; ++col) {
            QTableWidgetItem *item = new QTableWidgetItem(parts[col]);
            item->setTextAlignment(Qt::AlignCenter);
            ui->tableWidget_info_Kepler->setItem(row, col, item);
        }
    }
    statusBar()->showMessage(QString("Успешно загружен файл: %1").arg(filePath_InitParam_Kepler));
}

//<======================< Кнопка "Добавить (Кеплеровы элементы)" >==========================>

void MainWindow::on_pushButton_Add_Kepler_clicked()
{
    if (highlightAndCheckFields(requiredForAdd_Kepler)) {
        return;
    }

    int row = ui->tableWidget_info_Kepler->rowCount();
    ui->tableWidget_info_Kepler->insertRow(row);

    // Добавление в таблицу
    for (int i = 0; i < requiredForAdd_Kepler.size(); ++i)
    {
        ui->tableWidget_info_Kepler->setItem(row, i, new QTableWidgetItem(requiredForAdd_Kepler[i]->text()));
    }

    for (int col = 0; col < requiredForAdd_Kepler.size(); ++col) {
        ui->tableWidget_info_Kepler->item(row, col)->setTextAlignment(Qt::AlignCenter);
    }

    labels_Kepler.push_back(ui->lineEdit_ObjectName->text().toStdString());

    body_system_Keplerian.push_back(
        {ui->lineEdit_ObjectName_Kepler->text().toStdString(), ui->lineEdit_a->text().toDouble(), 
        ui->lineEdit_i->text().toDouble(), ui->lineEdit_omega->text().toDouble(), ui->lineEdit_Omega->text().toDouble(),
        ui->lineEdit_M_Kepler->text().toDouble(), ui->lineEdit_ecc->text().toDouble(), ui->lineEdit_period_days->text().toDouble(),
        ui->lineEdit_T_p_year->text().toDouble()
        }
    );

    ui->tableWidget_info_Kepler->scrollToItem(ui->tableWidget_info_Kepler->item(row, 0));

    for (QLineEdit *le : requiredForAdd_Kepler) {
        le->blockSignals(true);
        le->clear();
        le->setStyleSheet("");
        le->blockSignals(false);
    }

    syncPlotSettingsTable();
}
//<======================< !Кнопка "Добавить (Кеплеровы элементы)" >==========================>

//<======================< Кнопка "Удалить (Кеплеровы элементы)" >==========================>

void MainWindow::on_deleteButton_Kepler_clicked()
{
    QItemSelectionModel *select = ui->tableWidget_info_Kepler->selectionModel();
    QModelIndexList selectedRows = select->selectedRows();

    if (selectedRows.isEmpty()) return;

    std::sort(selectedRows.begin(), selectedRows.end(), [](const QModelIndex &a, const QModelIndex &b) {
        return a.row() > b.row();
    });

    // 3. Удаляем из вектора и из модели
    QAbstractItemModel *model = ui->tableWidget_info_Kepler->model();
    
    for (const QModelIndex &index : selectedRows) {
        int row = index.row();

        if (row < body_system_Keplerian.size()) {
            labels_Kepler.erase(labels_Kepler.begin() + row);
            body_system_Keplerian.erase(body_system_Keplerian.begin() + row);
        }

        model->removeRow(row);
    }

    syncPlotSettingsTable();
}

//<======================< !Кнопка "Удалить (Кеплеровы элементы)" >==========================>


//<======================< Функция для изменения имени файла, если есть с таким же именем >==========================>

QString getUniqueFilePath(QString fullPath) {
    QFileInfo fileInfo(fullPath);
    
    // Если файла нет, возвращаем исходный путь
    if (!fileInfo.exists()) {
        return fullPath;
    }

    // Разделяем путь на составляющие: папка, имя без расширения и само расширение
    QString path = fileInfo.absolutePath(); // Напр: C:/Results
    QString baseName = fileInfo.baseName(); // Напр: simulation_results
    QString extension = fileInfo.completeSuffix(); // Напр: txt
    
    int counter = 1;
    QString newPath = fullPath;

    while (QFile::exists(newPath)) {
        // Формируем новое имя: C:/Results/simulation_results(1).txt
        newPath = QString("%1/%2 (%3).%4")
                    .arg(path)
                    .arg(baseName)
                    .arg(counter)
                    .arg(extension);
        counter++;
    }

    return newPath;
}

//<======================< !Функция для изменения имени файла, если есть с таким же именем >==========================>

//<======================< Кнопка "Рассчитать" >==========================>

void MainWindow::on_pushButtonCalculate_clicked()
{
    startCalculation(false); 
}

//<======================< !Кнопка "Рассчитать" >==========================>

//<======================< Кнопка "Рассчитать и визуализировать" >==========================>

void MainWindow::on_pushButtonCalculate_Visual_clicked()
{
    startCalculation(true); 
}

//<======================< !Кнопка "Рассчитать и визуализировать" >==========================>

//<======================< Кнопка "Выбрать папку" >==========================>
void MainWindow::on_pushButton__Folder_clicked()
{
    QString dir = QFileDialog::getExistingDirectory(
        this, 
        "Выберите папку для сохранения результатов",
        QDir::homePath(), 
        QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks
    );

    if (!dir.isEmpty()) {
        pathFolder = dir.toStdString() + "/";
        updatePathLabels();
    }else return;

    QFontMetrics metrics(ui->lbl_path_Folder->font());
    QString elidedText = metrics.elidedText(dir, Qt::ElideMiddle, ui->lbl_path_Folder->width());
    ui->lbl_path_Folder->setText(elidedText);
}
//<======================< !Кнопка "Выбрать папку" >==========================>


//<======================< Кнопка "Выбрать файл" >==========================>
void MainWindow::on_pushButton_FileVisual_clicked() {
    QString filePath_Visual = QFileDialog::getOpenFileName(this, "Выберите файл", QDir::homePath(), "Данные (*.txt *.csv)");

    if (filePath_Visual.isEmpty()) return;

    if (!filePath_Visual.isEmpty()) {
        pathFileVisual = filePath_Visual.toLocal8Bit().toStdString(); 
        updatePathLabels();

        QFontMetrics metrics(ui->lbl_path_FileVisual->font());
        QString elidedPath = metrics.elidedText(filePath_Visual, Qt::ElideMiddle, ui->lbl_path_FileVisual->width());
        
        ui->lbl_path_FileVisual->setText(elidedPath);
        ui->lbl_path_FileVisual->setToolTip(filePath_Visual); 
    }

    // Читаем первую строку файла (имена объектов)
    QFile file(filePath_Visual);
    if (file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        QTextStream in(&file);
        QString firstLine = in.readLine();
        file.close();

        QStringList objectNames = firstLine.split(" ", Qt::SkipEmptyParts);

        // Если первая строка содержит цифры (старый формат), не заполняем таблицу имен
        bool isNumeric = false;
        if (!objectNames.isEmpty()) objectNames[0].toDouble(&isNumeric);

        if (!isNumeric && !objectNames.isEmpty()) {
            ui->tableWidget_PlotSettings->setRowCount(0);
            
            std::vector<QColor> defaultColors = {Qt::yellow, Qt::gray, Qt::red, Qt::blue, Qt::magenta, Qt::cyan};

            for (int i = 0; i < objectNames.size(); ++i) {
                int row = ui->tableWidget_PlotSettings->rowCount();
                ui->tableWidget_PlotSettings->insertRow(row);

                // Чекбокс видимости
                QTableWidgetItem* checkItem = new QTableWidgetItem();
                checkItem->setFlags(Qt::ItemIsUserCheckable | Qt::ItemIsEnabled);
                checkItem->setCheckState(Qt::Checked);
                ui->tableWidget_PlotSettings->setItem(row, 0, checkItem);

                // Имя из файла
                ui->tableWidget_PlotSettings->setItem(row, 1, new QTableWidgetItem(objectNames[i]));

                // Цвет
                QTableWidgetItem* colorItem = new QTableWidgetItem();
                QColor col = (i < defaultColors.size()) ? defaultColors[i] : QColor::fromHsv((i * 50) % 360, 200, 200);
                colorItem->setBackground(col);
                ui->tableWidget_PlotSettings->setItem(row, 2, colorItem);
            }
        }
    }
}
//<======================< !Кнопка "Выбрать файл" >==========================>


//<======================< Кнопка "Визуализировать" >==========================>
void MainWindow::on_pushButtonVisual_clicked()
{
    visualizeData(false);
}
//<======================< !Кнопка "Визуализировать" >==========================>

//<======================< Кнопка "Применить" >==========================>

void MainWindow::on_pushButtonApplySettings_clicked()
{
    visualizeData(true);
}

//<======================< !Кнопка "Применить" >==========================>

void MainWindow::visualizeData(bool keepScale)
{
    if (pathFileVisual.empty()) {
        QMessageBox::warning(this, "Внимание", "Пожалуйста, сначала выберите файл для визуализации!");
        return;
    }

    QString qPath = QString::fromLocal8Bit(pathFileVisual.c_str());
    QFile file(qPath); 
    
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        QMessageBox::critical(this, "Ошибка", "Не удалось открыть файл для визуализации!");
        return;
    }

    PlotSettings settings;
    settings.title = ui->lineEdit_PlotTitle->text();
    settings.projection = ui->comboBox_PlotProjection->currentIndex();
    settings.theme = ui->comboBox_PlotTheme->currentIndex();
    settings.fontSize = ui->spinBox_PlotFontSize->value();
    
    settings.showAxisLabels = ui->checkBox_PlotAxisLabels->isChecked(); 
    settings.showLegend = ui->checkBox_PlotLegend->isChecked();
    settings.showGrid = ui->checkBox_PlotGrid->isChecked();
    
    settings.keepScale = keepScale; 

    for (int i = 0; i < ui->tableWidget_PlotSettings->rowCount(); ++i) {
        PlotSettings::BodyStyle bs;
        QTableWidgetItem* checkItem = ui->tableWidget_PlotSettings->item(i, 0);
        bs.visible = checkItem ? (checkItem->checkState() == Qt::Checked) : true;

        QTableWidgetItem* nameItem = ui->tableWidget_PlotSettings->item(i, 1);
        bs.name = nameItem ? nameItem->text() : QString("Тело %1").arg(i + 1);

        QTableWidgetItem* colorItem = ui->tableWidget_PlotSettings->item(i, 2);
        bs.color = colorItem ? colorItem->background().color() : Qt::white;

        settings.bodyStyles.push_back(bs);
    }

    std::vector<Body> dummy_bodies;
    std::vector<double> dummy_t, dummy_e, dummy_p;
    Visual vis(dummy_bodies, dummy_t, dummy_e, dummy_p);
    
    vis.PlotOrbits(ui->customPlot, pathFileVisual, settings);
}

//<======================< Кнопка "Синхронизировать с таблицей данных" >======================>

// void MainWindow::on_pushButton_SyncBodies_clicked()
// {
//     ui->tableWidget_PlotSettings->setRowCount(0); // Очищаем таблицу

//     // Собираем имена тел из списков (проверяем Декартовы и Кеплеровы)
//     std::vector<std::string> allLabels = labels_XYZ;
//     if (allLabels.empty()) {
//         allLabels = labels_Kepler;
//     }

//     if (allLabels.empty()) {
//         QMessageBox::information(this, "Инфо", "Нет добавленных тел для синхронизации.");
//         return;
//     }

//     // Базовые цвета для красивой отрисовки
//     std::vector<QColor> defaultColors = {
//         Qt::yellow, Qt::gray, Qt::red, Qt::blue, 
//         QColor(255, 140, 0), QColor(210, 180, 140), 
//         Qt::magenta, Qt::cyan, QColor(0, 0, 139)
//     };

//     for (size_t i = 0; i < allLabels.size(); ++i) {
//         int row = ui->tableWidget_PlotSettings->rowCount();
//         ui->tableWidget_PlotSettings->insertRow(row);

//         // Колонка 0: Видимость (Checkbox)
//         QTableWidgetItem* checkItem = new QTableWidgetItem();
//         checkItem->setFlags(Qt::ItemIsUserCheckable | Qt::ItemIsEnabled);
//         checkItem->setCheckState(Qt::Checked);
//         ui->tableWidget_PlotSettings->setItem(row, 0, checkItem);

//         // Колонка 1: Имя в легенде
//         QTableWidgetItem* nameItem = new QTableWidgetItem(QString::fromStdString(allLabels[i]));
//         ui->tableWidget_PlotSettings->setItem(row, 1, nameItem);

//         // Колонка 2: Цвет орбиты
//         QTableWidgetItem* colorItem = new QTableWidgetItem();
//         QColor col = (i < defaultColors.size()) ? defaultColors[i] : QColor::fromHsv((i * 50) % 360, 200, 200);
//         colorItem->setBackground(col);
//         ui->tableWidget_PlotSettings->setItem(row, 2, colorItem);
//     }
// }

//<======================< !Кнопка "Синхронизировать..." >======================>


//<======================< Кнопка "Авто-масштаб" >======================>

void MainWindow::on_pushButton_PlotRescale_clicked()
{
    ui->customPlot->rescaleAxes();
    ui->customPlot->replot();
}

//<======================< !Кнопка "Авто-масштаб" >======================>


//<======================< Кнопка "Сохранить график (PNG)" >======================>

void MainWindow::on_pushButton_PlotSavePNG_clicked()
{
    // Открываем диалог для выбора места сохранения
    QString fileName = QFileDialog::getSaveFileName(this, 
                                                    "Сохранить график", 
                                                    QDir::homePath() + "/plot.png", 
                                                    "Images (*.png)");
    if (!fileName.isEmpty()) {
        // Проверяем, добавил ли пользователь расширение .png
        if (!fileName.endsWith(".png", Qt::CaseInsensitive)) {
            fileName += ".png";
        }
        
        // Сохраняем в высоком разрешении (Full HD)
        bool success = ui->customPlot->savePng(fileName, 1080, 1080, 2.0, 100);
        
        if (success) {
            statusBar()->showMessage("График успешно сохранен: " + fileName, 5000);
        } else {
            QMessageBox::warning(this, "Ошибка", "Не удалось сохранить график.");
        }
    }
}

//<======================< !Кнопка "Сохранить график (PNG)" >======================>

// Общая логика
void MainWindow::startCalculation(bool autoVisualize)
{
    if (highlightAndCheckFields(requiredForStart)) return;

    double t_interval = ui->lineEdit_SimulationTime->text().toDouble() * year_to_sec;
    double t_step = ui->lineEdit_TimeStep->text().toDouble() * day_to_sec;
    int methodIdx = ui->comboBox_Methods->currentIndex();
    int relIdx = ui->comboBox_PN_correction->currentIndex();
    bool usePN = ui->radioButton_yesEffect->isChecked();
    int tabIdx = ui->tabWidget->currentIndex();

    if (tabIdx == 0 && ui->tableWidget_info_dekart->rowCount() == 0) {
        QMessageBox::warning(this, "Ошибка ввода", "Нет начальных параметров моделирования!");
        return;
    }

    QString baseName = (usePN && relIdx == 0) ? 
        QString("BodySystem_RK%1_1PN.txt").arg((methodIdx == 0) ? 2 : (methodIdx == 1) ? 4 : 6) : 
        QString("BodySystem_RK%1.txt").arg((methodIdx == 0) ? 2 : (methodIdx == 1) ? 4 : 6);
    
    // Получаем полный путь в QString
    QString fullPath = getUniqueFilePath(QString::fromStdString(pathFolder) + baseName);
    
    filename = fullPath.toLocal8Bit().toStdString();
    
    ui->pushButtonCalculate->setEnabled(false);
    ui->pushButtonCalculate_Visual->setEnabled(false);
    if (progressBar) progressBar->setVisible(true);
    statusBar()->showMessage("Идет расчет системы...");

    std::vector<Body> initial_bodies = body_system_XYZ;
    std::vector<std::string> labels = labels_XYZ;
    std::string out_file = filename; // Здесь уже лежит Local8Bit

    QFutureWatcher<void> *watcher = new QFutureWatcher<void>(this);
    
    connect(watcher, &QFutureWatcher<void>::finished, this, [this, watcher, out_file, autoVisualize]() {
        ui->pushButtonCalculate->setEnabled(true);
        ui->pushButtonCalculate_Visual->setEnabled(true);
        if (progressBar) progressBar->setVisible(false);
        
        if (QFile::exists(QString::fromLocal8Bit(out_file.c_str()))) {
            statusBar()->showMessage("Расчет окончен! Файл успешно создан!", 5000);
            
            if (autoVisualize) {
                pathFileVisual = out_file;
                ui->lbl_path_FileVisual->setText(QString::fromLocal8Bit(out_file.c_str()));
                
                visualizeData(false);
                ui->toolBox->setCurrentIndex(1);
            }
        } else {
            statusBar()->showMessage("Ошибка: Файл не был создан!", 5000);
        }
        watcher->deleteLater();
    });

    QFuture<void> future = QtConcurrent::run([initial_bodies, out_file, t_interval, t_step, usePN, methodIdx, labels]() mutable {
        BodySystem system(initial_bodies);
        std::vector<std::string> labels_copy = labels; 
        
        if (usePN) {
            if (methodIdx == 0) system.SaveToFile_RK2(out_file, t_interval, t_step, labels_copy, true);
            else if (methodIdx == 1) system.SaveToFile_RK4(out_file, t_interval, t_step, labels_copy, true);
            else if (methodIdx == 2) system.SaveToFile_RK6(out_file, t_interval, t_step, labels_copy, true);
        } else {
            if (methodIdx == 0) system.SaveToFile_RK2(out_file, t_interval, t_step, labels_copy, false);
            else if (methodIdx == 1) system.SaveToFile_RK4(out_file, t_interval, t_step, labels_copy, false);
            else if (methodIdx == 2) system.SaveToFile_RK6(out_file, t_interval, t_step, labels_copy, false);
        }
    });

    watcher->setFuture(future);
}


//<======================<        ВКЛАДКА         >======================>
//<======================< "Гравитационные волны" >======================>

void MainWindow::on_gw_btn_calculate_clicked() 
{
    gsl_set_error_handler_off();
    
    // Сбор и валидация данных
    double m1 = ui->lineEdit_m1->text().toDouble(); 
    double m2 = ui->lineEdit_m2->text().toDouble();
    double f0_hz = ui->lineEdit_f0->text().toDouble(); 
    int sys_type = ui->gw_combo_systemType->currentIndex();

    if (m1 <= 0 || m2 <= 0 || f0_hz <= 0) {
        statusBar()->showMessage("Ошибка: массы и начальная частота должны быть положительными!");
        return;
    }

    // Базовая настройка
    EOBParameters *pars = nullptr;
    EOBParameters_alloc(&pars);
    
    // Устанавливаем тип системы ДО вызова дефолтов
    EOBParameters_defaults(sys_type, pars);
    pars->binary = sys_type; // 0: BBH, 1: BNS, 2: BHNS

    // Установка масс
    pars->M = m1 + m2;
    pars->q = (m1 >= m2) ? (m1 / m2) : (m2 / m1); // Гарантируем q >= 1
    
    // Перевод частоты в геометрические единицы
    double M_sec = pars->M * MSUN_S; 
    pars->initial_frequency = f0_hz * M_sec;

    // Спины
    pars->chi1z = ui->lineEdit_s1z->text().toDouble();
    pars->chi2z = ui->lineEdit_s2z->text().toDouble();
    pars->chi1 = pars->chi1z; 
    pars->chi2 = pars->chi2z;
    
    // Приливные параметры 
    pars->LambdaAl2 = ui->lineEdit_LambdaAl2->text().toDouble();
    pars->LambdaBl2 = ui->lineEdit_LambdaBl2->text().toDouble();
    
    // Настройка физики приливов в зависимости от типа
    if (sys_type == 1) { // BNS
        pars->use_tidal = TIDES_TEOBRESUM; 
        pars->use_tidal_gravitomagnetic = 1;
        pars->use_lambdaell_fits = 1;
    } else if (sys_type == 2) { // BHNS
        pars->use_tidal = TIDES_TEOBRESUM_BHNS; 
        pars->use_tidal_gravitomagnetic = 1;
        pars->use_lambdaell_fits = 1;
    } else { // BBH
        pars->use_tidal = 0;
        pars->use_tidal_gravitomagnetic = 0;
    }

    // Общие настройки
    pars->distance = ui->lineEdit_dist->text().toDouble();
    pars->use_geometric_units = 1;  
    pars->domain = 0; // Time domain
    pars->dt = 0.5;   // Адаптивный шаг по времени
    
    // Настройки ODE
    pars->ode_stop_radius = 0.0;
    pars->ode_stop_afterNdt = 100;
    pars->ode_stop_after_peak = 1; 
    pars->ode_tmax = 2e6;
    pars->ode_abstol = 1e-13;
    pars->ode_reltol = 1e-11;
    pars->srate_interp = 16384;

    // Расчет внутренних коэффициентов (nu, kappa, X1, X2)
    EOBParameters *old_pars = EOBPars;
    EOBPars = pars; 
    
    // Вызов с нулями (default_choice=0, firstcall=0) согласно сигнатуре в TEOBResumSPars.c
    if (eob_set_params(0, 0) != OK) {
        statusBar()->showMessage("Ошибка: eob_set_params не смог рассчитать коэффициенты.");
        EOBPars = old_pars;
        EOBParameters_free(pars);
        return;
    }

    // Отключаем вывод в файлы, так как забираем данные напрямую из оперативной памяти
    pars->output_dynamics = 0;
    pars->output_hpc = 0;
    pars->output_multipoles = 0;
    pars->output_nqc = 0;

    ui->gw_btn_calculate->setEnabled(false);
    if (progressBar) progressBar->setVisible(true);
    statusBar()->showMessage("Запуск моделирования TEOBResumS...");

    // Получаем имена объектов ПЕРЕД запуском потока
    QString name1 = ui->gw_lineEdit_name1->text().isEmpty() ? "Объект 1" : ui->gw_lineEdit_name1->text();
    QString name2 = ui->gw_lineEdit_name2->text().isEmpty() ? "Объект 2" : ui->gw_lineEdit_name2->text();

    // Подготовка контейнеров для результатов
    struct EOBResults {
        Waveform *hpc = nullptr;
        WaveformFD *hfpc = nullptr;
        Waveform_lm *hmodes = nullptr;
        WaveformFD_lm *hfmodes = nullptr;
        Dynamics *dynf = nullptr;
        Waveform_lm *hTmodes = nullptr;
        Waveform_lm *hTmmodes = nullptr;
        Waveform_lm *hT0modes = nullptr;
        WaveformFD_lm *hfTmodes = nullptr;

        void freeAll() {
            if (hpc) Waveform_free(hpc);
            if (hfpc) WaveformFD_free(hfpc);
            if (hmodes) Waveform_lm_free(hmodes);
            if (hfmodes) WaveformFD_lm_free(hfmodes);
            if (dynf) Dynamics_free(dynf);
            if (hTmodes) Waveform_lm_free(hTmodes);
            if (hTmmodes) Waveform_lm_free(hTmmodes);
            if (hT0modes) Waveform_lm_free(hT0modes);
            if (hfTmodes) WaveformFD_lm_free(hfTmodes);
        }
    };

    auto *results = new EOBResults();
    QFutureWatcher<int> *watcher = new QFutureWatcher<int>(this);
    
    connect(watcher, &QFutureWatcher<int>::finished, this, [=, this]() {
        ui->gw_btn_calculate->setEnabled(true);
        if (progressBar) progressBar->setVisible(false);
        
        int status = watcher->result();

        if (status == OK && results->hpc && results->hpc->size > 1) {
            last_t.clear(); last_hp.clear(); last_hx.clear();
            last_t.reserve(results->hpc->size);
            last_hp.reserve(results->hpc->size);
            last_hx.reserve(results->hpc->size);

            for (int i = 0; i < results->hpc->size; ++i) {
                last_t.append(results->hpc->time[i]);
                last_hp.append(results->hpc->real[i]); 
                last_hx.append(results->hpc->imag[i]); 
            }

            // Извлечение траекторий из dynf
            if (results->dynf && results->dynf->size > 0) {
                last_x1.clear(); last_y1.clear(); last_x2.clear(); last_y2.clear();
                double M_tot = m1 + m2;
                for (int i = 0; i < results->dynf->size; ++i) {
                    double r = results->dynf->data[EOB_RAD][i];
                    double phi = results->dynf->data[EOB_PHI][i];
                    double x_rel = r * std::cos(phi);
                    double y_rel = r * std::sin(phi);
                    last_x1.append((m2 / M_tot) * x_rel);
                    last_y1.append((m2 / M_tot) * y_rel);
                    last_x2.append(-(m1 / M_tot) * x_rel);
                    last_y2.append(-(m1 / M_tot) * y_rel);
                }
            }

            if (!vis) {
                std::vector<double> empty;
                vis = new Visual(body_system_XYZ, empty, empty, empty);
            }
            
            int currentFontSize = ui->gw_spin_fontSize->value();
            int tailValue = ui->gw_slider_tail->value();


            vis->PlotGWWaveform(ui->gw_plot_waveform, last_t, last_hp, last_hx, currentFontSize);
            vis->PlotGWDynamics(ui->gw_plot_dynamics, last_x1, last_y1, last_x2, last_y2, tailValue, name1, name2, currentFontSize);
            
            statusBar()->showMessage(QString("Расчет завершен успешно. Точек: %1").arg(results->hpc->size));
        } else {
            statusBar()->showMessage(QString("Ошибка: симуляция прервана (код %1). Проверьте параметры.").arg(status));
        }

        // Очистка памяти
        results->freeAll(); 
        delete results;     
        
        EOBParameters_free(pars);
        if (EOBPars == pars) EOBPars = old_pars;
        
        watcher->deleteLater();
    });

    QFuture<int> future = QtConcurrent::run([results]() {
        extern NQCdata *NQC;
        NQC = nullptr; 
        
        return EOBRun(
            &(results->hpc), &(results->hfpc), &(results->hmodes), 
            &(results->hfmodes), &(results->dynf), &(results->hTmodes), 
            &(results->hTmmodes), &(results->hT0modes), &(results->hfTmodes), 
            0, 1
        );
    });

    watcher->setFuture(future);
} 

//<======================< Кнопка "Сохранить график гравитационных волн" >======================>

// void MainWindow::on_pushButton_save_GW_clicked()
// {
//     QString fileName = QFileDialog::getSaveFileName(this, 
//                                                     "Сохранить график", 
//                                                     QDir::homePath() + "/waveform.pdf", 
//                                                     "PDF файлы (*.pdf)");
//     if (!fileName.isEmpty()) {
//         // Проверяем, добавил ли пользователь расширение .png
//         if (!fileName.endsWith(".pdf", Qt::CaseInsensitive)) {
//             fileName += ".pdf";
//         }
        
//         // Сохраняем в высоком разрешении (Full HD) с увеличенным скейлингом для качества
//         bool success = ui->gw_plot_waveform->savePdf(fileName, ui->gw_plot_waveform->width(), ui->gw_plot_waveform->height());
        
//         if (success) {
//             statusBar()->showMessage("График успешно сохранен: " + fileName, 5000);
//         } else {
//             QMessageBox::warning(this, "Ошибка", "Не удалось сохранить график.");
//         }
//     }
// }

void MainWindow::on_pushButton_save_GW_clicked()
{
    QString fileName = QFileDialog::getSaveFileName(this, 
                                                    "Сохранить график", 
                                                    QDir::homePath() + "/waveform.png", 
                                                    "Image (*.png)");
    if (!fileName.isEmpty()) {
        // Проверяем, добавил ли пользователь расширение .png
        if (!fileName.endsWith(".png", Qt::CaseInsensitive)) {
            fileName += ".png";
        }
        
        // Сохраняем в высоком разрешении (Full HD) с увеличенным скейлингом для качества
        bool success = ui->gw_plot_waveform->savePng(fileName, 1080, 1080, 2.0, 100);
        
        if (success) {
            statusBar()->showMessage("График успешно сохранен: " + fileName, 5000);
        } else {
            QMessageBox::warning(this, "Ошибка", "Не удалось сохранить график.");
        }
    }
}

//<======================< !Кнопка "Сохранить график гравитационных волн" >======================>

//<======================< Кнопка "Сохранить график динамики слияния" >======================>

void MainWindow::on_pushButton_save_Dynamic_GW_clicked()
{
    QString fileName = QFileDialog::getSaveFileName(this, 
                                                    "Сохранить график", 
                                                    QDir::homePath() + "/dynamic_merge.png", 
                                                    "Images (*.png)");
    if (!fileName.isEmpty()) {
        // Проверяем, добавил ли пользователь расширение .png
        if (!fileName.endsWith(".png", Qt::CaseInsensitive)) {
            fileName += ".png";
        }
        
        // Сохраняем в высоком разрешении (Full HD)
        bool success = ui->gw_plot_dynamics->savePng(fileName, 1080, 1080, 2.0, 100);
        
        if (success) {
            statusBar()->showMessage("График успешно сохранен: " + fileName, 5000);
        } else {
            QMessageBox::warning(this, "Ошибка", "Не удалось сохранить график.");
        }
    }
}

//<======================< !Кнопка "Сохранить график динамики слияния" >======================>

void MainWindow::on_gw_combo_systemType_currentIndexChanged(int index)
{
    bool isBNS = (index == 1);
    bool isBHNS = (index == 2);

    ui->lineEdit_LambdaAl2->setEnabled(isBNS); 
    ui->lineEdit_LambdaBl2->setEnabled(isBNS || isBHNS); 

    if (!isBNS) ui->lineEdit_LambdaAl2->setText("0.0");
    if (!isBNS && !isBHNS) ui->lineEdit_LambdaBl2->setText("0.0");
}

void MainWindow::on_tableWidget_PlotSettings_cellClicked(int row, int col)
{
    const int COLOR_COLUMN = 2; 
    if (col == COLOR_COLUMN) {
        // Блокируем сигналы, чтобы окно не открылось дважды от случайного даблклика
        ui->tableWidget_PlotSettings->blockSignals(true);
        
        QTableWidgetItem *item = ui->tableWidget_PlotSettings->item(row, col);
        if (item) {
            QColor initialColor = item->background().color();
            QColor newColor = QColorDialog::getColor(initialColor, this, "Выберите цвет орбиты");
            
            if (newColor.isValid()) {
                item->setBackground(newColor);
                item->setText(""); 
                // Сразу обновляем график
                on_pushButtonApplySettings_clicked(); 
            }
        }
        ui->tableWidget_PlotSettings->blockSignals(false);
    }
}

void MainWindow::on_dekart_cell_changed(int row, int col)
{
    if (row < 0 || row >= static_cast<int>(body_system_XYZ.size())) return;

    ui->tableWidget_info_dekart->blockSignals(true);
    
    QTableWidgetItem* item = ui->tableWidget_info_dekart->item(row, col);
    if (!item) {
        ui->tableWidget_info_dekart->blockSignals(false);
        return;
    }

    QString newText = item->text();
    bool ok = true;
    double newValue = 0.0;
    
    if (col != 0) { 
        newValue = newText.toDouble(&ok);
    }

    switch (col) {
        case 0: 
            labels_XYZ[row] = newText.toStdString();
            syncPlotSettingsTable();
            break;
        case 1: 
            if (ok && newValue >= 0) body_system_XYZ[row].m = newValue;
            else item->setText(QString::number(body_system_XYZ[row].m));
            break;
        case 2: if (ok) body_system_XYZ[row].x = newValue; else item->setText(QString::number(body_system_XYZ[row].x)); break;
        case 3: if (ok) body_system_XYZ[row].y = newValue; else item->setText(QString::number(body_system_XYZ[row].y)); break;
        case 4: if (ok) body_system_XYZ[row].z = newValue; else item->setText(QString::number(body_system_XYZ[row].z)); break;
        case 5: if (ok) body_system_XYZ[row].vx = newValue; else item->setText(QString::number(body_system_XYZ[row].vx)); break;
        case 6: if (ok) body_system_XYZ[row].vy = newValue; else item->setText(QString::number(body_system_XYZ[row].vy)); break;
        case 7: if (ok) body_system_XYZ[row].vz = newValue; else item->setText(QString::number(body_system_XYZ[row].vz)); break;
    }

    ui->tableWidget_info_dekart->blockSignals(false);
    
    if (!ok && col != 0) {
        statusBar()->showMessage("Ошибка: введите число!", 3000);
    } else {
        statusBar()->showMessage("Данные тела обновлены", 2000);
    }
}

// Таблица настроек
void MainWindow::syncPlotSettingsTable() {

    ui->tableWidget_PlotSettings->setRowCount(0);

    std::vector<std::string> allLabels;
    if (ui->tabWidget->currentIndex() == 0) {
        allLabels = labels_XYZ;
    } else {
        allLabels = labels_Kepler;
    }

    // Стандартные цвета для инициализации
    std::vector<QColor> defaultColors = {
        Qt::yellow, Qt::gray, Qt::red, Qt::blue, 
        QColor(255, 140, 0), QColor(210, 180, 140), 
        Qt::magenta, Qt::cyan, QColor(0, 0, 139)
    };

    for (size_t i = 0; i < allLabels.size(); ++i) {
        int row = ui->tableWidget_PlotSettings->rowCount();
        ui->tableWidget_PlotSettings->insertRow(row); 

        // Колонка видимости
        QTableWidgetItem* checkItem = new QTableWidgetItem();
        checkItem->setFlags(Qt::ItemIsUserCheckable | Qt::ItemIsEnabled);
        checkItem->setCheckState(Qt::Checked);
        ui->tableWidget_PlotSettings->setItem(row, 0, checkItem); 

        // Колонка легенды
        QTableWidgetItem* nameItem = new QTableWidgetItem(QString::fromStdString(allLabels[i]));
        ui->tableWidget_PlotSettings->setItem(row, 1, nameItem); 

        // Колонка цвета
        QTableWidgetItem* colorItem = new QTableWidgetItem();
        QColor col = (i < defaultColors.size()) ? defaultColors[i] : Qt::white;
        colorItem->setBackground(col);
        ui->tableWidget_PlotSettings->setItem(row, 2, colorItem); 
    }
}

// Динамическое сокращение путей
void MainWindow::updatePathLabels() {

    if (!pathFolder.empty()) {
        QString fullPath = QString::fromStdString(pathFolder);
        if (fullPath.endsWith('/')) fullPath.chop(1); 
        
        QFontMetrics metrics(ui->lbl_path_Folder->font());
        QString elidedText = metrics.elidedText(fullPath, Qt::ElideMiddle, ui->lbl_path_Folder->width());
        ui->lbl_path_Folder->setText(elidedText);
        ui->lbl_path_Folder->setToolTip(fullPath); // Подсказка с полным путем
    }

    // Для файла визуализации
    if (!pathFileVisual.empty()) {
        QString fullPath = QString::fromLocal8Bit(pathFileVisual.c_str());
        QFontMetrics metrics(ui->lbl_path_FileVisual->font());
        QString elidedText = metrics.elidedText(fullPath, Qt::ElideMiddle, ui->lbl_path_FileVisual->width());
        ui->lbl_path_FileVisual->setText(elidedText);
        ui->lbl_path_FileVisual->setToolTip(fullPath); // Подсказка с полным путем
    }
}

void MainWindow::resizeEvent(QResizeEvent *event) {
    QMainWindow::resizeEvent(event);
    updatePathLabels(); 
}

void MainWindow::on_gw_slider_tail_valueChanged(int value) 
{
    if (ui->gw_lbl_tail) {
        ui->gw_lbl_tail->setText(QString("Отображение: %1%").arg(value));
    }
    
    if (last_x1.isEmpty()) return;

    if (!vis) {
        std::vector<double> empty;
        vis = new Visual(body_system_XYZ, empty, empty, empty);
    }
    QString name1 = ui->gw_lineEdit_name1->text().isEmpty() ? "Объект 1" : ui->gw_lineEdit_name1->text();
    QString name2 = ui->gw_lineEdit_name2->text().isEmpty() ? "Объект 2" : ui->gw_lineEdit_name2->text();

    // Считываем текущий размер шрифта из интерфейса
    int currentFontSize = ui->gw_spin_fontSize->value();

    vis->PlotGWDynamics(ui->gw_plot_dynamics, last_x1, last_y1, last_x2, last_y2, value, name1, name2, currentFontSize);
    vis->PlotGWWaveform(ui->gw_plot_waveform, last_t, last_hp, last_hx, currentFontSize);
}

void MainWindow::on_gw_chk_showWave_stateChanged(int state) {
    if (ui->gw_plot_waveform->graphCount() > 0) {
        ui->gw_plot_waveform->graph(0)->setVisible(state == Qt::Checked);
        ui->gw_plot_waveform->replot(); 
    }
}

void MainWindow::on_gw_combo_waveColor_currentIndexChanged(const QString &color) {
    if (ui->gw_plot_waveform->graphCount() > 0) {
        ui->gw_plot_waveform->graph(0)->setPen(QPen(QColor(color)));
        ui->gw_plot_waveform->replot(); 
    }
}

void MainWindow::on_gw_spin_lineWidth_valueChanged(int value) {
    if (ui->gw_plot_waveform->graphCount() > 0) {
        QPen pen = ui->gw_plot_waveform->graph(0)->pen();
        pen.setWidth(value);
        ui->gw_plot_waveform->graph(0)->setPen(pen);
        ui->gw_plot_waveform->replot(); 
    }
}

void MainWindow::on_gw_combo_plotLayout_currentIndexChanged(int index) {
    switch(index) {
        case 0: // Оба графика
            ui->gw_plot_dynamics->setVisible(true);
            ui->gw_plot_waveform->setVisible(true);
            ui->gw_plot_splitter->setSizes(QList<int>({200, 200})); 
            break;
        case 1: // Только динамика
            ui->gw_plot_dynamics->setVisible(true);
            ui->gw_plot_waveform->setVisible(false);
            break;
        case 2: // Только волна
            ui->gw_plot_dynamics->setVisible(false);
            ui->gw_plot_waveform->setVisible(true);
            break;
    }
}

void MainWindow::on_btnBuildStudyPlot_clicked()
{
    gsl_set_error_handler_off();
    const double totalM = 60.0; 
    QVector<double> massRatios, maxAmplitudes;

    ui->gw_plot_waveform->clearGraphs();
    QCPGraph *studyGraph = ui->gw_plot_waveform->addGraph();
    studyGraph->setPen(QPen(Qt::blue, 2)); // Сделаем синим для отличия

    // Настройка осей (теперь всё будет подписано)
    QFont font = ui->gw_plot_waveform->xAxis->labelFont();
    int currentFontSize = ui->gw_spin_fontSize->value();
    font.setPointSize(currentFontSize);
    ui->gw_plot_waveform->xAxis->setLabelFont(font);
    ui->gw_plot_waveform->yAxis->setLabelFont(font);
    ui->gw_plot_waveform->xAxis->setTickLabelFont(font);
    ui->gw_plot_waveform->yAxis->setTickLabelFont(font);
    if (ui->gw_plot_waveform->legend) ui->gw_plot_waveform->legend->setFont(font);
    
    ui->gw_plot_waveform->xAxis->setLabel("Отношение масс (m2/m1)");
    ui->gw_plot_waveform->yAxis->setLabel("Максимальная амплитуда |h|");
    ui->gw_plot_waveform->xAxis->setRange(0, 1); // Твое требование по диапазону

    EOBParameters *old_pars = EOBPars;

    // Цикл по отношению масс от 0.1 до 1.0 
    for (double r = 0.1; r <= 1.0; r += 0.05) 
    {
        EOBParameters *pars = nullptr;
        EOBParameters_alloc(&pars);
        
        int sys_type = ui->gw_combo_systemType->currentIndex();
        EOBParameters_defaults(sys_type, pars);
        
        pars->binary = sys_type;
        pars->M = totalM;
        
        pars->q = 1.0 / r; 
        
        double M_sec = pars->M * 4.9254910255435759e-06; 
        pars->initial_frequency = 20.0 * M_sec; 

        pars->chi1z = 0.0;
        pars->chi2z = 0.0;
        pars->distance = 100.0;
        pars->domain = 0; 
        pars->use_geometric_units = 1;
        pars->output_dynamics = 0;
        pars->output_hpc = 0;

        EOBPars = pars;
        if (eob_set_params(0, 0) != 0) {
            EOBParameters_free(pars);
            continue;
        }

        Waveform *hpc = nullptr;
        WaveformFD *hfpc = nullptr;
        Waveform_lm *hmodes = nullptr;
        WaveformFD_lm *hfmodes = nullptr;
        Dynamics *dynf = nullptr;
        Waveform_lm *hTmodes = nullptr;
        Waveform_lm *hTmmodes = nullptr;
        Waveform_lm *hT0modes = nullptr;
        WaveformFD_lm *hfTmodes = nullptr;

        extern NQCdata *NQC;
        NQC = nullptr;

        int status = EOBRun(&hpc, &hfpc, &hmodes, &hfmodes, &dynf, 
                            &hTmodes, &hTmmodes, &hT0modes, &hfTmodes, 0, 1);

        if (status == 0 && hpc && hpc->size > 0) 
        {
            double currentMaxAmp = 0.0;
            for (int i = 0; i < hpc->size; ++i) {
                double hp = hpc->real[i];
                double hx = hpc->imag[i];
                double amp = std::sqrt(hp*hp + hx*hx);
                if (amp > currentMaxAmp) currentMaxAmp = amp;
            }
            
            // Сохраняем r (0.1...1.0) для оси Ox
            massRatios.append(r); 
            maxAmplitudes.append(currentMaxAmp);
        }

        if (hpc) Waveform_free(hpc);
        if (hfpc) WaveformFD_free(hfpc);
        if (hmodes) Waveform_lm_free(hmodes);
        if (hfmodes) WaveformFD_lm_free(hfmodes);
        if (dynf) Dynamics_free(dynf);
        if (hTmodes) Waveform_lm_free(hTmodes);
        if (hTmmodes) Waveform_lm_free(hTmmodes);
        if (hT0modes) Waveform_lm_free(hT0modes);
        if (hfTmodes) WaveformFD_lm_free(hfTmodes);

        EOBParameters_free(pars);
        QCoreApplication::processEvents();
    }

    EOBPars = old_pars;

    studyGraph->setData(massRatios, maxAmplitudes);
    ui->gw_plot_waveform->rescaleAxes();
    ui->gw_plot_waveform->replot();
}