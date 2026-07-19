cmake_minimum_required(VERSION 3.24)

# Компиляторы
set(CMAKE_CXX_COMPILER "C:/Qt/Tools/mingw1310_64/bin/g++.exe")
set(CMAKE_C_COMPILER "C:/Qt/Tools/mingw1310_64/bin/gcc.exe")

project(NIR_Project LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

# --- ФЛАГИ СТАТИЧЕСКОЙ КОМПИЛЯЦИИ ---
if(WIN32)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static -static-libgcc -static-libstdc++")
    
    set(VCPKG_TARGET_TRIPLET "x64-mingw-static")
    set(VCPKG_HOST_TRIPLET "x64-mingw-static")
endif()
# ------------------------------------

set(CMAKE_PREFIX_PATH "C:/Qt-6.9.2-ninja")

find_package(Qt6 REQUIRED COMPONENTS Core Gui Widgets Concurrent PrintSupport OpenGL)
find_package(GSL REQUIRED)

# Путь к исходникам TEOBResumS
set(TEOB_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../teobresums/src")

# Список файлов TEOBResumS
set(TEOB_SOURCES
    "${TEOB_SRC_DIR}/SpecialFuns.c"
    "${TEOB_SRC_DIR}/TEOBResumS.c"
    "${TEOB_SRC_DIR}/TEOBResumSDynamics.c"
    "${TEOB_SRC_DIR}/TEOBResumSFits.c"
    "${TEOB_SRC_DIR}/TEOBResumSFlux.c"
    "${TEOB_SRC_DIR}/TEOBResumSInitialConditions.c"
    "${TEOB_SRC_DIR}/TEOBResumSMetric.c"
    "${TEOB_SRC_DIR}/TEOBResumSOMP.c"
    "${TEOB_SRC_DIR}/TEOBResumSPostAdiab.c"
    "${TEOB_SRC_DIR}/TEOBResumSUtils.c"
    "${TEOB_SRC_DIR}/TEOBResumSPars.c"
    "${TEOB_SRC_DIR}/TEOBResumSWaveform.c"
)

# Сборка библиотеки TEOBResumS
add_library(teobresums STATIC ${TEOB_SOURCES})

# target_include_directories(teobresums PUBLIC "${TEOB_SRC_DIR}" ${GSL_INCLUDE_DIRS})
target_link_libraries(teobresums PUBLIC GSL::gsl GSL::gslcblas)

# Для статики важно, чтобы макросы соответствовали
target_compile_definitions(teobresums PRIVATE NOTMAIN _GNU_SOURCE)

# Настройка основного приложения
qt_add_executable(myapp WIN32
    main.cpp
    algorithm.cpp
    system.cpp
    visualisation.cpp
    MainWindow.cpp
    qcustomplot.cpp
    qcustomplot.h
    algorithm.h
    system.h
    visualisation.h
    struct.h
    MainWindow.h
)

target_include_directories(myapp PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}" "${TEOB_SRC_DIR}" ${GSL_INCLUDE_DIRS})

target_link_libraries(myapp PUBLIC 
    GSL::gsl 
    GSL::gslcblas 
    teobresums 
    Qt6::Core 
    Qt6::Gui 
    Qt6::Widgets 
    Qt6::Concurrent 
    Qt6::PrintSupport
    Qt6::OpenGL
)