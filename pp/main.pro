QT       += core gui
greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport

greaterThan(QT_MAJOR_VERSION, 4): CONFIG += c++11
lessThan(QT_MAJOR_VERSION, 5): QMAKE_CXXFLAGS += -std=c++11 -O3 -mtune=native -march=native -mfpmath=both

TARGET = main
TEMPLATE = app

SOURCES += main.cpp\
           mainwindow.cpp \
           dynamical_system.cpp \
           sys_func.cpp \
           qcustomplot.cpp

HEADERS += mainwindow.h \
           dynamical_system.hpp \
           qcustomplot.h

FORMS   += mainwindow.ui