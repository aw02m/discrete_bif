#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include "dynamical_system.hpp"
#include "ui_mainwindow.h"
#include "qcustomplot.h"
#include <QMainWindow>
#include <QTimer>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow {
  Q_OBJECT

public:
  explicit MainWindow(QWidget *parent = 0);
  ~MainWindow();

  void ppsetup();
  void ppcalc(QCustomPlot *customPlot);

  void setupPlayground(QCustomPlot *customPlot);

private slots:
  void ppSlot();
  void mousePress(QMouseEvent *event);

private:
  Ui::MainWindow *ui;
  QTimer dataTimer;
  QCPGraph *trajectory;
  dynamical_system ds;

protected:
  void keyPressEvent(QKeyEvent *event);
  // void mousePressEvent(QMouseEvent *event);
};

#endif // MAINWINDOW_H
