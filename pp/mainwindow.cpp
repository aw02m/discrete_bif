#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QDebug>
#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
#include <QDesktopWidget>
#endif
#include <QMessageBox>
#include <QMetaEnum>
#include <QScreen>
#include <iostream>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent), ui(new Ui::MainWindow),
      ds(QApplication::arguments().at(1).toStdString()) {
  ui->setupUi(this);
  setGeometry(400, 250, 768, 768);

  ppsetup();
}

void MainWindow::ppsetup() {
  ppcalc(ui->customPlot);
  setWindowTitle("Phase Portrait : pp");
  statusBar()->clearMessage();
  ui->customPlot->replot();
}

void MainWindow::ppcalc(QCustomPlot *customPlot) {
  trajectory = new QCPGraph(customPlot->xAxis, customPlot->yAxis);

  trajectory->setLineStyle(QCPGraph::lsNone);
  trajectory->setScatterStyle(QCPScatterStyle(
      QCPScatterStyle::ssDisc, QColor(0, 0, 0, 32), QColor(0, 0, 0, 32), 4));

  customPlot->xAxis->setLabel("axis[0]");
  customPlot->yAxis->setLabel("axis[1]");
  customPlot->xAxis->setRange(ds.xrange[0], ds.xrange[1]);
  customPlot->yAxis->setRange(ds.yrange[0], ds.yrange[1]);

  customPlot->setMouseTracking(true);

  connect(&dataTimer, SIGNAL(timeout()), this, SLOT(ppSlot()));
  connect(ui->customPlot, SIGNAL(mousePress(QMouseEvent *)), this,
          SLOT(mousePress(QMouseEvent *)));
  dataTimer.start(0); // Interval 0 means to refresh as fast as possible
}

void MainWindow::ppSlot() {
  double secs =
      QCPAxisTickerDateTime::dateTimeToKey(QDateTime::currentDateTime());
  Eigen::VectorXd x = ds.last_state;

  for (unsigned int i = 0; i < ds.max_plot; i++) {
    x = ds.function(x);
    ds.QCPGsol << QCPGraphData(x(ds.axis[0]), x(ds.axis[1]));
  }
  trajectory->data()->set(ds.QCPGsol, true);
  ui->customPlot->replot();
  ds.last_state = x;
  ds.QCPGsol.clear();

  // calculate frames per second:
  double key = secs;
  static double lastFpsKey;
  static int frameCount;
  ++frameCount;
  if (key - lastFpsKey > 2) // average fps over 2 seconds
  {
    ui->statusBar->showMessage(
        QString("%1 FPS, Total Data points: %2")
            .arg(frameCount / (key - lastFpsKey), 0, 'f', 0)
            .arg(trajectory->data()->size()),
        0);
    lastFpsKey = key;
    frameCount = 0;
  }
}

void MainWindow::mousePress(QMouseEvent *event) {
  if (event->button() == Qt::LeftButton) {
    ds.last_state(ds.axis[0]) =
        ui->customPlot->xAxis->pixelToCoord(event->pos().x());
    ds.last_state(ds.axis[1]) =
        ui->customPlot->yAxis->pixelToCoord(event->pos().y());
    // ds.QCPGsol.clear();
  }
}

void MainWindow::keyPressEvent(QKeyEvent *event) {
  static unsigned int param_index = 0;
  static Eigen::IOFormat Comma(8, 0, ", ", "\n", "[", "]");
  switch (event->key()) {
  case Qt::Key_Left:
    if (param_index == 0) {
      param_index = ds.p.size() - 1;
    } else {
      param_index--;
    }
    std::cout << "parameter index changed : " << param_index << std::endl;
    break;
  case Qt::Key_Right:
    if (param_index == ds.p.size() - 1) {
      param_index = 0;
    } else {
      param_index++;
    }
    std::cout << "parameter index changed : " << param_index << std::endl;
    break;
  case Qt::Key_Up:
    ds.p(param_index) += ds.dparams[param_index];
    std::cout << param_index << ":" << ds.p.transpose().format(Comma)
              << std::endl;
    break;
  case Qt::Key_Down:
    ds.p(param_index) -= ds.dparams[param_index];
    std::cout << param_index << ":" << ds.p.transpose().format(Comma)
              << std::endl;
    break;
  case Qt::Key_Space:
    std::cout << "param : " << ds.p.transpose().format(Comma) << std::endl;
    std::cout << "state : " << ds.last_state.transpose().format(Comma)
              << std::endl;
    break;
  case Qt::Key_X:
    if (ds.axis[0] == ds.xdim - 1) {
      ds.axis[0] = 0;
    } else {
      ds.axis[0]++;
    }
    std::cout << "x-axis changed : [" << ds.axis[0] << ", " << ds.axis[1] << "]"
              << std::endl;
    ds.QCPGsol.clear();
    break;
  case Qt::Key_Y:
    if (ds.axis[1] == ds.xdim - 1) {
      ds.axis[1] = 0;
    } else {
      ds.axis[1]++;
    }
    std::cout << "y-axis changed : [" << ds.axis[0] << ", " << ds.axis[1] << "]"
              << std::endl;
    ds.QCPGsol.clear();
    break;
  case Qt::Key_W:
    ds.save_json();
    std::cout << "wrote json output >>" << ds.json["json_out_path"]
              << std::endl;
    break;
  case Qt::Key_Q:
    QApplication::exit();
    break;
  }
}

void MainWindow::setupPlayground(QCustomPlot *customPlot) {
  Q_UNUSED(customPlot)
}

MainWindow::~MainWindow() { delete ui; }