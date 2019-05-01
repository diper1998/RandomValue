#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "magic.h"
#include "ui_mainwindow.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

    magic WIZARD;
    int numbGraph = 0;
    int numbHistogram = 0;
    int numb = 0;
    void MainWindow::WriteTable(int i, double X, double Y);
    void MainWindow::WriteTable_1(int i, double X, double Y, double Z, double Q);
    void MainWindow::WriteTable_2(int i, double X, double Y, double Z, double Q, double V, double W, double O, double L, double H);
    void MainWindow::WriteTable_3(int i, double X, double Y, double Z, double Q);
    QTableWidgetItem* MainWindow::CreateItem(QVariant it);
    void MainWindow::PaintGraph(double* myX, double* myY, int N, const int nGraph, QString nameGraph);
    void MainWindow::PaintGraph_1(double* myX, double* myY, int N, const int nGraph, QString nameGraph);


private slots:
    void on_BUTTON_clicked();

    void on_pushButton_CLEAN_clicked();



    void on_pushButton_THISRV_clicked();

private:
    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
