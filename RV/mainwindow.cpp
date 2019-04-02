#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <iostream>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    this->setWindowTitle("Random");

    ui->GRAPH->setInteraction(QCP::iRangeZoom, true); //удаление и приближение
    ui->GRAPH->setInteraction(QCP::iRangeDrag, true); //перетаскивание

    ui->GRAPH->xAxis->setLabel("X");
    ui->GRAPH->yAxis->setLabel("Y");

    ui->GRAPH_1->setInteraction(QCP::iRangeZoom, true); //удаление и приближение
    ui->GRAPH_1->setInteraction(QCP::iRangeDrag, true); //перетаскивание

    ui->GRAPH_1->xAxis->setLabel("X");
    ui->GRAPH_1->yAxis->setLabel("Y");
}

MainWindow::~MainWindow()
{


    delete ui;
}

void MainWindow::on_BUTTON_clicked()
{
    WIZARD.n = ui->lineEdit_n->text().toInt();
    WIZARD.k = ui->lineEdit_k->text().toInt();
    WIZARD.a = ui->lineEdit_a->text().toDouble();
    WIZARD.sigma = ui->lineEdit_Sigma->text().toDouble();
    double begin = ui->lineEdit_Begin->text().toDouble();
    double end = ui->lineEdit_End->text().toDouble();

    WIZARD.SetMassifs();


    srand(time(0));
    ui->TABLE->setRowCount(0);
    ui->TABLE_1->setRowCount(0);

    for(int i = 0; i < WIZARD.n; i++){
        ui->TABLE->insertRow(i);
        WIZARD.rV = (double)rand()/RAND_MAX;
                while(WIZARD.rV == 0){
                WIZARD.rV = (double)rand()/RAND_MAX;
    }
        WIZARD.SetRandomValue(i,WIZARD.rV);
        WriteTable(i,i,WIZARD.randomValue[i]);

    }


    ui->TABLE->sortItems(1);

    std::sort(WIZARD.randomValue, WIZARD.randomValue+(WIZARD.n));

    WIZARD.SetAverage();
    WIZARD.SetDispersion();
    WIZARD.SetRange();
    WIZARD.SetMedian();

    WIZARD.SetDistibution();
    WIZARD.SetSampleDistibution();
    WIZARD.SetD();
    WIZARD.SetDS();
    WIZARD.SetHistogram();
    WIZARD.SetZ();
    WIZARD.SetDensity();
    WIZARD.SetDensityH();
    WIZARD.SetQ();
   // WIZARD.SetR();
    //WIZARD.SetFR();
    WIZARD.SetDen();
    WIZARD.SetCheck(begin, end);
    WIZARD.SetFR();

     WIZARD.SetTmp();

    for(int i = 0; i < WIZARD.k; i++){
        ui->TABLE_1->insertRow(i);
        WriteTable_1(i,WIZARD.z[i], WIZARD.f[i], WIZARD.h[i], WIZARD.densityH[i]);

    }

    for(int i = 0; i < WIZARD.k; i++){
        ui->TABLE_3->insertRow(i);
        WriteTable_3(i, WIZARD.mCheck[i][0], WIZARD.mCheck[i][1], WIZARD.q[i]);

    }

/*
    for(int i = 0; i < WIZARD.k; i++){
        ui->TABLE_1->insertRow(i);
        WriteTable_1(i,WIZARD.z[i], WIZARD.f[i], WIZARD.h[i], WIZARD.densityH[i]);

    }

*/

    ui->label_Range->setText("R^ (Range) = " + QVariant(WIZARD.range).toString());
    ui->label_D->setText("D (max|Fj-F^j|) = " + QVariant(WIZARD.D).toString());
    ui->label_DS->setText("|D - S^2| = " + QVariant(WIZARD.DS).toString());
    ui->label_maxDensityH->setText("max|fj - nj/n*dj| = " + QVariant(WIZARD.maxDensityH).toString());
    ui->label_R->setText("R0 = " + QVariant(WIZARD.R).toString());
    ui->label_FR->setText("F(R) = " + QVariant(WIZARD.FR).toString());

    if(WIZARD.FR < WIZARD.a){
        ui->label_Hypothesis->setText("Hypothesis: unacceptable");
    } else {
        ui->label_Hypothesis->setText("Hypothesis: accepted");
    }


    PaintGraph(WIZARD.randomValueD, WIZARD.distribution_, 1000, numbGraph, "F"+QVariant(numbGraph).toString());
    numbGraph++;
     PaintGraph(WIZARD.randomValue_, WIZARD.sampleDistribution_, WIZARD.n*2, numbGraph, "F^"+QVariant(numbGraph).toString());
  //  PaintGraph(WIZARD.randomValue, WIZARD.sampleDistribution, WIZARD.n, numbGraph, "F^"+QVariant(numbGraph).toString());
    numbGraph++;

    ui->TABLE_2->insertRow(numb);
    WriteTable_2(numb, WIZARD.average, pow(M_PI/2,0.5)*WIZARD.sigma, abs( pow(M_PI/2,0.5)*WIZARD.sigma-WIZARD.average),
                 WIZARD.dispersion, (2-M_PI/2)*WIZARD.sigma*WIZARD.sigma, abs(WIZARD.dispersion - (2-M_PI/2)*WIZARD.sigma*WIZARD.sigma),
                 WIZARD.median, WIZARD.sigma*pow(log(4),0.5), abs(WIZARD.median - WIZARD.sigma*pow(log(4), 0.5)));
    numb++;




    //PaintGraph_1(WIZARD.randomValue, WIZARD.rectangle, WIZARD.n, numbHistogram, "H"+QVariant(numbHistogram).toString());
    //numbHistogram++;
    PaintGraph_1(WIZARD.tmp, WIZARD.tmpF, WIZARD.k*2, numbHistogram, "H"+QVariant(numbHistogram).toString());
    numbHistogram++;

    PaintGraph_1(WIZARD.randomValueD, WIZARD.density, 1000, numbHistogram, "H"+QVariant(numbHistogram).toString());
    numbHistogram++;



}


QTableWidgetItem* MainWindow::CreateItem(QVariant it){
    QTableWidgetItem *item = new QTableWidgetItem();
    item->setText(it.toString());
    return item;
}

void MainWindow::WriteTable(int i, double X, double Y){


    QVariant qX = X;
    QVariant qY = Y;


    ui->TABLE->setItem(i, 0, CreateItem(qX));
    ui->TABLE->setItem(i, 1, CreateItem(qY));


}

void MainWindow::WriteTable_1(int i, double X, double Y, double Z, double Q){


    QVariant qX = X;
    QVariant qY = Y;
    QVariant qZ = Z;
    QVariant qQ = Q;


    ui->TABLE_1->setItem(i, 0, CreateItem(qX));
    ui->TABLE_1->setItem(i, 1, CreateItem(qY));
    ui->TABLE_1->setItem(i, 2, CreateItem(qZ));
    ui->TABLE_1->setItem(i, 3, CreateItem(qQ));


}


void MainWindow::WriteTable_3(int i, double X, double Y, double Z){


    QVariant qX = X;
    QVariant qY = Y;
    QVariant qZ = Z;



    ui->TABLE_3->setItem(i, 0, CreateItem(qX));
    ui->TABLE_3->setItem(i, 1, CreateItem(qY));
    ui->TABLE_3->setItem(i, 2, CreateItem(qZ));



}


void MainWindow::WriteTable_2(int i, double X, double Y, double Z, double Q, double V, double W, double O, double L, double H){


    QVariant qX = X;
    QVariant qY = Y;
    QVariant qZ = Z;
    QVariant qQ = Q;
    QVariant qV = V;
    QVariant qW = W;
    QVariant qO = O;
    QVariant qL = L;
    QVariant qH = H;


    ui->TABLE_2->setItem(i, 0, CreateItem(qX));
    ui->TABLE_2->setItem(i, 1, CreateItem(qY));
    ui->TABLE_2->setItem(i, 2, CreateItem(qZ));
    ui->TABLE_2->setItem(i, 3, CreateItem(qQ));
    ui->TABLE_2->setItem(i, 4, CreateItem(qV));
    ui->TABLE_2->setItem(i, 5, CreateItem(qW));
    ui->TABLE_2->setItem(i, 6, CreateItem(qO));
    ui->TABLE_2->setItem(i, 7, CreateItem(qL));
    ui->TABLE_2->setItem(i, 8, CreateItem(qH));


}


void MainWindow::PaintGraph(double* myX, double* myY, int N, const int nGraph, QString nameGraph){


    QVector<double> X(N);
    QVector<double> Y(N);

    for(int i = 0; i < N; i++){
        X[i] = myX[i];
        Y[i] = myY[i];
    }

    ui->GRAPH->legend->setVisible(true);
    ui->GRAPH->legend->setFont(QFont("Helvetica", 9));
    ui->GRAPH->legend->setRowSpacing(-3);

    QVector<QCPScatterStyle::ScatterShape> shapes;
    shapes << QCPScatterStyle::ssCross;
    shapes << QCPScatterStyle::ssPlus;
    shapes << QCPScatterStyle::ssCircle;
    shapes << QCPScatterStyle::ssDisc;
    shapes << QCPScatterStyle::ssSquare;
    shapes << QCPScatterStyle::ssDiamond;
    shapes << QCPScatterStyle::ssStar;
    shapes << QCPScatterStyle::ssTriangle;
    shapes << QCPScatterStyle::ssTriangleInverted;
    shapes << QCPScatterStyle::ssCrossSquare;
    shapes << QCPScatterStyle::ssPlusSquare;
    shapes << QCPScatterStyle::ssCrossCircle;
    shapes << QCPScatterStyle::ssPlusCircle;
    shapes << QCPScatterStyle::ssPeace;
    shapes << QCPScatterStyle::ssCustom;

    QPen pen;


    pen.setColor(QColor(qSin(nGraph%14*0.3)*100+100, qSin(nGraph%14*0.6+0.7)*100+100, qSin(nGraph%14*0.4+0.6)*100+100));


    ui->GRAPH->addGraph();
    ui->GRAPH->graph(nGraph)->setData(X, Y);

    ui->GRAPH->graph(nGraph)->setPen(pen);
    ui->GRAPH->graph(nGraph)->setName(nameGraph);
    ui->GRAPH->graph(nGraph)->setLineStyle(QCPGraph::lsLine);


    ui->GRAPH->graph(nGraph)->setPen(pen);
    ui->GRAPH->graph(nGraph)->setName(nameGraph);
    ui->GRAPH->graph(nGraph)->setLineStyle(QCPGraph::lsLine);
/*
    // set scatter style:
    if (shapes.at(numbGraph%14) != QCPScatterStyle::ssCustom)
     {
       ui->GRAPH->graph(numbGraph)->setScatterStyle(QCPScatterStyle(shapes.at(numbGraph%14), 10));
     }
     else
     {
       QPainterPath customScatterPath;
       for (int i=0; i<3; ++i)
         customScatterPath.cubicTo(qCos(2*M_PI*i/3.0)*9, qSin(2*M_PI*i/3.0)*9, qCos(2*M_PI*(i+0.9)/3.0)*9, qSin(2*M_PI*(i+0.9)/3.0)*9, 0, 0);
           ui->GRAPH->graph(numbGraph)->setScatterStyle(QCPScatterStyle(customScatterPath, QPen(Qt::black, 0), QColor(40, 70, 255, 50), 10));
     }
*/

     ui->GRAPH->graph(nGraph)->rescaleAxes();
     ui->GRAPH->replot();

   }

void MainWindow::PaintGraph_1(double* myX, double* myY, int N, const int nGraph, QString nameGraph){


    QVector<double> X(N);
    QVector<double> Y(N);

    for(int i = 0; i < N; i++){
        X[i] = myX[i];
        Y[i] = myY[i];
    }

    ui->GRAPH_1->legend->setVisible(true);
    ui->GRAPH_1->legend->setFont(QFont("Helvetica", 9));
    ui->GRAPH_1->legend->setRowSpacing(-3);

    QVector<QCPScatterStyle::ScatterShape> shapes;
    shapes << QCPScatterStyle::ssCross;
    shapes << QCPScatterStyle::ssPlus;
    shapes << QCPScatterStyle::ssCircle;
    shapes << QCPScatterStyle::ssDisc;
    shapes << QCPScatterStyle::ssSquare;
    shapes << QCPScatterStyle::ssDiamond;
    shapes << QCPScatterStyle::ssStar;
    shapes << QCPScatterStyle::ssTriangle;
    shapes << QCPScatterStyle::ssTriangleInverted;
    shapes << QCPScatterStyle::ssCrossSquare;
    shapes << QCPScatterStyle::ssPlusSquare;
    shapes << QCPScatterStyle::ssCrossCircle;
    shapes << QCPScatterStyle::ssPlusCircle;
    shapes << QCPScatterStyle::ssPeace;
    shapes << QCPScatterStyle::ssCustom;

    QPen pen;


    pen.setColor(QColor(qSin(nGraph%14*0.3)*100+100, qSin(nGraph%14*0.6+0.7)*100+100, qSin(nGraph%14*0.4+0.6)*100+100));


    ui->GRAPH_1->addGraph();
    ui->GRAPH_1->graph(nGraph)->setData(X, Y);

    ui->GRAPH_1->graph(nGraph)->setPen(pen);
    ui->GRAPH_1->graph(nGraph)->setName(nameGraph);
    ui->GRAPH_1->graph(nGraph)->setLineStyle(QCPGraph::lsLine);

    ui->GRAPH_1->graph(nGraph)->setPen(pen);
    ui->GRAPH_1->graph(nGraph)->setName(nameGraph);
    ui->GRAPH_1->graph(nGraph)->setLineStyle(QCPGraph::lsLine);
/*
    // set scatter style:
    if (shapes.at(numbGraph%14) != QCPScatterStyle::ssCustom)
     {
       ui->GRAPH->graph(numbGraph)->setScatterStyle(QCPScatterStyle(shapes.at(numbGraph%14), 10));
     }
     else
     {
       QPainterPath customScatterPath;
       for (int i=0; i<3; ++i)
         customScatterPath.cubicTo(qCos(2*M_PI*i/3.0)*9, qSin(2*M_PI*i/3.0)*9, qCos(2*M_PI*(i+0.9)/3.0)*9, qSin(2*M_PI*(i+0.9)/3.0)*9, 0, 0);
           ui->GRAPH->graph(numbGraph)->setScatterStyle(QCPScatterStyle(customScatterPath, QPen(Qt::black, 0), QColor(40, 70, 255, 50), 10));
     }
*/

     ui->GRAPH_1->graph(nGraph)->rescaleAxes();
     ui->GRAPH_1->replot();

   }

void MainWindow::on_pushButton_CLEAN_clicked()
{
    ui->TABLE->setRowCount(0);
    ui->TABLE_1->setRowCount(0);
    ui->TABLE_2->setRowCount(0);
    ui->TABLE_3->setRowCount(0);
    ui->GRAPH->clearGraphs();
    ui->GRAPH->replot();

    ui->GRAPH_1->clearGraphs();
    ui->GRAPH_1->replot();

    numbGraph = 0;
    numbHistogram = 0;
    numb = 0;


}
