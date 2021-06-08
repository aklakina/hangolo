#include "mainwindow.h"
#include "ui_mainwindow.h"


MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    scene = new QGraphicsScene();
    ui->graphicsView->setScene(scene);
    scene->setSceneRect(ui->graphicsView->rect());
    scene->setBackgroundBrush(QBrush(Qt::white));
    ui->graphicsView->setSceneRect(scene->sceneRect());
    ui->graphicsView->setAlignment(Qt::AlignTop | Qt::AlignLeft);
    ui->graphicsView->setScene(scene);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::hang(float freq) {
    double fele=ui->graphicsView->height()/2;
    scene->addLine(0,fele,ui->graphicsView->width(),fele,QPen(Qt::black));
    double freq_lenth=std::max(guitar_freq[0],freq);
    double graph_length=(double)(ui->graphicsView->width())*0.80;

    if (freq_lenth>500) {
        graph_length=(double)(ui->graphicsView->width());
        scene->addLine(ui->graphicsView->width()/2,ui->graphicsView->height(),ui->graphicsView->width()/2,0,QPen(Qt::red));
        int suruseg=50;
        for (auto i=0;i<suruseg;i++) {
            double x=graph_length/suruseg*i;
            if (i!=50) {
                scene->addLine(x,fele+15,x,fele-15,QPen(Qt::black));
            }
            text=new QGraphicsTextItem();
            text->setPos(x-10,fele-5+(25*pow(-1,i)));
            text->setPlainText(QString::number(freq-suruseg/2+i)+"Hz");
            text->setScale(0.8);
            scene->addItem(text);
        }
    } else {
        double x_1=graph_length/freq_lenth*freq+(double)(ui->graphicsView->width())*0.1;
        scene->addLine(x_1,ui->graphicsView->height(),x_1,0,QPen(Qt::red));
        for (auto i:guitar_freq) {
            double x=graph_length/freq_lenth*i+(double)(ui->graphicsView->width())*0.1;
            scene->addLine(x,ui->graphicsView->height(),x,0,QPen(Qt::black));
        }
    }
}

void MainWindow::on_pushButton_clicked()
{
    scene->clear();
    hang(ui->doubleSpinBox->value());
    /*scene->addLine(0,ui->graphicsView->height()/2,ui->graphicsView->width(),ui->graphicsView->height()/2,QPen(Qt::black));
    scene->addRect(100,100,50,50,QPen(Qt::black));*/
}
