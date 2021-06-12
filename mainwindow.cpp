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
    scene->clear();
    double fele=ui->graphicsView->height()/2;
    scene->addLine(0,fele,ui->graphicsView->width(),fele,QPen(Qt::black));
    double freq_lenth=std::max(guitar_freq[0],freq);
    double graph_length=(double)(ui->graphicsView->width())*0.80;
    graph_length=(double)(ui->graphicsView->width());
    scene->addLine(ui->graphicsView->width()/2,ui->graphicsView->height(),ui->graphicsView->width()/2,0,QPen(Qt::red));
    int k=0;
    float delta_x=graph_length / (zoom*2) * (steps*2)-graph_length / (zoom*2) * steps;
    if (delta_x<minimal_space) {
        if (steps==1) {
            steps+=4;
        } else {
            steps+=5;
        }
    } else if (delta_x>maximal_space) {
        if (steps!=5 && steps>1) {
            steps-=5;
        } else {
            steps=1;
        }
    }
    for (auto i=0.0f;i<=zoom*2;i+=steps) {
        double x = graph_length / (zoom*2) * i;
        if (i != zoom) {
            scene->addLine(x, fele + 15, x, fele - 15, QPen(Qt::black));
        }
        text = new QGraphicsTextItem();
        text->setPos(x - 10, fele - 5 + (25 * pow(-1, k)));
        text->setPlainText(QString::number(freq - zoom + i) + "Hz");
        text->setScale(0.8);
        scene->addItem(text);
        k++;
    }
    k=0;
    for (auto i:guitar_freq) {
        if (i >= freq - zoom && i <= freq + zoom) {
            double x = graph_length / (zoom*2) * (i-freq+zoom);
            scene->addLine(x, fele + 55, x, fele - 55, QPen(Qt::green));
            text = new QGraphicsTextItem();
            text->setPos(x - 20, fele - 5 + (80 * pow(-1, k)));
            text->setPlainText(QString::number(i) + "Hz");
            text->setScale(0.8);
            scene->addItem(text);
            k++;
        }
    }/*
    } else {
        double x_1=graph_length/freq_lenth*freq+(double)(ui->graphicsView->width())*0.1;
        scene->addLine(x_1,ui->graphicsView->height(),x_1,0,QPen(Qt::red));
        for (auto i:guitar_freq) {
            double x=graph_length/freq_lenth*i+(double)(ui->graphicsView->width())*0.1;
            scene->addLine(x,ui->graphicsView->height(),x,0,QPen(Qt::black));
        }
    }*/
}

void MainWindow::wheelEvent(QWheelEvent *event)
{
    if (event->angleDelta().y()>0) { //zoom in
        zoom-=event->angleDelta().y()/100;
        hang(last_freq);
    } else if (event->angleDelta().y()<0) {
        zoom-=event->angleDelta().y()/100;
        hang(last_freq);
    } //zoom out
}

void MainWindow::on_pushButton_clicked()
{
    float freq=ui->doubleSpinBox->value();
    hang(freq);
    last_freq=freq;
    ui->frequency_visual->setText(QString::number(freq));
    float closest=guitar_freq[0];
    for (auto i:guitar_freq) {
        if (abs(freq-i)<abs(closest-freq)) {
            closest=i;
        }
    }
    ui->label_4->setText(QString::number(closest));
    /*scene->addLine(0,ui->graphicsView->height()/2,ui->graphicsView->width(),ui->graphicsView->height()/2,QPen(Qt::black));
    scene->addRect(100,100,50,50,QPen(Qt::black));*/
}
