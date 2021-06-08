#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <qgraphicsscene.h>
#include <QGraphicsTextItem>
#include <QDebug>
#include <QString>
#include <math.h>

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

public slots:
    void hang(float freq);

private slots:
    void on_pushButton_clicked();

private:
    Ui::MainWindow *ui;
    QGraphicsScene *scene;
    QGraphicsTextItem *text;
    const float guitar_freq[6]={
        329.63,
        246.94,
        196,
        146.83,
        110.0,
        82.41
    };
};
#endif // MAINWINDOW_H
