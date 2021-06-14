#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <qgraphicsscene.h>
#include <QGraphicsTextItem>
#include <QDebug>
#include <QString>
#include <math.h>
#include <QWheelEvent>
#include <QFile>
#include <QFileDialog>
#include "AudioFile.h"
#include "./FreqDetect/FFT_for_cpp.h"
#include "./FreqDetect/FFT_for_cpp_emxAPI.h"
#include "./FreqDetect/FFT_for_cpp_emxutil.h"
#include "./FreqDetect/FFT_for_cpp_types.h"
#include <windows.h>
#include <mmsystem.h>
#include <conio.h>
#include <stdlib.h>
#include <QListWidgetItem>

#ifndef OLDCPP
   #include <iostream>
   using namespace std;
#else
   #include <iostream.h>
#endif

#include "json.hpp"
//using namespace peti;
using json=nlohmann::json;

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE


enum display_form {
    dist,
    dB,
    sine
};

enum jelek {
    sinic,

};

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();


public slots:
    void hang(float freq);

    void wheelEvent(QWheelEvent * event);

private slots:
    void on_pushButton_clicked();

    void on_DisplayOption_currentRowChanged(int currentRow);

    void on_actionPlay_audio_file_triggered();

private:
    Ui::MainWindow *ui;
    QGraphicsScene *scene;
    QGraphicsTextItem *text;
    json options;
    const float guitar_freq[6]={
        329.63,
        246.94,
        196,
        146.83,
        110.0,
        82.41
    };
    float zoom=25;
    float last_freq;
    const float minimal_space=30,maximal_space=70;
    int steps=1;
    QString last_path=".";
    display_form disp=dist;
    emxArray_real_T * dBs;
    emxArray_real_T * freqs;
    AudioFile<float> file;
    vector<double> sample_data;
    float closest;
    bool first_button=false,blocker=false;
};
#endif // MAINWINDOW_H
