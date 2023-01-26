#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QFileDialog>
#include <QMainWindow>

#include "subdivision/subdivider.h"

namespace Ui {
class MainWindow;
}

/**
 * @brief The MainWindow class represents the main window.
 */
class MainWindow : public QMainWindow {
    Q_OBJECT

    public:
        explicit MainWindow(QWidget *parent = nullptr);
        ~MainWindow() override;


    private slots:
        void on_LoadOBJ_pressed();
        void on_MeshPresetComboBox_currentTextChanged(const QString &meshName);
        void on_SubdivSteps_valueChanged(int subdivLevel);
        void on_sharpnessSliderValue_valueChanged(int value);

        void timeout();


        void on_decimalSharpnessSpinBox_valueChanged(double arg1);

        void on_edgeSelectionRadioButton_toggled(bool checked);

private:
        void importOBJ(const QString &fileName);

        Ui::MainWindow *ui;
        Subdivider *subdivider;
};

#endif  // MAINWINDOW_H
