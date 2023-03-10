#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QTimer>

#include "initialization/meshinitializer.h"
#include "initialization/objfile.h"
#include "subdivision/catmullclarksubdivider.h"
#include "subdivision/subdivider.h"


MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent), ui(new Ui::MainWindow) {
    ui->setupUi(this);
    ui->MeshGroupBox->setEnabled(ui->MainDisplay->settings.modelLoaded);
    ui->sharpnessSettings->setEnabled(ui->MainDisplay->settings.modelLoaded);
    QTimer *timer = new QTimer(this);
    connect(timer, &QTimer::timeout, this, QOverload<>::of(&MainWindow::timeout));
    timer->start(500);
}

MainWindow::~MainWindow() {
    delete ui;
    ui->MainDisplay->settings.meshes.clear();
    ui->MainDisplay->settings.meshes.squeeze();
}

/**
 * @brief MainWindow::importOBJ Imports an obj file and adds the constructed
 * half-edge to the collection of meshes.
 * @param fileName Path of the .obj file.
 */
void MainWindow::importOBJ(const QString& fileName) {
    OBJFile newModel = OBJFile(fileName);
    ui->MainDisplay->settings.meshes.clear();
    ui->MainDisplay->settings.meshes.squeeze();

    if (newModel.loadedSuccessfully()) {
        MeshInitializer meshInitializer;
        ui->MainDisplay->settings.meshes.append(meshInitializer.constructHalfEdgeMesh(newModel));
        ui->MainDisplay->updateBuffers(ui->MainDisplay->settings.meshes[0]);
        ui->MainDisplay->settings.modelLoaded = true;
    } else {
        qDebug() << "Model not loaded ";
        ui->MainDisplay->settings.modelLoaded = false;
    }
    ui->MeshGroupBox->setEnabled(ui->MainDisplay->settings.modelLoaded);
    ui->sharpnessSettings->setEnabled(ui->MainDisplay->settings.modelLoaded);

    ui->lcdNumber->display(0);
    ui->SubdivSteps->setValue(0);
    ui->MainDisplay->update();
}

void MainWindow::on_LoadOBJ_pressed() {
    QString filename = QFileDialog::getOpenFileName(this, "Import OBJ File", "../", tr("Obj Files (*.obj)"));
    importOBJ(filename);
}

void MainWindow::on_MeshPresetComboBox_currentTextChanged(const QString& meshName) {
    importOBJ(":/models/" + meshName + ".obj");
}

void MainWindow::on_SubdivSteps_valueChanged(int value) {
    ui->MainDisplay->settings.subDivValue = value;
    Subdivider* subdivider = new CatmullClarkSubdivider();
    for (int k = ui->MainDisplay->settings.meshes.size() - 1; k < value; k++) {
        ui->MainDisplay->settings.meshes.append(subdivider->subdivide(ui->MainDisplay->settings.meshes[k]));
    }
    ui->MainDisplay->updateBuffers(ui->MainDisplay->settings.meshes[value]);
    delete subdivider;
}

void MainWindow::on_sharpnessSliderValue_valueChanged(int intSharpness) {
    double decSharpness = ui->MainDisplay->settings.decSharpnessSelectedHE;
    double sharpness = intSharpness + decSharpness;

    ui->lcdNumber->display(sharpness);
    ui->MainDisplay->settings.intSharpnessSelectedHE = intSharpness;
    // Sharpness for current half-edge
    ui->MainDisplay->settings.selectedHE->sharpness = sharpness;
    // Sharpness for twin half-edge
    if (!ui->MainDisplay->settings.selectedHE->isBoundaryEdge()) {
        ui->MainDisplay->settings.selectedHE->twin->sharpness = sharpness;
    }
    update();
}

void MainWindow::timeout() {
    if (ui->MainDisplay->settings.edgeSlectionEnabled && ui->MainDisplay->settings.isEdgeSelected) {

        ui->sharpnessSliderValue->setEnabled(true);
        ui->decimalSharpnessSpinBox->setEnabled(true);
        double intSharpness = floor(ui->MainDisplay->settings.selectedHE->sharpness);
        double decSharpness = ui->MainDisplay->settings.selectedHE->sharpness - intSharpness;

        ui->MainDisplay->settings.intSharpnessSelectedHE = intSharpness;
        ui->MainDisplay->settings.decSharpnessSelectedHE = decSharpness;

        ui->lcdNumber->display(ui->MainDisplay->settings.selectedHE->sharpness);
        ui->decimalSharpnessSpinBox->setValue(decSharpness);
        ui->sharpnessSliderValue->setValue(ui->MainDisplay->settings.selectedHE->sharpness);
    } else {
        ui->sharpnessSliderValue->setEnabled(false);
        ui->decimalSharpnessSpinBox->setEnabled(false);
    }

    if (ui->SubdivSteps->value() > 0 ) {
        ui->sharpnessSliderValue->setEnabled(false);
        ui->decimalSharpnessSpinBox->setEnabled(false);
    } else if (ui->SubdivSteps->value() == 0 && ui->MainDisplay->settings.edgeSlectionEnabled
               && ui->MainDisplay->settings.isEdgeSelected) {
        ui->sharpnessSliderValue->setEnabled(true);
        ui->decimalSharpnessSpinBox->setEnabled(true);
    }
}

void MainWindow::on_decimalSharpnessSpinBox_valueChanged(double decSharpness) {
    double integerSharpness = ui->MainDisplay->settings.intSharpnessSelectedHE;
    double sharpness = integerSharpness + decSharpness;
    ui->lcdNumber->display(sharpness);
    ui->MainDisplay->settings.decSharpnessSelectedHE = decSharpness;

    // Sharpness for current half-edge
    ui->MainDisplay->settings.selectedHE->sharpness = sharpness;
    // Sharpness for twin half-edge
    if (!ui->MainDisplay->settings.selectedHE->isBoundaryEdge()) {
        ui->MainDisplay->settings.selectedHE->twin->sharpness = sharpness;
    }
    update();
}

void MainWindow::on_edgeSelectionRadioButton_toggled(bool checked) {
    ui->MainDisplay->settings.edgeSlectionEnabled = checked;
    ui->MainDisplay->updateBuffers(ui->MainDisplay->settings.meshes[ui->SubdivSteps->value()]);
    update();
}

