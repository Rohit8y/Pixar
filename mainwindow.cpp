#include "mainwindow.h"
#include "ui_mainwindow.h"

#include "initialization/meshinitializer.h"
#include "initialization/objfile.h"
#include "subdivision/catmullclarksubdivider.h"
#include "subdivision/subdivider.h"


MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent), ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    ui->MeshGroupBox->setEnabled(ui->MainDisplay->settings.modelLoaded);
}

MainWindow::~MainWindow()
{
    delete ui;

    meshes.clear();
    meshes.squeeze();
}

/**
 * @brief MainWindow::importOBJ Imports an obj file and adds the constructed
 * half-edge to the collection of meshes.
 * @param fileName Path of the .obj file.
 */
void MainWindow::importOBJ(const QString& fileName) {
  OBJFile newModel = OBJFile(fileName);
  meshes.clear();
  meshes.squeeze();

  if (newModel.loadedSuccessfully()) {
    MeshInitializer meshInitializer;
    meshes.append(meshInitializer.constructHalfEdgeMesh(newModel));
    ui->MainDisplay->updateBuffers(meshes[0]);
    ui->MainDisplay->settings.modelLoaded = true;
  } else {
    ui->MainDisplay->settings.modelLoaded = false;
  }

  ui->MeshGroupBox->setEnabled(ui->MainDisplay->settings.modelLoaded);
  ui->SubdivSteps->setValue(0);
  ui->MainDisplay->update();
}

void MainWindow::on_LoadOBJ_pressed() {
  QString filename = QFileDialog::getOpenFileName(
      this, "Import OBJ File", "../", tr("Obj Files (*.obj)"));
  importOBJ(filename);
}

void MainWindow::on_MeshPresetComboBox_currentTextChanged(
    const QString& meshName) {
  importOBJ(":/models/" + meshName + ".obj");
}
void MainWindow::on_SubdivSteps_valueChanged(int value)
{
    Subdivider* subdivider = new CatmullClarkSubdivider();
    for (int k = meshes.size() - 1; k < value; k++) {
      meshes.append(subdivider->subdivide(meshes[k]));
    }
    ui->MainDisplay->updateBuffers(meshes[value]);
    delete subdivider;
}

void MainWindow::on_HideMeshCheckBox_toggled(bool arg1)
{

}


void MainWindow::on_TessellationCheckBox_toggled(bool arg1)
{

}


