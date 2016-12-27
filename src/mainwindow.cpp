//* Author: Dmitry Anisimov (c) 2016-2017
//* With any questions please contact: danston@ymail.com

// Todo list:
// 1. Finish loading obj mesh.
// 2. Export materials when saving the mesh.
// 3. Split saving mesh from materials and move it to the mesh class.

// STL includes.
#include <fstream>

// Local includes.
#include "mainwindow.hpp"
#include "ui_mainwindow.h"

// Qt includes.
#include <QFileDialog>

// warpit namespace
namespace warpit {

    // Constructor.
    MainWindow::MainWindow(QWidget *parent) :
        QMainWindow(parent),
        ui(new Ui::MainWindow),
        program(Program::instance()),
        defaultPath(PROJECT_PATH)
    {
        ui->setupUi(this);                                               // set up ui
        ui->subdivisionStepsEdit->setValidator(new QIntValidator(0,25)); // limit input values in the subdivision box to [0,25]
        viewer = static_cast<QTViewer *>(ui->viewerWidget);              // cast the viewer

        // Update statistics.
        connect(viewer, SIGNAL(updateStatistics(const int, const double)), this, SLOT(updateStatistics(const int, const double)));

        // Use program without interface.
        // QString path = QString::fromStdString("/path_to/warpit/polygons/sources/cow.poly");
        // program.run(path, 100, 8, false);
    }

    // Update statistics.
    void MainWindow::updateStatistics(const int numFaces, const double morphingTime)
    {
        if(program.getProgramMode() != program.INITIAL) {

            QString text = QString::number(numFaces);
            ui->facesLabel->setText("Faces: " + text);
            text = QString::number(morphingTime);
            ui->morphingTimeLabel->setText("Time: " + text + " sec.");
        }
    }

    // Destructor.
    MainWindow::~MainWindow()
    {
        // Delete all interface stuff.
        delete ui;
    }

    // Draw grid in the viewer.
    void MainWindow::on_showGridBox_clicked(const bool gridStatus)
    {
        viewer->setGridStatus(gridStatus);
        viewer->update();
    }

    // Draw mesh.
    void MainWindow::on_showMeshBox_clicked(const bool meshStatus)
    {
        viewer->setMeshStatus(meshStatus);
        viewer->update();
    }

    // Apply vertex adjsutment.
    void MainWindow::on_preprocessMeshBox_clicked(const bool preprocessingStatus)
    {
        program.setPreprocessingStatus(preprocessingStatus);
        viewer->update();
    }

    // Choose the default = the subdivision method to morph the image.
    void MainWindow::on_defaultMorphButton_clicked()
    {
        program.toggleMorphingMethod(true);
        viewer->update();
    }

    // Choose the barycentric method to morph the image.
    void MainWindow::on_barycentricMorphButton_clicked()
    {
        program.toggleMorphingMethod(false);
        viewer->update();
    }

    // Increase subdivision steps.
    void MainWindow::on_upStepsButton_released()
    {
        int numSubSteps = program.getNumSubSteps();

        ++numSubSteps;
        if(numSubSteps > program.getMaxNumSubSteps()) return;

        program.setNumSubSteps(numSubSteps);

        QString text = QString::number(numSubSteps);

        ui->stepsResultLabel->setText(": [ " + text + " ]");

        program.subdivideUp();

        const int numFaces = program.halfedges(numSubSteps).size() / 3;
        const double morphingTime = program.getMorphingTime();
        updateStatistics(numFaces, morphingTime);

        viewer->update();
    }

    // Decrease subdivision steps.
    void MainWindow::on_downStepsButton_released()
    {
        int numSubSteps = program.getNumSubSteps();

        --numSubSteps;
        if(numSubSteps < 0) return;

        program.setNumSubSteps(numSubSteps);

        QString text = QString::number(numSubSteps);

        ui->stepsResultLabel->setText(": [ " + text + " ]");

        program.subdivideDown();

        const int numFaces = program.halfedges(numSubSteps).size() / 3;
        const double morphingTime = 0.0;
        updateStatistics(numFaces, morphingTime);

        viewer->update();
    }

    // Change max number of subdivision steps.
    void MainWindow::on_subdivisionStepsEdit_textChanged()
    {
        const QString text = ui->subdivisionStepsEdit->text();
        const int maxNumSubSteps = text.toInt();
        program.setMaxNumSubSteps(maxNumSubSteps);

        ui->maxStepsResultLabel->setText(": [ " + text + " ]");
    }

    // Morphing.
    void MainWindow::on_startMorphingBox_clicked(const bool morphingStatus)
    {
        if(morphingStatus) {
            setMorphedStateStart();
            program.morph();

            const int numFaces = program.halfedges(program.getNumSubSteps()).size() / 3;
            const double morphingTime = program.getMorphingTime();
            updateStatistics(numFaces, morphingTime);
        }
        else {
            setMorphedStateEnd();

            program.triMesh().clear();
            program.triangulate();

            const int numFaces = program.triMesh().numFaces();
            const double morphingTime = program.getMorphingTime();
            updateStatistics(numFaces, morphingTime);
        }

        viewer->update();
    }

    // Set all options to default.
    void MainWindow::setDefaultState()
    {
        program.clear();

        program.setProgramMode(program.INITIAL);
        ui->startMorphingBox->setEnabled(false);

        ui->morphingGroupBox->setEnabled(false);
        ui->defaultMorphButton->setChecked(true);
        ui->barycentricMorphButton->setChecked(false);
        ui->defaultMorphButton->setEnabled(false);
        ui->barycentricMorphButton->setEnabled(false);

        program.setPreprocessingStatus(false);
        ui->preprocessMeshBox->setChecked(false);
        ui->preprocessMeshBox->setEnabled(false);

        program.setMaxNumSubSteps(0);
        ui->maxStepsResultLabel->setText(": [ 0 ]");
        ui->subdivisionStepsEdit->setText("0");
        ui->subdivisionStepsEdit->setEnabled(false);

        program.setNumSubSteps(0);
        ui->stepsResultLabel->setText(": [ 0 ]");
        ui->downStepsButton->setEnabled(false);
        ui->upStepsButton->setEnabled(false);

        viewer->setMeshStatus(false);
        ui->showMeshBox->setChecked(false);
        ui->showMeshBox->setEnabled(false);

        viewer->setGridStatus(true);
        ui->showGridBox->setChecked(true);
        ui->showGridBox->setEnabled(true);

        ui->actionLoadTexture->setEnabled(true);

        ui->actionLoadSource->setEnabled(false);
        ui->actionLoadTarget->setEnabled(false);
        ui->actionLoadMesh->setEnabled(false);

        ui->actionSaveSource->setEnabled(false);
        ui->actionSaveTarget->setEnabled(false);
        ui->actionSaveMesh->setEnabled(false);
    }

    // Set all options to the unmorphed state.
    void MainWindow::setUnmorphedState()
    {
        program.setProgramMode(program.PREPROCESSING);

        ui->startMorphingBox->setEnabled(true);

        ui->morphingGroupBox->setEnabled(true);
        ui->defaultMorphButton->setEnabled(true);
        ui->barycentricMorphButton->setEnabled(true);

        ui->preprocessMeshBox->setEnabled(true);
        ui->subdivisionStepsEdit->setEnabled(true);
        ui->showMeshBox->setEnabled(true);

        ui->actionLoadSource->setEnabled(true);
        ui->actionSaveSource->setEnabled(true);
        ui->actionSaveMesh->setEnabled(true);
    }

    // Set all options to the morphed state - at the beginning.
    void MainWindow::setMorphedStateStart()
    {
        program.setProgramMode(program.MORPHING);

        ui->startMorphingBox->setText("Stop");

        ui->preprocessMeshBox->setEnabled(false);
        ui->subdivisionStepsEdit->setEnabled(false);

        ui->morphingGroupBox->setEnabled(false);

        ui->downStepsButton->setEnabled(true);
        ui->upStepsButton->setEnabled(true);

        ui->actionLoadTexture->setEnabled(false);
        ui->actionLoadSource->setEnabled(false);
        ui->actionLoadTarget->setEnabled(true);
        ui->actionLoadMesh->setEnabled(false);

        ui->actionSaveSource->setEnabled(false);
        ui->actionSaveTarget->setEnabled(true);
    }

    // Set all options to the morphed state - at the end.
    void MainWindow::setMorphedStateEnd()
    {
        program.setProgramMode(program.PREPROCESSING);

        ui->startMorphingBox->setText("Start");

        ui->morphingGroupBox->setEnabled(true);

        ui->preprocessMeshBox->setEnabled(true);
        ui->subdivisionStepsEdit->setEnabled(true);

        ui->downStepsButton->setEnabled(false);
        ui->upStepsButton->setEnabled(false);

        ui->stepsResultLabel->setText(": [ 0 ]");
        program.setNumSubSteps(0);

        ui->actionLoadTexture->setEnabled(true);
        ui->actionLoadSource->setEnabled(true);
        ui->actionLoadTarget->setEnabled(false);
        ui->actionLoadMesh->setEnabled(false);

        ui->actionSaveSource->setEnabled(true);
        ui->actionSaveTarget->setEnabled(false);

        program.setCurrMaxNumSubSteps(0);
    }

    // Load texture.
    void MainWindow::on_actionLoadTexture_triggered()
    {
        QString fileName = QFileDialog::getOpenFileName(this, tr("Load texture"), defaultPath, tr("Textures (*.png)"));

        // Open a texture with .png extension.
        if(!fileName.isNull() && fileName.contains(".png", Qt::CaseSensitive)) {
            // Load texture.
            setDefaultState();
            viewer->loadTexture(fileName.toStdString());
            setUnmorphedState();

            const int numFaces =  program.triMesh().numFaces();
            const double morphingTime = program.getMorphingTime();
            updateStatistics(numFaces, morphingTime);

            viewer->update();
        }
    }

    // Load source polygon.
    void MainWindow::on_actionLoadSource_triggered()
    {
        QString fileName = QFileDialog::getOpenFileName(this, tr("Load source"), defaultPath, tr("Polygons (*.poly)"));

        // Open a polygon with .poly extension.
        if(!fileName.isNull() && fileName.contains(".poly", Qt::CaseSensitive)) {
            // Load polygon.
            program.sourcePolygon().load(fileName.toStdString());
            program.cleanVirtualVertices();
            program.triangulate();

            const int numFaces =  program.triMesh().numFaces();
            const double morphingTime = program.getMorphingTime();
            updateStatistics(numFaces, morphingTime);

            viewer->update();
        }
    }

    // Load target polygon.
    void MainWindow::on_actionLoadTarget_triggered()
    {
        QString fileName = QFileDialog::getOpenFileName(this, tr("Load target"), defaultPath, tr("Polygons (*.poly)"));

        // Open a polygon with .poly extension.
        if(!fileName.isNull() && fileName.contains(".poly", Qt::CaseSensitive)) {
            // Load polygon.
            program.targetPolygon().load(fileName.toStdString());
            if(program.targetPolygon().numVertices() != program.sourcePolygon().numVertices()) {
                std::cout << "\nERROR: Number of vertices is not consistent with the source polygon!\n" << std::endl;
                exit(0);
            }

            program.subdivide();

            const int numFaces =  program.halfedges(program.getNumSubSteps()).size() / 3;
            const double morphingTime = program.getMorphingTime();
            updateStatistics(numFaces, morphingTime);

            viewer->update();
        }
    }

    // Load mesh with .obj extension.
    void MainWindow::on_actionLoadMesh_triggered()
    {
        QString fileName = QFileDialog::getOpenFileName(this, tr("Load mesh"), defaultPath, tr("Meshes (*.obj)"));

        // Open mesh with .obj extension.
        if(!fileName.isNull() && fileName.contains(".obj", Qt::CaseSensitive)) {
            // Load mesh.
            program.triMesh().load(fileName.toStdString());

            const int numFaces =  program.triMesh().numFaces();
            const double morphingTime = program.getMorphingTime();
            updateStatistics(numFaces, morphingTime);

            viewer->update();
        }
    }

    // Save the source polygon.
    void MainWindow::on_actionSaveSource_triggered()
    {
        QString fileName = QFileDialog::getSaveFileName(this, tr("Save source"), defaultPath, tr("Polygons (*.poly)"));

        // Save the source polygon with .poly extension.
        if(!fileName.isNull() && fileName.contains(".poly", Qt::CaseSensitive))
            program.sourcePolygon().save(fileName.toStdString());
    }

    // Save the target polygon.
    void MainWindow::on_actionSaveTarget_triggered()
    {
        QString fileName = QFileDialog::getSaveFileName(this, tr("Save target"), defaultPath, tr("Polygons (*.poly)"));

        // Save the target polygon with .poly extension.
        if(!fileName.isNull() && fileName.contains(".poly", Qt::CaseSensitive))
            program.targetPolygon().save(fileName.toStdString());
    }

    // Save the mesh with .obj extension.
    void MainWindow::on_actionSaveMesh_triggered()
    {
        QString fileName = QFileDialog::getSaveFileName(this, tr("Save mesh"), defaultPath, tr("Meshes (*.obj)"));

        // Save the mesh with .obj extension.
        if(!fileName.isNull() && fileName.contains(".obj", Qt::CaseSensitive)) {

            if(program.getProgramMode() == program.PREPROCESSING) {
                saveMesh(fileName, program.triMesh().vertices(), program.triMesh().halfedges(), false);
                return;
            }

            if(program.getProgramMode() == program.MORPHING) {
                const int numSubSteps = program.getNumSubSteps();
                saveMesh(fileName, program.vertices(numSubSteps), program.halfedges(numSubSteps), true);
                return;
            }
        }
    }

    // Save mesh with .obj extension.
    void MainWindow::saveMesh(const QString &fileName, const std::vector<Vertex> &vertices, const std::vector<Halfedge> &halfedges, const bool morphed)
    {
        QFileInfo info(fileName);

        std::string baseName = info.baseName().toStdString();

        std::string str = fileName.toStdString();

        if(str.empty()) return;

        std::ofstream saveFile(str.c_str(), std::ios_base::out);

        if(!saveFile) {
            std::cout << "\nError saving .obj file!" << std::endl;
            exit(EXIT_FAILURE);
        }

        const int numV  = vertices.size();
        const int numHE = halfedges.size();
        const int numF  = numHE / 3;

        if(numF == 0)
            saveFile << "You have not triangulated your polygon before saving it as .obj file!" << std::endl;
        else {
            saveFile << "# OBJ" << std::endl;
            saveFile << "# Vertices: " << numV << std::endl;
            saveFile << "mtllib " << baseName << ".mtl\n";

            if(morphed) {
                for(int i = 0; i < numV; ++i) {
                    saveFile << "v " << vertices[i].pm().x << " "
                                     << vertices[i].pm().y << " "
                                     << 0.0 << std::endl;

                    saveFile << "vt " << vertices[i].p().x << " "
                                      << vertices[i].p().y << "\n";
                }
            } else {
                Vector2d minB, maxB;
                program.sourcePolygon().boundingBox(minB, maxB);

                Vector2d tex;
                for(int i = 0; i < numV; ++i) {
                    tex = vertices[i].p();
                    program.texCoord(tex, minB, maxB);

                    saveFile << "v " << vertices[i].p().x << " "
                                     << vertices[i].p().y << " "
                                     << 0.0 << std::endl;

                    saveFile << "vt " << tex.x << " "
                                      << tex.y << "\n";
                }
            }

            saveFile << "# Faces: " << numF << std::endl;
            saveFile << "usemtl texture\n";

            int idxE = 0;
            for(int i = 0; i < numF; ++i) {
                saveFile << "f ";
                for(int j = 0; j < 3; ++j) {
                    const int id = halfedges[halfedges[idxE++].prev].dest;
                    saveFile << id + 1 << "/" << id + 1 << "/ ";
                }
                saveFile << std::endl;
            }

            saveFile << "# End of file";
        }
        saveFile.close();

        std::string path = info.path().toStdString() + "/" + baseName + ".mtl";

        saveFile.open(path.c_str());

        std::string map_Kd = "map_Kd " + baseName + ".png";

        saveFile << "newmtl texture\n";
        saveFile << "Ka 1 1 1\n";
        saveFile << "Kd 1 1 1\n";
        saveFile << "Ks 1 1 1\n";
        saveFile <<  map_Kd;

        saveFile.close();

        viewer->getTexture()->getImage().save(QString::fromStdString(info.path().toStdString() + "/" + baseName + ".png"));
    }

} // namespace warpit
