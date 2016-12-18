//* Authors: Dmitry Anisimov
//* mainwindow.hpp
//* With any questions please contact: danston@ymail.com

#ifndef MAINWINDOW_HPP
#define MAINWINDOW_HPP

// Local includes.
#include "qtviewer.hpp"

// Qt includes.
#include <QMainWindow>

namespace Ui {
    class MainWindow;
}

// warpit namespace.
namespace warpit {

    // Main window class.
    class MainWindow : public QMainWindow
    {
        Q_OBJECT

    public:
        // Constructor.
        explicit MainWindow(QWidget *parent = 0); // add "explicit" to the constructor declaration to prevent implicit conversions

        // Destructor.
        ~MainWindow();

    private slots:
        // Load routines.
        void on_actionLoadTexture_triggered();
        void on_actionLoadSource_triggered();
        void on_actionLoadTarget_triggered();
        void on_actionLoadMesh_triggered();

        // Save routines.
        void on_actionSaveSource_triggered();
        void on_actionSaveTarget_triggered();
        void on_actionSaveMesh_triggered();

        void on_showGridBox_clicked(const bool gridStatus); // draw grid in the viewer
        void on_showMeshBox_clicked(const bool meshStatus); // draw mesh in the viewer

        // Change max number of subdivision steps.
        void on_subdivisionStepsEdit_textChanged();

        // Apply mesh preprocessing.
        void on_preprocessMeshBox_clicked(const bool preprocessingStatus);

        // Morphing.
        void on_startMorphingBox_clicked(const bool morphingStatus);

        // Increase and decrease number of subdivision steps.
        void on_upStepsButton_released();
        void on_downStepsButton_released();

        // Update the statistics slot.
        void updateStatistics(const int numFaces, const double morphingTime);

        // Choose the default = the subdivision method to morph the image.
        void on_defaultMorphButton_clicked();

        // Choose the barycentric method to morph the image.
        void on_barycentricMorphButton_clicked();

    private:
        // Pointers and references.
        Ui::MainWindow *ui;      // create a pointer to the interface class
        QTViewer       *viewer;  // create a pointer to the viewer class
        Program        &program; // create a reference to the instance of the main program class

        // Default path to the project.
        QString defaultPath;

        // Set all options to default.
        void setDefaultState();

        // Set all options to the unmorphed state.
        void setUnmorphedState();

        // Set all options to the morphed state.
        void setMorphedStateStart();
        void setMorphedStateEnd();

        // Save mesh.
        void saveMesh(const QString &fileName, const std::vector<Vertex> &vertices, const std::vector<Halfedge> &halfedges, const bool morphed);
    };

} // namespace warpit

#endif // MAINWINDOW_HPP
