//* Author: Dmitry Anisimov (c) 2016-2017
//* With any questions please contact: danston@ymail.com

#ifndef QTVIEWER_HPP
#define QTVIEWER_HPP

// Local includes.
#include "program.hpp"

// Qt includes.
#include <QGLWidget>
#include <QPainter>
#include <QMouseEvent>

// warpit namespace
namespace warpit {

    // QT viewer.
    class QTViewer: public QGLWidget
    {
        Q_OBJECT

    public:
        // Constructor.
        explicit QTViewer(QWidget* parent = 0, const QGLWidget * shareWidget = 0, Qt::WindowFlags f = 0);

        // Destructor.
        ~QTViewer();

        // Show mesh.
        inline void setMeshStatus(const bool newMeshStatus) { meshStatus = newMeshStatus; }

        // Show grid.
        inline void setGridStatus(const bool newGridStatus) { gridStatus = newGridStatus; }

        // Load new texture.
        void loadTexture(const std::string &pathToTexture);

        // Get texture.
        inline const Texture* getTexture() const
        {
            return texture;
        }

    signals:
        // Update statistics.
        void updateStatistics(const int numFaces, const double morphingTime);

    protected:
        // Initialize GL.
        void initializeGL();

        // Paint event.
        void paintEvent(QPaintEvent* event);

        // Mouse events.
        void mousePressEvent(QMouseEvent* event);   // press event
        void mouseMoveEvent(QMouseEvent* event);    // move event
        void mouseReleaseEvent(QMouseEvent* event); // release event
        void wheelEvent(QWheelEvent* event);        // wheel event

    private:
        typedef QGLWidget superClass; // super class
        Program &program;             // create a reference to the main program class

        bool meshStatus;   // mesh status: true - show mesh in the viewer, false - do not show
        double zoomFactor; // viewer zoom factor
        bool gridStatus;   // gird status: true - show grid in the viewer, false - do not show

        QPoint buttonDown; // stores position of the point in the moment of press event
        QPoint buttonMove; // stores position of the point in the moment of move event

        Vertex *selectedVertex; // currently selected vertex

        Texture *texture; // stores a texture

        // Draw a regular grid in the viewer.
        void drawGrid(QPainter &painter, std::vector<std::vector<QPoint> > &grid) const;

        // Draw all the base/selected/action vertices.
        void drawVertices(QPainter &painter, Polygon &polygon);

        // Draw the polygon.
        void drawPolygon(QPainter &painter, Polygon &polygon);

        // Draw stuff before the morphing.
        void drawUnmorphed(QPainter &painter);
        void drawUnmorphedMesh(QPainter &painter);
        void drawUnmorphedTexture(QPainter &painter);

        // Draw stuff after the morphing.
        void drawMorphed(QPainter &painter);
        void drawMorphedMesh(QPainter &painter);
        void drawMorphedTexture(QPainter &painter);

        // Add a new vertex to the source polygon.
        void addVertex(const QPoint &buttonUp);

        // Delete a vertex from the source polygon.
        void deleteVertex(const QPoint &buttonUp);

        // Create a polygon from the texture.
        void createPolygonFromTexture();

        // Map world coordinates to screen coordinates.
        inline QPoint world2screen(const Vector2d &v, const double width, const double height) const
        {
            return QPoint(v.x * zoomFactor + 0.5 * width, -v.y * zoomFactor + 0.5 * height);
        }

        // Map screen coordinates to world coordinates.
        inline Vector2d screen2world(const QPoint &point, const double width, const double height) const
        {
            return Vector2d((point.x() - 0.5 * width) / zoomFactor, (-point.y() + 0.5 * height) / zoomFactor);
        }
    };

} // namespace warpit

#endif // QTVIEWER_HPP
