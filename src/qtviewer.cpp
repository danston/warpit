//* Authors: Dmitry Anisimov
//* qtviewer.cpp
//* With any questions please contact: danston@ymail.com

// Todo list:
// 1. Handle self-intersected polygons.
// 2. Add selection, moving, etc.

// Local includes.
#include "qtviewer.hpp"

// Qt includes.
#include <QPen>
#include <QApplication>

// warpit namespace.
namespace warpit {

    // Constructor.
    QTViewer::QTViewer(QWidget *parent, const QGLWidget *shareWidget, Qt::WindowFlags f) :
        superClass(parent, shareWidget, f),
        program(Program::instance()),
        meshStatus(false),
        zoomFactor(1000.0),
        gridStatus(true),
        selectedVertex(NULL),
        texture(NULL)
    {
        setAutoFillBackground(false);
        setMouseTracking(true);
    }

    // Destructor.
    QTViewer::~QTViewer()
    {
        delete texture;
    }

    // Initialize GL.
    void QTViewer::initializeGL()
    {
        CHECK_GL_VALID_STATE();

        glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
        glEnable(GL_COLOR_MATERIAL);
        glEnable(GL_TEXTURE_2D);
        glShadeModel(GL_SMOOTH);
        glClearColor(1, 1, 1, 1);

        CHECK_GL_VALID_STATE();
    }

    // Load new texture.
    void QTViewer::loadTexture(const std::string &pathToTexture)
    {
        delete texture;

        texture = new Texture(pathToTexture);
        texture->setTextureType(GL_TEXTURE_2D);

        createPolygonFromTexture();

        program.triangulate();
    }

    // Create a polygon from the texture.
    void QTViewer::createPolygonFromTexture()
    {
        const double texWidth  = texture->width();
        const double texHeight = texture->height();

        program.sourcePolygon() = Polygon(texWidth, texHeight);
        program.sourcePolygon().fitToSize();
        program.sourcePolygon().setBarycenter();

        program.targetPolygon() = program.sourcePolygon();
    }

    // Paint event.
    void QTViewer::paintEvent(QPaintEvent* event)
    {
        // Set defaults.
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        QPainter painter;
        painter.begin(this);

        // Draw the default white background.
        painter.setRenderHint(QPainter::Antialiasing);
        painter.fillRect(event->rect(), QBrush(Qt::white));

        // Draw a regular grid.
        if(gridStatus) {
            program.createRegularGrid(width(), height());
            drawGrid(painter, program.gridPoints());
        }

        // Draw texture, mesh, polygon, and vertices.
        switch(program.getProgramMode())
        {
            case 0:
                break;

            case 1:
                drawUnmorphed(painter); // unmorphed image
                break;

            case 2:
                drawMorphed(painter); // morphed image
                break;
        }

        painter.end();
    }

    // Draw stuff before the morphing.
    void QTViewer::drawUnmorphed(QPainter &painter)
    {
        // Draw the unmorphed texture.
        drawUnmorphedTexture(painter);

        // Draw the source polygon.
        drawPolygon(painter, program.sourcePolygon());

        // Draw the mesh.
        if(meshStatus == true) drawUnmorphedMesh(painter);

        // Draw all vertices.
        drawVertices(painter, program.sourcePolygon());
    }

    // Draw stuff after the morphing.
    void QTViewer::drawMorphed(QPainter &painter)
    {
        // Draw the morphed texture.
        drawMorphedTexture(painter);

        // Draw the target polygon.
        drawPolygon(painter, program.targetPolygon());

        // Draw the morphed mesh.
        if(meshStatus == true) drawMorphedMesh(painter);

        // Draw all vertices.
        drawVertices(painter, program.targetPolygon());
    }

    // Draw the unmorphed texture.
    void QTViewer::drawUnmorphedTexture(QPainter &painter)
    {
        assert(!program.triMesh().empty());

        // Compute bounding box of the source polygon.
        Vector2d minB, maxB;
        program.sourcePolygon().boundingBox(minB, maxB);

        // Set defaults.
        const int numHE = program.triMesh().numHalfedges();
        const int numF  = numHE / 3;

        const std::vector<Vertex>   &vertices  = program.triMesh().vertices();
        const std::vector<Halfedge> &halfedges = program.triMesh().halfedges();

        const double w = width();
        const double h = height();

        // Draw texture.
        painter.beginNativePainting();
        texture->bind();
        glBegin(GL_TRIANGLES);

        int idxE = 0;
        for(int i = 0; i < numF; ++i) {
            for(int j = 0; j < 3; ++j) {
                const int id = halfedges[halfedges[idxE++].prev].dest;
                const QPoint pos = world2screen(vertices[id].p(), w, h);

                Vector2d tex = vertices[id].p();
                program.texCoord(tex, minB, maxB);

                glTexCoord2d(tex.x, tex.y);
                glVertex2d(pos.x(), pos.y());
            }
        }

        glEnd();
        texture->unbind();
        painter.endNativePainting();
    }

    // Draw the morphed texture.
    void QTViewer::drawMorphedTexture(QPainter &painter)
    {
        const int numSubSteps = program.getNumSubSteps();

        const std::vector<Vertex>   &vertices  = program.vertices(numSubSteps);
        const std::vector<Halfedge> &halfedges = program.halfedges(numSubSteps);

        assert(vertices.size() != 0);
        assert(halfedges.size() != 0);

        const int numHE = halfedges.size();
        const int numF  = numHE / 3;

        const double w = width();
        const double h = height();

        painter.beginNativePainting();
        texture->bind();
        glBegin(GL_TRIANGLES);

        int idxE = 0;
        for(int i = 0; i < numF; ++i) {
            for(int j = 0; j < 3; ++j) {
                const int id = halfedges[halfedges[idxE++].prev].dest;
                const QPoint pos = world2screen(vertices[id].pm(), w, h);

                glTexCoord2d(vertices[id].p().x, vertices[id].p().y);
                glVertex2d(pos.x(), pos.y());
            }
        }

        glEnd();
        texture->unbind();
        painter.endNativePainting();
    }

    // Draw a regular grid in the viewer.
    void QTViewer::drawGrid(QPainter &painter, std::vector< std::vector<QPoint> > &grid) const
    {
        const QColor color("#E9E9E9");

        painter.setPen(QPen(color, 1.0));

        const int gridSize = grid.size();

        QPolygonF square(4);

        for(int i = 0; i < gridSize - 1; ++i) {
            for(int j = 0; j < gridSize - 1; ++j) {
                square.replace(0, QPointF(grid[i][j]));
                square.replace(1, QPointF(grid[i][j+1]));
                square.replace(2, QPointF(grid[i+1][j+1]));
                square.replace(3, QPointF(grid[i+1][j]));
                painter.drawPolygon(square);

                if((i + j) % 2 == 0) {
                    QPainterPath path;
                    path.addPolygon(square);
                    painter.fillPath(path,QBrush(color));
                }
            }
        }
    }

    // The drawing routine for a polygon.
    void QTViewer::drawPolygon(QPainter &painter, Polygon &polygon)
    {
        assert(!polygon.empty());

        const int numV = polygon.numVertices();
        const std::vector<Vertex> &vertices = polygon.vertices();

        const double w = width();
        const double h = height();

        QColor color("#7D7D7D");

        painter.setPen(QPen(color, 3.0));

        QPoint pos = world2screen(vertices[0].p(), w, h);

        for(int i = 1; i < numV; ++i) {
            painter.drawLine(pos, world2screen(vertices[i].p(), w, h));
            pos = world2screen(vertices[i].p(), w, h);
        }

        painter.drawLine(pos, world2screen(vertices[0].p(), w, h));
    }

    // The drawing routine for the unmorphed mesh.
    void QTViewer::drawUnmorphedMesh(QPainter &painter)
    {
        assert(!program.triMesh().empty());

        const int numF = program.triMesh().numFaces();

        const std::vector<Vertex> &vertices = program.triMesh().vertices();
        const std::vector<Face>   &faces    = program.triMesh().faces();

        const double w = width();
        const double h = height();

        const QColor color("#7D7D7D");

        painter.setPen(QPen(color, 2.0));

        painter.beginNativePainting();

        glColor3f(0.0f, 0.0f, 0.0f);
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

        glBegin(GL_TRIANGLES);

        for(int i = 0; i < numF; ++i) {
            for(int j = 0; j < 3; ++j) {
                const QPoint pos = world2screen(vertices[faces[i].v[j]].p(), w, h);
                glVertex2d(pos.x(), pos.y());
            }
        }

        glEnd();

        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        glColor3f(1.0f, 1.0f, 1.0f);

        painter.endNativePainting();
    }

    // The drawing routine for the morphed mesh.
    void QTViewer::drawMorphedMesh(QPainter &painter)
    {
        const int numSubSteps = program.getNumSubSteps();

        const std::vector<Vertex>   &vertices  = program.vertices(numSubSteps);
        const std::vector<Halfedge> &halfedges = program.halfedges(numSubSteps);

        assert(vertices.size() != 0);
        assert(halfedges.size() != 0);

        const int numHE = halfedges.size();
        const int numF  = numHE / 3;

        const double w = width();
        const double h = height();

        painter.beginNativePainting();

        glColor3f(0.0f, 0.0f, 0.0f);
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

        glBegin(GL_TRIANGLES);

        int idxE = 0;
        for(int i = 0; i < numF; ++i) {
            for(int j = 0; j < 3; ++j) {
                const int id = halfedges[halfedges[idxE++].prev].dest;
                const QPoint pos = world2screen(vertices[id].pm(), w, h);
                glVertex2d(pos.x(), pos.y());
            }
        }

        glEnd();

        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        glColor3f(1.0f, 1.0f, 1.0f);

        painter.endNativePainting();
    }

    // Draw polygon, virtual, and selected vertices.
    void QTViewer::drawVertices(QPainter &painter, Polygon &polygon)
    {
        // Set defaults.
        const int numV = polygon.numVertices();
        const std::vector<Vertex> &polyVertices = polygon.vertices();

        const double w = width();
        const double h = height();

        const QPoint vRadius(4.0, 4.0);

        QColor color("#C10026");

        painter.setPen(QPen(Qt::black, 1.0));

        // Draw vertices of the polygon. In red.
        painter.setBrush(color);
        for(int i = 0; i < numV; ++i) {
            const QPoint pos = world2screen(polyVertices[i].p(), w, h);
            painter.drawEllipse(QRect(pos - vRadius, pos + vRadius));
        }

        // Draw virtual vertices. In blue.
        if(program.getProgramMode() == program.PREPROCESSING && program.numVirtVertices() > 0) {

            const int numVV = program.numVirtVertices();
            const std::vector<Vertex> &virtualVertices = program.virtualVertices();

            color.setNamedColor("#0098C1");
            painter.setBrush(color);

            for(int i = 0; i < numVV; ++i) {
                const QPoint pos = world2screen(virtualVertices[i].p(), w, h);
                painter.drawEllipse(QRect(pos - vRadius, pos + vRadius));
            }
        }

        // Draw currently selected vertex - polygon or virtual. In green.
        if(selectedVertex != NULL) {
            color.setNamedColor("#00C114");
            painter.setBrush(color);

            const QPoint pos = world2screen(selectedVertex->p(), w, h);
            painter.drawEllipse(QRect(pos - vRadius, pos + vRadius));
        }
    }

    // The mouse press event.
    void QTViewer::mousePressEvent(QMouseEvent* event)
    {
        if(program.getProgramMode() == program.INITIAL) return;

        const int numV = program.sourcePolygon().numVertices();

        const double w = width();
        const double h = height();

        QPoint vertexTolerance(5.0, 5.0);

        buttonDown = event->pos();

        // Before morphing.
        if(program.getProgramMode() == program.PREPROCESSING) {
            std::vector<Vertex> &vertices = program.sourcePolygon().vertices();

            for(int i = 0; i < numV; ++i) {
                QPoint diff = buttonDown - world2screen(vertices[i].p(), w, h);
                if((std::fabs(diff.x()) < vertexTolerance.x()) && (std::fabs(diff.y()) < vertexTolerance.y())) {

                    buttonDown = world2screen(vertices[i].p(), w, h);
                    selectedVertex = &vertices[i];

                    update();
                    return;
                }
            }
        }

        // After morphing.
        if(program.getProgramMode() == program.MORPHING) {
            std::vector<Vertex> &vertices = program.targetPolygon().vertices();

            for(int i = 0; i < numV; ++i) {
                QPoint diff = buttonDown - world2screen(vertices[i].p(), w, h);
                if((std::fabs(diff.x()) < vertexTolerance.x()) && (std::fabs(diff.y()) < vertexTolerance.y())) {

                    buttonDown = world2screen(vertices[i].p(), w, h);
                    selectedVertex = &vertices[i];

                    update();
                    return;
                }
            }
            update();
            return;
        }

        // Virtual vertices.
        const int numVV = program.numVirtVertices();
        std::vector<Vertex> &virtualVertices = program.virtualVertices();

        for(int i = 0; i < numVV; ++i) {
            QPoint diff = buttonDown - world2screen(virtualVertices[i].p(), w, h);
            if((std::fabs(diff.x()) < vertexTolerance.x()) && (std::fabs(diff.y()) < vertexTolerance.y())) {

                buttonDown = world2screen(virtualVertices[i].p(), w, h);
                selectedVertex = &virtualVertices[i];

                update();
                break;
            }
        }
    }

    // The mouse move event.
    void QTViewer::mouseMoveEvent(QMouseEvent* event)
    {
        if(program.getProgramMode() == program.INITIAL) return;

        const double w = width();
        const double h = height();

        buttonMove = event->pos();

        // Move the vertex.
        if(event->buttons() == Qt::LeftButton && selectedVertex != NULL && !QApplication::keyboardModifiers().testFlag(Qt::ShiftModifier)) {

            if(buttonMove.x() < 0.0)     buttonMove.setX(0.0);
            if(buttonMove.x() > w - 1.0) buttonMove.setX(w - 1.0);
            if(buttonMove.y() < 0.0)     buttonMove.setY(0.0);
            if(buttonMove.y() > h - 1.0) buttonMove.setY(h - 1.0);

            selectedVertex->p() = screen2world(buttonMove, w, h);

            if(program.getProgramMode() == program.PREPROCESSING) program.triangulate();

            // Morphing.
            if(program.getProgramMode() == program.MORPHING) {

                program.subdivide();

                const int      numFaces   = program.halfedges(program.getNumSubSteps()).size() / 3;
                const double morphingTime = program.getMorphingTime();

                emit updateStatistics(numFaces, morphingTime);
            }
        }

        update();
    }

    // The mouse release event.
    void QTViewer::mouseReleaseEvent(QMouseEvent* event)
    {
        if(program.getProgramMode() == program.INITIAL) return;

        if(program.getProgramMode() == program.MORPHING) {
            selectedVertex = NULL;
            update();
            return;
        }

        // Add, move, or delete a vertex.
        if(selectedVertex != NULL) {
            if(event->button() == Qt::RightButton) { // Delete
                deleteVertex(event->pos());
            }
            else if(event->button() == Qt::LeftButton) { // Move

                Polygon &polygon = program.sourcePolygon();
                if(!polygon.contains(*selectedVertex))
                    polygon.projectVertexOnClosestEdge(*selectedVertex);
            }

            selectedVertex = NULL;
        } else if(event->button() == Qt::RightButton) addVertex(event->pos()); // Add

        // Remove all virtual vertices that lie outside the source polygon.
        program.cleanVirtualVertices();

        // Retriangulate the polygon.
        program.triangulate();

        // Update statistics.
        const int numFaces = program.triMesh().numFaces();
        const double morphingTime = program.getMorphingTime();

        emit updateStatistics(numFaces, morphingTime);

        // Udate viewer.
        update();
    }

    // Delete a vertex.
    void QTViewer::deleteVertex(const QPoint &buttonUp)
    {
        QPoint vertexTolerance(5.0, 5.0);

        std::vector<Vertex> &virtualVertices = program.virtualVertices();
        std::vector<Vertex> &vertices        = program.sourcePolygon().vertices();

        QPoint diff = buttonDown - buttonUp;
        if((std::fabs(diff.x()) < vertexTolerance.x()) && (std::fabs(diff.y()) < vertexTolerance.y())) {
            std::vector<Vertex>::iterator it = std::find(vertices.begin(), vertices.end(), *selectedVertex);
            if(it != vertices.end()) {
                if(vertices.size() > 4) vertices.erase(it);
            }
            else {
                it = std::find(virtualVertices.begin(), virtualVertices.end(), *selectedVertex);
                if(it != virtualVertices.end()) virtualVertices.erase(it);
                else assert(false);
            }
        }
    }

    // Add a new vertex.
    void QTViewer::addVertex(const QPoint &buttonUp)
    {
        const double w = width();
        const double h = height();

        const int numV = program.sourcePolygon().numVertices();
        std::vector<Vertex> &vertices = program.sourcePolygon().vertices();

        Vector2d newPosition(screen2world(buttonUp, w, h));

        // Add a new virtual vertex.
        if(QApplication::keyboardModifiers().testFlag(Qt::ShiftModifier)) {
            Vertex newVertex;
            newVertex.p() = newPosition;

            if(program.sourcePolygon().findEdgeIndex(newVertex) != -1)
               program.sourcePolygon().projectVertexOnClosestEdge(newVertex);

            program.virtualVertices().push_back(newVertex);

        } else if(numV < 2) {
            Vertex newVertex;
            newVertex.p() = newPosition;
            vertices.push_back(newVertex);
        }
        else {
            // Find the closest end vertex.
            double minDistance = (vertices[0].p() - newPosition).length();

            std::vector<Vertex>::iterator closest = vertices.begin();

            // Find the closest edge.
            for(std::vector<Vertex>::iterator iter = vertices.begin(); iter != vertices.end(); ++iter) {

                std::vector<Vertex>::iterator iterPlus = iter + 1;
                if(iterPlus == vertices.end()) iterPlus = vertices.begin();

                // Distance to line through an edge.
                double distance = std::fabs(((*iterPlus).p() - newPosition) % ((*iter).p() - newPosition)) / ((*iterPlus).p() - (*iter).p()).length();

                // If we are not "over" the edge, take the distance to the end vertex.
                Vector2d tmp((*iterPlus).p().x + (*iterPlus).p().y - (*iter).p().y, (*iterPlus).p().y - (*iterPlus).p().x + (*iter).p().x);

                if(((*iterPlus).p() - newPosition) % (tmp - newPosition) > 0.0) distance = ((*iterPlus).p() - newPosition).length();

                tmp = Vector2d((*iter).p().x + (*iterPlus).p().y - (*iter).p().y, (*iter).p().y - (*iterPlus).p().x + (*iter).p().x);

                if(((*iter).p() - newPosition) % (tmp - newPosition) < 0.0) {
                    distance = ((*iter).p() - newPosition).length();
                    iterPlus = iter;
                }

                if(distance < minDistance) {
                    minDistance = distance;
                    closest = iterPlus;
                }
            }

            Vertex newVertex;
            newVertex.p() = newPosition;
            vertices.insert(closest, newVertex);
        }
    }

    // The mouse wheel event.
    void QTViewer::wheelEvent(QWheelEvent* event)
    {
        if(program.getProgramMode() == program.INITIAL) return;

        const double speedFactor = 40.0;
        const double numDegrees  = -event->delta() / 8.0;
        const double numSteps    = numDegrees / 15.0;
        const double factor      = 1.125 * numSteps;

        if(event->delta() > 0.0)
            zoomFactor += (factor + speedFactor);
        else
            zoomFactor -= (factor + speedFactor);

        if(zoomFactor < 0.0) zoomFactor = 0.0;

        update();
    }

} // namespace warpit
