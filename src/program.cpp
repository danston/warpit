//* Author: Dmitry Anisimov (c) 2016-2017
//* With any questions please contact: danston@ymail.com

// Todo list:
// 1. Change x and y to x() and y() in the Vector2D class.
// 2. Optimize space consumption.

// STL includes.
#include <time.h>

// Local includes.
#include "program.hpp"

// Qt includes.
#include <QFileInfo>

// warpit namespace
namespace warpit {

    // The default constructor.
    Program::Program() :
        _sourcePolygon(),
        _targetPolygon(),
        _triMesh(),
        _vertices(),
        _halfedges(),
        _virtualVertices(),
        _gridPoints(),
        numSubSteps(0),
        maxNumSubSteps(0),
        currMaxSteps(0),
        programMode(INITIAL),
        preprocessingStatus(false),
        defaultMorphingMethod(true),
        eps(1.0e-14),
        morphingTime(0.0),
        initMesh()
    {

    }

    // Remove all virtual vertices that lie outside the source polygon.
    void Program::cleanVirtualVertices()
    {
        const int numVV = numVirtVertices();

        for(int i = numVV - 1; i >= 0; --i)
            if(_sourcePolygon.findEdgeIndex(_virtualVertices[i]) == -1 && !_sourcePolygon.contains(_virtualVertices[i]))
                _virtualVertices.erase(_virtualVertices.begin() + i);
    }

    // Triangulate the source polygon.
    void Program::triangulate()
    {
        const int numVV = numVirtVertices();

        Triangulator triangulator(_sourcePolygon);
        if(numVV > 0) {
             triangulator.generate(initMesh);

             Polygon tmpPolygon = _sourcePolygon;
             std::vector<Vertex> tmpVertices;

             for(int i = 0; i < numVV; ++i) {
                 const int boundaryIndex = tmpPolygon.findEdgeIndex(_virtualVertices[i]);
                 if(boundaryIndex >= 0) tmpPolygon.insert(_virtualVertices[i], boundaryIndex);
                 else tmpVertices.push_back(_virtualVertices[i]);
             }

             triangulator.setPolygon(tmpPolygon);
             triangulator.setVirtualVertices(tmpVertices);

             triangulator.generate(_triMesh);
        }
        else triangulator.generate(_triMesh);
    }

    // Initiate the mesh.
    void Program::initiateMesh()
    {
        if(numVirtVertices() > 0)
            fixVirtualCoordinates();

        setMeshInitialData();
    }

    // Preprocess the mesh.
    void Program::preprocessMesh()
    {
        if(preprocessingStatus == false) return;

        _triMesh.preprocess();
    }

    // Subdivide the mesh.
    void Program::subdivide()
    {
        setTargetInitialData();

        if(defaultMorphingMethod) morphWithSubdivisionWeights();
        else morphWithBarycentricCoordinates();

        currMaxSteps = numSubSteps;

        // compareMorphingMethods();
    }

    // Set initial data for a target mesh.
    void Program::setTargetInitialData()
    {
        const int tmp = numSubSteps;
        numSubSteps = 0;
        morphWithBarycentricCoordinates(); // no need for both methods, see below

        // use instead:
        // const int numV = _vertices[numSubSteps].size();
        // for(int i = 0; i < numV; ++i) _vertices[numSubSteps][i].pm() = _vertices[numSubSteps][i].p();
        //

        numSubSteps = tmp;
    }

    // Subdivide the mesh up.
    void Program::subdivideUp()
    {
        if(numSubSteps > currMaxSteps) {
            if(defaultMorphingMethod) {
                const clock_t startClock = clock();
                _triMesh.subdivideMorphed(_halfedges[currMaxSteps], _halfedges[numSubSteps],
                                          _vertices[currMaxSteps] , _vertices[numSubSteps]);
                const clock_t stopClock  = clock();
                morphingTime = double(stopClock - startClock) / CLOCKS_PER_SEC;

                // compareMorphingMethods();

                currMaxSteps = numSubSteps;
                return;

            } else {
                morphWithBarycentricCoordinates();
                currMaxSteps = numSubSteps;
                return;
            }
        }

        morphingTime = 0.0;
    }

    // Subdivide the mesh down.
    void Program::subdivideDown()
    {
        if(!defaultMorphingMethod) morphWithBarycentricCoordinates();
    }

    // Morph the mesh.
    void Program::morph()
    {
        _targetPolygon = _sourcePolygon;

        _halfedges.clear(); // need for barycentric method

        // TIME START
        // const clock_t startClock = clock();

        _halfedges.resize(maxNumSubSteps + 1); // no need for barycentric method
        _vertices.resize(maxNumSubSteps + 1);

        initiateMesh();

        preprocessMesh();

        _halfedges[0] = _triMesh.halfedges(); // no need for barycentric method
        _vertices[0]  = _triMesh.vertices();

        for(int i = 1; i < maxNumSubSteps + 1; ++i) {
            _triMesh.subdivide();
            _halfedges[i] = _triMesh.halfedges(); // no need for barycentric method
            _vertices[i]  = _triMesh.vertices();
        }

        // const clock_t stopClock  = clock();
        // morphingTime = double(stopClock - startClock) / CLOCKS_PER_SEC;
        // std::cout << morphingTime << std::endl;
        // TIME END

        // checkBarycentricProperties();

        _triMesh.clear();

        updateTextureCoordinates();
    }

    // Update texture coordinates for each subdivision step.
    void Program::updateTextureCoordinates()
    {
        Vector2d minB, maxB;
        _sourcePolygon.boundingBox(minB, maxB);

        const int size = _vertices.size();

        for(int i = 0; i < size; ++i) {
            const int numV = _vertices[i].size();
            for(int j = 0; j < numV; ++j) texCoord(_vertices[i][j].p(), minB, maxB);
        }
    }

    // Set initial data to the mesh.
    void Program::setMeshInitialData()
    {
        const int numVV = numVirtVertices();

        if(numVV == 0) {
            for(int i = 0; i < _targetPolygon.numVertices(); ++i)
                _triMesh.vertices()[i].pm() = _targetPolygon.vertices()[i].p();

            _triMesh.initializeBarycentricCoordinates(); // no need for subdivision method

            return;
        }

        const int numV = _triMesh.numVertices();
        const int numC = _sourcePolygon.numVertices();

        for(int i = 0; i < numV; ++i) {
            const int vInd = _sourcePolygon.findClosestVertex(_triMesh.vertices()[i]);
            if(vInd != -1) {
                _triMesh.vertices()[i].b().resize(numC, 0.0);
                _triMesh.vertices()[i].b()[vInd] = 1.0;

                Vector2d morphed;
                for(int j = 0; j < numC; ++j) morphed += _targetPolygon.vertices()[j].p() * _triMesh.vertices()[i].b()[j];
                _triMesh.vertices()[i].pm() = morphed;
            }
        }

        for(int i = 0; i < numVV; ++i) {
            int vertInd = -1;
            for(int j = 0; j < _triMesh.numVertices(); ++j)
                if(_triMesh.vertices()[j].p() == _virtualVertices[i].p())
                { vertInd = j; break; }

            assert(vertInd != -1);

            Vector2d morphed;
            for(int j = 0; j < numC; ++j) morphed += _targetPolygon.vertices()[j].p() * _virtualVertices[i].b()[j];
            _triMesh.vertices()[vertInd].pm() = morphed;
            _triMesh.vertices()[vertInd].b()  = _virtualVertices[i].b();
        }
    }

    // Morph with subdivision weights.
    void Program::morphWithSubdivisionWeights()
    {
        const clock_t startClock = clock();

        for(int i = 0; i < numSubSteps; ++i)
            _triMesh.subdivideMorphed(_halfedges[i], _halfedges[i+1], _vertices[i], _vertices[i+1]);

        const clock_t stopClock  = clock();
        morphingTime = double(stopClock - startClock) / CLOCKS_PER_SEC;
    }

    // Morph with barycentric coordinates.
    void Program::morphWithBarycentricCoordinates()
    {
        const int numV = _vertices[numSubSteps].size();
        const int numC = _sourcePolygon.numVertices();

        const clock_t startClock = clock();

        for(int i = 0; i < numV; ++i) {
            Vector2d morphed;
            for(int j = 0; j < numC; ++j) morphed += _targetPolygon.vertices()[j].p() * _vertices[numSubSteps][i].b()[j];
            _vertices[numSubSteps][i].pm() = morphed;
        }

        const clock_t stopClock  = clock();
        morphingTime = double(stopClock - startClock) / CLOCKS_PER_SEC;
    }

    // Fix barycentric coordinates for virtual vertices.
    void Program::fixVirtualCoordinates()
    {
        const int numVV = numVirtVertices();

        assert(numVV > 0);

        initMesh.initializeBarycentricCoordinates();

        const int numF  = initMesh.numFaces();
        const int numC  = initMesh.numVertices();

        const std::vector<Vertex> &vertices = initMesh.vertices();
        const std::vector<Face>   &faces    = initMesh.faces();

        for(int i = 0; i < numVV; ++i) {

            int count = 0;
            for(int j = 0; j < numF; ++j) {

                std::vector<double> lambda;
                computeTriangleCoordinates(vertices[faces[j].v[0]], vertices[faces[j].v[1]], vertices[faces[j].v[2]], i, lambda);

                if(lambda[0] > -eps && lambda[1] > -eps && lambda[2] > -eps) {
                    _virtualVertices[i].b().resize(numC);

                    for(int k = 0; k < numC; ++k) {
                        _virtualVertices[i].b()[k] = std::fabs(lambda[0]) * vertices[faces[j].v[0]].b()[k] +
                                                     std::fabs(lambda[1]) * vertices[faces[j].v[1]].b()[k] +
                                                     std::fabs(lambda[2]) * vertices[faces[j].v[2]].b()[k] ;
                    } break;
                } else ++count;
            }
            assert(count != numF);
        }
    }

    // Compute triangle coordinates.
    void Program::computeTriangleCoordinates(const Vertex &v1, const Vertex &v2, const Vertex &v3, const int vertInd, std::vector<double> &coordinates) const
    {
        coordinates.resize(3);

        const double area_second = 0.5 * crossProduct(v2, v3, _virtualVertices[vertInd]);
        const double area_third  = 0.5 * crossProduct(v3, v1, _virtualVertices[vertInd]);

        const double total_area          = 0.5 * crossProduct(v1, v2, v3);
        const double inverted_total_area = 1.0 / total_area;

        coordinates[0] = area_second * inverted_total_area;
        coordinates[1] = area_third  * inverted_total_area;
        coordinates[2] = 1.0 - coordinates[0] - coordinates[1];
    }

    // Compute cross product.
    double Program::crossProduct(const Vertex &v1, const Vertex &v2, const Vertex &v3) const
    {
        const double x1 = v1.p().x; const double y1 = v1.p().y;
        const double x2 = v2.p().x; const double y2 = v2.p().y;
        const double x3 = v3.p().x; const double y3 = v3.p().y;

        return x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2);
    }

    // Create a regular grid.
    void Program::createRegularGrid(const double width, const double height, const int gridSize)
    {
        _gridPoints.clear();
        _gridPoints.resize(gridSize + 1);

        const double offset = 2.0;

        const double dx = (width  - 1.7 * offset) / gridSize;
        const double dy = (height - 1.7 * offset) / gridSize;

        for(int i = 0; i <= gridSize; ++i)
            for(int j = 0; j <= gridSize; ++j)
                _gridPoints[i].push_back(QPoint(i * dx + offset, j * dy + offset));
    }

    // Compare barycentric and subdivision morphing methods.
    // NOTE: Works only when the default morphing method is the subdivision one!
    void Program::compareMorphingMethods()
    {
        const int numV = _vertices[numSubSteps].size();
        const int numC = _sourcePolygon.numVertices();

        for(int i = 0; i < numV; ++i) {
            const Vector2d vms = _vertices[numSubSteps][i].pm();
            Vector2d vmb;
            for(int j = 0; j < numC; ++j) vmb += _targetPolygon.vertices()[j].p() * _vertices[numSubSteps][i].b()[j];

            if(vms != vmb) std::cout << "\nERROR: Points do not coincide with error = " << vms - vmb << "; " << std::endl;
        }
        std::cout << "\n";
    }

    // Check properties of barycentric coordinates.
    void Program::checkBarycentricProperties()
    {
        const bool partitionOfUnity = _triMesh.partitionOfUnity();
        const bool linearPrecision  = _triMesh.linearPrecision(_sourcePolygon.vertices());
        const bool boundedness      = _triMesh.boundedness();

        if(!partitionOfUnity || !linearPrecision || !boundedness)
            std::cout << "\nERROR: Barycentric properties: FAILED!\n" << std::endl;
    }

    // Run program without interface.
    void Program::run(QString &path, const int timesToSubdivide, const int maxSteps, bool defaultMethod)
    {
        QFileInfo info(path);
        const std::string name = info.baseName().toStdString();

        /// Source:

        // Load source polygon.
        sourcePolygon().load(path.toStdString());
        triangulate();

        /// Target:

        // Load target polygon.
        std::string str = "/sources/" + name + ".poly";
        path.remove(QString::fromStdString(str));

        str = "/targets/" + name + ".poly";
        path.append(QString::fromStdString(str));

        targetPolygon().load(path.toStdString());

        /// Options:

        // File name to save time.
        str = "polygons/targets/" + name + ".poly";
        path.remove(QString::fromStdString(str));

        str = name + ".txt";
        path.append(QString::fromStdString(str));

        std::string fileName = path.toStdString();

        // Vertex adjustment.
        setPreprocessingStatus(true);

        // Set maximum number of subdivision steps.
        setMaxNumSubSteps(maxSteps);

        // Choose barycentric morphing method.
        toggleMorphingMethod(defaultMethod);

        // Store average time.
        std::vector<double> averageTime(getMaxNumSubSteps(), 0.0);

        /// Subdivision:

        // Average time.
        for(int i = 0; i < timesToSubdivide; ++i) {

            // Preprocess subdivision.
            morph();

            for(int j = 0; j < getMaxNumSubSteps(); ++j) {
                setNumSubSteps(j+1); // no need for preprocess
                subdivideUp(); // no need for preprocess

                subdivide(); // no need for preprocess

                // setMaxNumSubSteps(j);
                // morph();
                // setNumSubSteps(0);
                // setCurrMaxNumSubSteps(0);
                // _vertices.clear();
                // _halfedges.clear();
                // triMesh().clear();
                // triangulate();

                averageTime[j] += getMorphingTime();
                // std::cout << "Time for level " << j + 1 << " is " << averageTime[j] << std::endl;
            }
            // std::cout << std::endl;

            std::cout << "\nIteration " << i << " is finished!" << std::endl;

            // Clear all stuff.
            setNumSubSteps(0);
            setCurrMaxNumSubSteps(0);

            _vertices.clear();
            _halfedges.clear();

            triMesh().clear();
            triangulate();
        }

        // Save time in a file.
        for(int i = 0; i < getMaxNumSubSteps(); ++i) averageTime[i] /= timesToSubdivide;

        std::ofstream saveFile(fileName.c_str(), std::ios_base::out);

        std::cout.precision(30);

        if(!saveFile) {
            std::cout << "\nError saving file!" << std::endl;
            exit(EXIT_FAILURE);
        }

        saveFile << "Time for " << (defaultMorphingMethod ? "default" : "barycentric") << " subdivision method (" << timesToSubdivide << " iterations)." << std::endl;
        for(int i = 0; i < getMaxNumSubSteps(); ++i)
            saveFile << averageTime[i] <<  std::endl;

        saveFile.close();

        // Status.
        std::cout << "\nSUCCESS\n" << std::endl;
        exit(EXIT_SUCCESS);
    }

} // namespace warpit
