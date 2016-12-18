//* Authors: Dmitry Anisimov
//* program.hpp
//* With any questions please contact: danston@ymail.com

#ifndef PROGRAM_HPP
#define PROGRAM_HPP

// Local includes.
#include "support/texture.hpp"
#include "support/triangulator.hpp"

// Qt includes.
#include <QPoint>

// warpit namespace.
namespace warpit {

    // The main program class.
    class Program
    {
    public:
        // Create an instance of the main program.
        static Program& instance()
        {
            static Program program = Program();
            return program;
        }

        // Return the source polygon.
        inline Polygon& sourcePolygon()
        {
            return _sourcePolygon;
        }

        // Return the const source polygon.
        inline const Polygon& sourcePolygon() const
        {
            return _sourcePolygon;
        }

        // Return the target polygon.
        inline Polygon& targetPolygon()
        {
            return _targetPolygon;
        }

        // Return the const target polygon.
        inline const Polygon& targetPolygon() const
        {
            return _targetPolygon;
        }

        // Return the triangle mesh.
        inline TriMesh& triMesh()
        {
            return _triMesh;
        }

        // Return the const triangle mesh.
        inline const TriMesh& triMesh() const
        {
            return _triMesh;
        }

        // Return the number of virtual vertices.
        inline int numVirtVertices() const
        {
             return _virtualVertices.size();
        }

        // Return virtual vertices.
        inline std::vector<Vertex>& virtualVertices()
        {
            return _virtualVertices;
        }

        // Return const virtual vertices.
        inline const std::vector<Vertex>& virtualVertices() const
        {
            return _virtualVertices;
        }

        // Return stored halfedges.
        inline std::vector<Halfedge>& halfedges(const int index)
        {
            return _halfedges[index];
        }

        // Return stored const halfedges.
        inline const std::vector<Halfedge>& halfedges(const int index) const
        {
            return _halfedges[index];
        }

        // Return vertices.
        inline std::vector<Vertex>& vertices(const int index)
        {
            return _vertices[index];
        }

        // Return const vertices.
        inline const std::vector<Vertex>& vertices(const int index) const
        {
            return _vertices[index];
        }

        // Remove all virtual vertices that lie outside the source polygon.
        void cleanVirtualVertices();

        // Clear the program data.
        void clear()
        {
             _triMesh.clear();
             _sourcePolygon.clear();
             _targetPolygon.clear();

             _vertices.clear();
             _halfedges.clear();

             _virtualVertices.clear();
        }

        // Triangulate the source polygon.
        void triangulate();

        // Subdivide the mesh.
        void subdivide();

        // Subdivide the mesh up.
        void subdivideUp();

        // Subdivide the mesh down.
        void subdivideDown();

        // Perform the morphing.
        void morph();

        // Set the number of subdivision steps to zero.
        inline void setNumSubStepsToDefault()
        {
            numSubSteps = 0;
        }

        // Set new number of subdivision steps.
        inline void setNumSubSteps(const int newNumSubSteps)
        {
            numSubSteps = newNumSubSteps;
        }

        // Set new current max number of subdivision steps.
        inline void setCurrMaxNumSubSteps(const int newCurrMaxSteps)
        {
            currMaxSteps = newCurrMaxSteps;
        }

        // Set maximum number of subdivision steps.
        inline void setMaxNumSubSteps(const int newMaxNumSubSteps)
        {
            maxNumSubSteps = newMaxNumSubSteps;
        }

        // Get current number of subdivision steps.
        inline int getNumSubSteps() const
        {
            return numSubSteps;
        }

        // Get maximum number of subdivision steps.
        inline int getMaxNumSubSteps() const
        {
            return maxNumSubSteps;
        }

        // Create a regular grid.
        void createRegularGrid(const double width, const double height, const int gridSize = 16);

        // Return points of a regular grid.
        inline std::vector<std::vector<QPoint> >& gridPoints()
        {
            return _gridPoints;
        }

        // Return constant points of a regular grid.
        inline const std::vector<std::vector<QPoint> >& gridPoints() const
        {
            return _gridPoints;
        }

        // Set new preprocessing status.
        inline void setPreprocessingStatus(const bool newPreprocessingStatus)
        {
            preprocessingStatus = newPreprocessingStatus;
        }

        // Enumeration with different program modes.
        enum ProgramMode { INITIAL, PREPROCESSING, MORPHING };

        // Get current program mode.
        inline ProgramMode getProgramMode() const
        {
            return programMode;
        }

        // Set new program mode.
        inline void setProgramMode(const ProgramMode newProgramMode)
        {
            programMode = newProgramMode;
        }

        // Change morphing method.
        inline void toggleMorphingMethod(const bool changeMethod)
        {
            defaultMorphingMethod = changeMethod;
        }

        // Return current morphing time.
        inline double getMorphingTime() const
        {
            return morphingTime;
        }

        // Create a texture coordinate from a native coordinate.
        inline void texCoord(Vector2d &v, const Vector2d &min, const Vector2d &max) const
        {
            v.x -= min.x; v.x /= (max.x - min.x);
            v.y -= min.y; v.y /= (max.y - min.y);
        }

        // Run program without interface.
        void run(QString &path, const int timesToSubdivide, const int maxSteps, bool defaultMethod);

    private:
        // The default constructor.
        Program();

        // Create a source polygon.
        Polygon _sourcePolygon;

        // Create a target polygon.
        Polygon _targetPolygon;

        // Create a triangle mesh.
        // In addition, create two supplementary data storages.
        // The first one stores vertices for each subdivision step up to max.
        // The second one stores the connectivity for each subdivision step up to max.
        TriMesh _triMesh;
        std::vector<std::vector<Vertex> >   _vertices;
        std::vector<std::vector<Halfedge> > _halfedges;

        // Stores virtual vertices.
        std::vector<Vertex> _virtualVertices;

        // Stores a regular grid for the viewer.
        std::vector< std::vector<QPoint> > _gridPoints;

        // Number of subdivision steps.
        int numSubSteps;

        // Maximum number of subdivision steps.
        int maxNumSubSteps;

        // Current maximum number of subdivision steps.
        int currMaxSteps;

        // Program mode.
        ProgramMode programMode;

        // Preprocessing status.
        bool preprocessingStatus;

        // Morphing method. True - the subdivision method, false - that one with barycentric coordinates.
        bool defaultMorphingMethod;

        // Epsilon.
        double eps;

        // Time to morph the mesh.
        double morphingTime;

        // Initial mesh.
        TriMesh initMesh;

        // Initiate the mesh.
        void initiateMesh();

        // Preprocess the mesh.
        void preprocessMesh();

        // Set mesh initial data.
        void setMeshInitialData();

        // Fix barycentric coordinates for virtual vertices.
        void fixVirtualCoordinates();

        // Morph the mesh with subdivision weights.
        void morphWithSubdivisionWeights();

        // Morph the mesh with barycentric coordinates.
        void morphWithBarycentricCoordinates();

        // Update texture coordinates for each subdivision step.
        void updateTextureCoordinates();

        // Set initial data for a target mesh.
        void setTargetInitialData();

        // Compute triangle barycentric coordinates.
        void computeTriangleCoordinates(const Vertex &v1, const Vertex &v2, const Vertex &v3, const int meshInd, std::vector<double> &coordinates) const;

        // Compute the cross product.
        double crossProduct(const Vertex &v1, const Vertex &v2, const Vertex &v3) const;

        // Test the equality between both above methods to morph the mesh.
        void compareMorphingMethods();

        // Check properties of barycentric coordinates.
        void checkBarycentricProperties();
    };

} // namespace warpit

#endif // PROGRAM_HPP
