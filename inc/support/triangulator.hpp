//* Authors: Teseo Schneider and Dmitry Anisimov (c) 2016-2017
//* With any questions please contact: danston@ymail.com

// This is the wrapper for the triangle library.

#ifndef TRIANGULATOR_HPP
#define TRIANGULATOR_HPP

#ifdef SINGLE
#define REAL float
#else
#define REAL double
#endif

// STL includes.
#include <vector>
#include <sstream>

// Local includes.
#include "trimesh.hpp"
#include "polygon.hpp"
#include "triangle.hpp"

// warpit namespace
namespace warpit {

    static const double SQRT3_OVER_4 = 0.43301270189222;

    // The triangulator class.
    class Triangulator
    {
    public:
        // Set the flags and boundaries.
        inline void allowBoundaryRefinement(const bool value)                      { _allowBoundaryRefinement = value; }
        inline void allowEdgeRefinement(const bool value)                          { _allowEdgeRefinement = value; }
        inline void constrainAngle(const bool value)                               { _constrainAngle = value; }
        inline void setPlanarGraph(const bool value)                               { _planarGraph = value; }
        inline void recordOrder(const bool value)                                  { _recordOrder = value; }
        inline void setMaxEdgeLength(const double value)                           { if(value > 0.0) _maxArea = value * value * SQRT3_OVER_4; else _maxArea = -1.0; }
        inline void allowSteinerPoints(const bool value)                           { _allowSteinerPoints = value; }
        inline void setPolygon(const Polygon &polygon)                             { _polygon = polygon; }
        inline void setVirtualVertices(const std::vector<Vertex> &virtualVertices) { _virtualVertices = virtualVertices; }
        inline void setEdgeTag(const bool value)                                   { _setEdgeTag = value; }

        // The default constructor.
        Triangulator()
        {
            resetFlags();
        }

        // Constructor.
        Triangulator(const Polygon &polygon, const double maxEdgeLength = -1.0, const bool refineBoundary = false)
        {
            resetFlags();
            setPolygon(polygon);

            allowBoundaryRefinement(refineBoundary);
            setMaxEdgeLength(maxEdgeLength);
        }

        // Generate the mesh.
        void generate(TriMesh &triMesh)
        {
            triangulateio out;
            triangulateio in;

            assignValues(in);
            triangulateValues(in, out);

            freeIn(in);

            setTriangulation(out, triMesh);
            freeOut(out);

            if(!_setEdgeTag) return;
        }

    private:
        // Global variables.
        typedef struct triangulateio triangulateio;

        Polygon _polygon;
        std::vector<Vertex> _virtualVertices;

        bool _allowBoundaryRefinement, _allowEdgeRefinement;
        bool _constrainAngle, _planarGraph;

        double _maxArea;
        int    _allowSteinerPoints;
        bool   _setEdgeTag, _recordOrder, _verbose;

        // Reset all flags to the default values.
        void resetFlags()
        {
            _allowBoundaryRefinement = true;
            _allowEdgeRefinement     = true;
            _constrainAngle          = true;
            _planarGraph             = true;
            _maxArea                 = -1.0;
            _allowSteinerPoints      = 0;
            _setEdgeTag              = false;
            _recordOrder             = false;
            _verbose                 = false;
        }

        // Set all the points and markers.
        void assignValues(triangulateio &in) const
        {
            const int numV = _polygon.numVertices();
            const int numP = numV + _virtualVertices.size();
            const int numE = numV;

            in.numberofpoints = numP;
            in.numberofpointattributes = _recordOrder ? 1 : 0;

            in.pointlist = (double *) malloc(in.numberofpoints * 2 * sizeof(double));
            in.pointmarkerlist = (int *) malloc(in.numberofpoints * sizeof(int));

            if(_recordOrder)
                in.pointattributelist = (double *) malloc(in.numberofpointattributes * in.numberofpoints * sizeof(double));

            in.numberofsegments = numE;
            in.segmentlist = (int *) malloc(in.numberofsegments * 2 * sizeof(int));
            in.segmentmarkerlist = (int *) malloc(in.numberofsegments * sizeof(int));

            if(_recordOrder) {
                for(int i = 0; i < in.numberofpoints; ++i)
                    in.pointattributelist[i] = i;
            }

            for(int i = 0; i < numV; ++i)
            {
                in.pointlist[2 * i]     = _polygon.vertices()[i].p().x;
                in.pointlist[2 * i + 1] = _polygon.vertices()[i].p().y;

                in.pointmarkerlist[i]   = 1;
                in.segmentmarkerlist[i] = 1;

                in.segmentlist[2 * i] = i;

                if(i == numV - 1) in.segmentlist[2 * i + 1] = 0;
                else in.segmentlist[2 * i + 1] = i + 1;
            }

            for(int i = numV; i < numP; ++i)
            {
                in.pointlist[2 * i]     = _virtualVertices[i - numV].p().x;
                in.pointlist[2 * i + 1] = _virtualVertices[i - numV].p().y;

                in.pointmarkerlist[i] = _polygon.belongsToBoundary(_virtualVertices[i - numV]) ? 1 : 0;
            }

            in.numberofholes   = 0;
            in.numberofregions = 0;
        }

        // Triangulate the set of input points.
        void triangulateValues(triangulateio &mid, triangulateio &out)
        {
            out.pointlist = (double *) NULL;
            out.pointmarkerlist = (int *) NULL;

            if(_recordOrder)
                out.pointattributelist = (double *) NULL;

            out.trianglelist = (int *) NULL;
            out.triangleattributelist = (double *) NULL;

            out.segmentlist = (int *) NULL;
            out.segmentmarkerlist = (int *) NULL;

            out.edgelist=(int *) NULL;
            out.edgemarkerlist = (int *) NULL;

            std::stringstream buf;
            buf.precision(100);
            buf.setf(std::ios::fixed, std::ios::floatfield);

            if(!_verbose)
                buf << "Q";
            else
                buf << "V";

            buf << "ez";

            if(!_allowEdgeRefinement)
                buf << "YY";
            else if(!_allowBoundaryRefinement)
                buf << "Y";

            if(_constrainAngle)
                buf << "q";

            if(_planarGraph)
                buf << "p";

            if(_allowSteinerPoints >= 0)
                buf << "S" << _allowSteinerPoints;

            if(_maxArea > 0)
                buf<< "a" << _maxArea;

            char* str = new char[buf.str().size() + 1];
            strcpy(str, buf.str().c_str());

            triangulate(str, &mid, &out, (triangulateio *) NULL);

            delete[] str;
        }

        // Set the obtained triangulation to the provided mesh and initialize this mesh.
        void setTriangulation(triangulateio &triangulation, TriMesh &triMesh)
        {
            triMesh.clear();

            std::vector<Vertex> vertices;
            std::vector<Face>   faces;

            Vertex newVertex;
            for(int i = 0; i < triangulation.numberofpoints; ++i) {
                newVertex.p() = Vector2d(triangulation.pointlist[i * 2], triangulation.pointlist[i * 2 + 1]);
                vertices.push_back(newVertex);
            }

            Face newFace;
            for(int i = 0; i < triangulation.numberoftriangles; ++i) {
                newFace.v[0] = triangulation.trianglelist[i * 3];
                newFace.v[1] = triangulation.trianglelist[i * 3 + 1];
                newFace.v[2] = triangulation.trianglelist[i * 3 + 2];

                faces.push_back(newFace);
            }

            triMesh.initialize(vertices, faces);
        }

        // Free memory in the input triangulation.
        void freeIn(triangulateio &in)
        {
            free(in.pointlist);
            free(in.pointmarkerlist);

            free(in.segmentlist);
            free(in.segmentmarkerlist);

            if(_recordOrder)
                free(in.pointattributelist);
        }

        // Free memory in the output triangulation.
        void freeOut(triangulateio &out)
        {
            free(out.pointlist);
            free(out.pointmarkerlist);
            if(_recordOrder)
                free(out.pointattributelist);

            free(out.trianglelist);
            free(out.triangleattributelist);

            free(out.segmentlist);
            free(out.segmentmarkerlist);

            free(out.edgelist);
            free(out.edgemarkerlist);
        }
    };

} // namespace warpit

#endif // TRIANGULATOR_HPP
