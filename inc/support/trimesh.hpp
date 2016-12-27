//* Authors: Kai Hormann and Dmitry Anisimov (c) 2016-2017
//* With any questions please contact: danston@ymail.com

#ifndef TRIMESH_HPP
#define TRIMESH_HPP

// STL includes.
#include <map>
#include <fstream>

// Local includes.
#include "face.hpp"
#include "vertex.hpp"
#include "halfedge.hpp"

// warpit namespace
namespace warpit {

    // Triangle mesh with the mixed data structure.
    class TriMesh
    {
    public:
        // Constructor.
        TriMesh() : vertex(), halfedge(), face(), eps(1.0e-14) { }

        // Check if mesh contains any vertices.
        inline bool empty() const
        {
            if(numVertices() == 0) return true;
            return false;
        }

        // Clear mesh.
        inline void clear()
        {
            vertex.clear();
            halfedge.clear();
            face.clear();
        }

        // Return vertices of the mesh.
        inline std::vector<Vertex>& vertices()
        {
            return vertex;
        }

        // Return const vertices of the mesh.
        inline const std::vector<Vertex>& vertices() const
        {
            return vertex;
        }

        // Return halfedges of the mesh.
        inline std::vector<Halfedge>& halfedges()
        {
            return halfedge;
        }

        // Return const halfedges of the mesh.
        inline const std::vector<Halfedge>& halfedges() const
        {
            return halfedge;
        }

        // Return faces of the mesh.
        inline std::vector<Face>& faces()
        {
            return face;
        }

        // Return const faces of the mesh.
        inline const std::vector<Face>& faces() const
        {
            return face;
        }

        // Return number of vertices in the mesh.
        inline int numVertices() const
        {
            return vertex.size();
        }

        // Return number of halfedges in the mesh.
        inline int numHalfedges() const
        {
            return halfedge.size();
        }

        // Return number of faces in the mesh.
        inline int numFaces() const
        {
            return face.size();
        }

        // Return number of edges in the mesh.
        inline int numEdges() const
        {
            return numVertices() + (numHalfedges() / 3) - 1;
        }

        // Epsilon used internally in the class.
        inline double epsilon() const
        {
            return eps;
        }

        // Set new user-defined epsilon = tolerance.
        inline void setTolerance(const double newTolerance)
        {
            eps = newTolerance;
        }

        // Subdivide mesh.
        inline int subdivide(const int timesToSubdivide = 1)
        {
            for(int i = 0; i < timesToSubdivide; ++i) subdivideMesh();
            return 0;
        }

        // Subdivide the morphed mesh.
        inline int subdivideMorphed(const std::vector<Halfedge> &halfedgesPrev,
                                    const std::vector<Halfedge> &halfedgesNext,
                                    const std::vector<Vertex>    &verticesPrev,
                                    std::vector<Vertex>          &verticesNext)
        {
            subdivideMorphedMesh(halfedgesPrev, halfedgesNext, verticesPrev, verticesNext);
            return 0;
        }

        // Initialize mesh.
        int initialize(const std::vector<Vertex> &vertices, const std::vector<Face> &faces)
        {
            // Set defaults.
            const int numV  = vertices.size();
            const int numF  = faces.size();
            const int numHE = 3 * numF;

            face.resize(numF);
            vertex.resize(numV);
            halfedge.resize(numHE);

            // Copy vertices.
            for(int i = 0; i < numV; ++i) {
                // Assertions.
                assert(vertices[i].val == 0);
                assert(vertices[i].out == -1);
                assert(vertices[i].type == 3);
                assert(vertices[i].alpha == -1.0);

                vertex[i] = vertices[i];
            }

            // Copy faces.
            for(int i = 0; i < numF; ++i) {
                // Assertions.
                assert(faces[i].f[0] == -1);
                assert(faces[i].f[1] == -1);
                assert(faces[i].f[2] == -1);
                assert(faces[i].f[3] == -1);

                face[i] = faces[i];
            }

            // Copy vertex indices of each triangle in the mesh.
            std::vector<int> vertexIndices;
            vertexIndices.resize(numHE);

            int ind = 0;
            for(int i = 0; i < numF; ++i)
            {
                vertexIndices[ind++] = faces[i].v[0];
                vertexIndices[ind++] = faces[i].v[1];
                vertexIndices[ind++] = faces[i].v[2];
            }

            int idxV = 0; int idxE = 0; int base;
            for(int i = 0; i < numF; i++) {
                base = idxE;
                for(int j = 0; j < 3; j++) {
                    halfedge[idxE].prev   = base + (j + 2) % 3;
                    halfedge[idxE].next   = base + (j + 1) % 3;
                    halfedge[idxE++].dest = vertexIndices[idxV++];
                }
            }

            // Build the halfedge connectivity.
            buildHEConnectivity();

            // Return a flag.
            return 0;
        }

        // Initialize barycentric coordinates.
        void initializeBarycentricCoordinates()
        {
            const int numV = numVertices();

            // Set unit coordinates to all initial vertices.
            for(int i = 0; i < numV; ++i) {
                vertex[i].b().resize(numV);
                for(int j = 0; j < numV; ++j) vertex[i].b()[j] = 0.0;
                vertex[i].b()[i] = 1.0;
            }
        }

        // Get a one ring neighbourhood around the vertex.
        void getRing(const int meshIndex, std::vector<int> &neighbours) const
        {
            int n = 0, prev, curr, next;

            // Find neighbours of the vertex.
            curr = next = vertex[meshIndex].out;

            assert(vertex[meshIndex].val > 0);
            neighbours.resize(vertex[meshIndex].val);
            do {
                neighbours[n] = halfedge[curr].dest;

                prev = halfedge[curr].prev;
                curr = halfedge[prev].neigh;
                n++;
            } while((curr >= 0) && (curr != next));

            // Add one more neighbour for boundary vertices.
            if(curr < 0) {
                curr          = halfedge[prev].prev;
                neighbours[n] = halfedge[curr].dest;
            }
        }

        // Create faces in the mesh.
        void createFaces()
        {
            const int numHE = numHalfedges();
            const int numF  = numHE / 3;

            face.resize(numF);

            assert(numHE != 0 && numF != 0);

            // Create triangular faces.
            int idxE = 0;
            for(int i = 0; i < numF; ++i)
                for(int j = 0; j < 3; ++j)
                    face[i].v[j] = halfedge[halfedge[idxE++].prev].dest;
        }

        // Preprocess this mesh.
        inline int preprocess()
        {
            // Do a ternary step.
            const int failure = ternaryStep();

            // Save mesh vertices and faces before the corner modification
            // in order to sample linearly barycentric coordinates later.
            std::vector<Vertex> prevVertex = vertex;
            std::vector<Face>   prevFace   = face;

            // Modify 1-ring neighbourhood around each concave corner.
            if(!failure) return modifyCorners(prevVertex, prevFace);

            return failure;
        }

        // Barycentric properties.

        // Partition of unity.
        bool partitionOfUnity() const
        {
            const int numV = numVertices();
            const int numC = vertex[0].b().size();

            double sum;
            bool state = true;

            for(int i = 0; i < numV; ++i) {
                sum = 0.0;
                for(int j = 0; j < numC; ++j) sum += vertex[i].b()[j];

                if(std::fabs(sum - 1.0) > eps) {
                    std::cerr << "WARNING: Partition of unity difference is " << std::fabs(sum - 1.0) << ";" << std::endl;
                    state = false;
                }
            }

            if(!state) std::cerr << "\n";
            return state;
        }

        // Linear precision.
        bool linearPrecision(const std::vector<Vertex> &initialVertex) const
        {
            const int numV = numVertices();
            const int numC = vertex[0].b().size();

            bool state = true;

            for(int i = 0; i < numV; ++i) {
                Vector2d sum;
                for(int j = 0; j < numC; ++j)
                    sum += vertex[i].b()[j]*initialVertex[j].p();

                if(std::fabs(sum.x - vertex[i].p().x) > eps || std::fabs(sum.y - vertex[i].p().y) > eps) {
                    std::cerr << "WARNING: Linear precision difference is (" << std::fabs(sum.x - vertex[i].p().x) << ", "
                                                                             << std::fabs(sum.y - vertex[i].p().y) << ");"
                                                                             << std::endl;
                    state = false;
                }
            }

            if(!state) std::cerr << "\n";
            return state;
        }

        // Boundedness.
        bool boundedness() const
        {
            const int numV = numVertices();
            const int numC = vertex[0].b().size();

            bool state = true;

            for(int i = 0; i < numV; ++i) {
                for(int j = 0; j < numC; ++j) {
                    if(vertex[i].b()[j] < 0.0 || vertex[i].b()[j] > 1.0) {
                        std::cerr << "WARNING: Found a value out of range [0,1] that is " << vertex[i].b()[j] << ";" << std::endl;
                        state = false;
                    }
                }
            }

            if(!state) std::cerr << "\n";
            return state;
        }

        // Load obj mesh from a file.
        void load(const std::string &fileName)
        {
            if(fileName.empty()) return;

            std::ifstream readFile(fileName.c_str(), std::ios_base::in);

            if(!readFile) {
                std::cout << "\nError reading .obj file!" << std::endl;
                exit(EXIT_FAILURE);
            }

            std::string dummy;

            double firstCoord, secondCoord, thirdCoord;
            int    firstIndex, secondIndex, thirdIndex;

            Vertex newVertex;
            Face newFace;

            while(!readFile.eof()) {
                readFile >> dummy;

                if(dummy == "v") {
                    readFile >> firstCoord >> secondCoord >> thirdCoord;

                    newVertex.p().x = firstCoord;
                    newVertex.p().y = secondCoord;

                    vertex.push_back(newVertex);
                }

                if(dummy == "f") {
                    readFile >> firstIndex >> secondIndex >> thirdIndex;

                    newFace.v[0] = firstIndex  - 1;
                    newFace.v[1] = secondIndex - 1;
                    newFace.v[2] = thirdIndex  - 1;

                    face.push_back(newFace);
                }
            }
            readFile.close();
        }

    private:
        // Mesh.
        std::vector<Vertex>   vertex;   // stores vertices  of the mesh
        std::vector<Halfedge> halfedge; // stores halfedges of the mesh
        std::vector<Face>     face;     // stores faces     of the mesh
        std::vector<Vector2d> morphed;  // stores morphed positions of the mesh

        // Epsilon.
        double eps;

        // Build the halfedge connectivity.
        void buildHEConnectivity()
        {
            const int numV  = numVertices();
            const int numHE = numHalfedges();

            typedef std::map<int,int> halfedgeList;
            std::vector<halfedgeList> halfedgeTable(numV);
            halfedgeList* destList;
            int source, dest;

            // Build the halfedge connectivity.
            for(int i = 0; i < numHE; ++i) {
                // Source and destination of the current halfedge.
                source = halfedge[halfedge[i].prev].dest;
                dest   = halfedge[i].dest;

                // Is halfedge from destination to source already in the edge table?
                destList = &(halfedgeTable[dest]);
                halfedgeList::iterator it = destList->find(source);
                if(it != destList->end()) {
                    halfedge[i].neigh = it->second;
                    halfedge[it->second].neigh = i;
                    destList->erase(it);
                }
                else {
                    // Put a halfedge in the edge table.
                    halfedgeTable[source].insert(std::make_pair(dest,i));
                }
            }

            // Determine valency and some outgoing halfedge for each vertex.
            // Mark and count the boundary vertices.
            Vertex* destV;
            for(int i = 0; i < numHE; ++i) {
                // Destination vertex of the current halfedge.
                destV = &(vertex[halfedge[i].dest]);
                // Increase valency of destination vertex.
                destV->val++;
                // Take next of the current halfedge as the outgoing halfedge.
                destV->out = halfedge[i].next;

                if(halfedge[i].neigh < 0) {
                    // NOTE: boundary vertices have one more neighbour than the outgoing edges.
                    destV->type = 2; // Flat vertices
                    destV->val++;
                }
            }

            // Get the "rightmost" outgoing halfedge for boundary vertices.
            for(int i = 0; i < numV; ++i) {
                if(vertex[i].type != 3) {
                    while (halfedge[vertex[i].out].neigh >= 0) {
                        // Move the outgoing halfedge "one to the right".
                        vertex[i].out = halfedge[halfedge[vertex[i].out].neigh].next;
                    }
                }
            }

            // Update type of each boundary vertex in the mesh.
            Vertex* prevV; Vertex* nextV;
            for(int i = 0; i < numHE; ++i) {
                if(halfedge[i].neigh < 0) {
                    destV = &(vertex[halfedge[i].dest]);

                    const int prev = destV->out;
                    prevV  = &(vertex[halfedge[prev].dest]);

                    const int next = halfedge[i].prev;
                    nextV  = &(vertex[halfedge[next].dest]);

                    defineBoundaryVertexType(prevV, destV, nextV);
                }
            }
        }

        // Define type of a vertex on the mesh boundary.
        inline void defineBoundaryVertexType(Vertex* prevV, Vertex* destV, Vertex* nextV) const
        {
            const double x1 = prevV->p().x; const double y1 = prevV->p().y;
            const double x2 = destV->p().x; const double y2 = destV->p().y;
            const double x3 = nextV->p().x; const double y3 = nextV->p().y;

            // Compute the following determinant:
            //
            //       | x1 y1 1 |
            // det = | x2 y2 1 |
            //       | x3 y3 1 |
            //
            // det = 0 if three points are collinear;
            // det > 0 if they create a concave corner;
            // det < 0 if they create a convex  corner;

            const double det = x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2);

            // Compute the external signed angle for each corner.
            const double dotProd  = (x1 - x2) * (x3 - x2) + (y1 - y2) * (y3 - y2);
            const double extAngle = atan2(det, dotProd); // positive in case of concave and negative in case of convex corner

            const double eps = epsilon();

            // std::cout << "det = " << det << std::endl;

            if(det > eps) { // CONCAVE corner
                destV->type  = 1;
                destV->alpha = extAngle; // external angle
            } else { // CONVEX corner
                destV->type  = 0;
                destV->alpha = 2.0 * M_PI + extAngle; // internal angle
            }

            // std::cout << "Corner: " << (destV->type ? "CONCAVE\n" : "CONVEX\n") << std::endl;
            // std::cout << "Angle: "  << destV->alpha << std::endl;
        }

        // Ternary subdivision step.
        int ternaryStep()
        {
            const int numF  = numFaces();
            const int numV  = numVertices();
            const int numHE = numHalfedges();
            const int numE  = numEdges();

            assert(numF != 0 && numV != 0 && numHE != 0 && numE != 0);

            std::vector<Vertex> tmpVertices(numV + numF + 2 * numE);
            std::vector<Face>   tmpFaces(9 * numF);

            // Create temporary vertices.

            // Old vertices.
            for(int i = 0; i < numV; ++i) {
                tmpVertices[i].p()  = vertex[i].p();
                tmpVertices[i].pm() = vertex[i].pm();
                tmpVertices[i].b()  = vertex[i].b();
            }

            // New face vertices.
            for(int i = 0; i < numF; ++i) createFaceBarycenter(face[i].v[0], face[i].v[1], face[i].v[2], tmpVertices[numV + i]);

            // New edge vertices.
            int count = 0;
            for(int i = 0; i < numHE; ++i) {

                const int i_ = halfedge[i].neigh;
                if(i > i_) {

                    const int i1 = halfedge[halfedge[i].prev].dest;
                    const int i2 = halfedge[i].dest;

                    assert(i1 != -1 && i2 != -1);

                    trisectEdge(i1, i2, tmpVertices[numV + numF + count], tmpVertices[numV + numF + count + 1]);
                    count += 2;
                }
            }

            // Create temporary faces.

            Vertex tmpV1, tmpV2;
            std::vector<int> neighbours(6);

            // Create new faces.
            count = 0;
            for(int i = 0; i < numF; ++i) {

                const int i1 = face[i].v[0];
                const int i2 = face[i].v[1];
                const int i3 = face[i].v[2];

                assert(i1 != -1 && i2 != -1 && i3 != -1);

                trisectEdge(i1, i2, tmpV1, tmpV2);
                neighbours[0] = findIndex(tmpVertices, tmpV1);
                neighbours[1] = findIndex(tmpVertices, tmpV2);

                trisectEdge(i2, i3, tmpV1, tmpV2);
                neighbours[2] = findIndex(tmpVertices, tmpV1);
                neighbours[3] = findIndex(tmpVertices, tmpV2);

                trisectEdge(i3, i1, tmpV1, tmpV2);
                neighbours[4] = findIndex(tmpVertices, tmpV1);
                neighbours[5] = findIndex(tmpVertices, tmpV2);

                for(int j = 0; j < 6; ++j) {
                    tmpFaces[count].v[0]   = numV + i;
                    tmpFaces[count].v[1]   = neighbours[j];
                    tmpFaces[count++].v[2] = neighbours[(j + 1) % 6];
                }

                tmpFaces[count].v[0] = face[i].v[0]; tmpFaces[count].v[1] = neighbours[0]; tmpFaces[count++].v[2] = neighbours[5];
                tmpFaces[count].v[0] = face[i].v[1]; tmpFaces[count].v[1] = neighbours[2]; tmpFaces[count++].v[2] = neighbours[1];
                tmpFaces[count].v[0] = face[i].v[2]; tmpFaces[count].v[1] = neighbours[4]; tmpFaces[count++].v[2] = neighbours[3];
            }

            // Clear old and initialize new mesh.
            clear();
            const int failure = initialize(tmpVertices, tmpFaces);

            if(failure) std::cerr << "ERROR: Mesh initialization after a ternary subdivision step failed!\n" << std::endl;

            // Copy back barycentric coordinates and some flags.
            if(!failure) {
                for(int i = 0; i < numVertices(); ++i) {
                    vertex[i].b()  = tmpVertices[i].b();
                    if(i > numV - 1 && vertex[i].type != 3) {
                        vertex[i].type  = 2;
                        vertex[i].alpha = -1.0;
                    }
                }
                return 0;
            }

            return failure;
        }

        // Create the barycenter point in a triangle.
        inline void createFaceBarycenter(const int i1, const int i2, const int i3, Vertex &barycenter) const
        {
            // Find the triangle barycenter.
            barycenter  = vertex[i1];
            barycenter += vertex[i2];
            barycenter += vertex[i3];

            barycenter *= (1.0 / 3.0);

            // Update its flags.
            barycenter.val   = 0;
            barycenter.out   = -1;
            barycenter.type  = 3;
            barycenter.alpha = -1.0;
        }

        // Trisect an edge.
        inline void trisectEdge(const int i1, const int i2, Vertex &newV1, Vertex &newV2) const
        {
            // Trisect the edge.
            const double third = 1.0 / 3.0;

            newV1  = vertex[i1];
            newV1 *= (1.0 - third);
            newV1 += third * vertex[i2];

            newV2  = vertex[i1];
            newV2 *= third;
            newV2 += (1.0 - third) * vertex[i2];

            // Update flags.
            newV1.val   = 0;    newV2.val  = 0;
            newV1.out   = -1;   newV2.out  = -1;
            newV1.type  = 3;    newV2.type = 3;
            newV1.alpha = -1.0; newV2.alpha = -1.0;
        }

        // Find vertex index.
        inline int findIndex(const std::vector<Vertex> &vertices, const Vertex &v) const
        {
            const int numV = vertices.size();

            for(int i = 0; i < numV; ++i)
                if(vertices[i].p() == v.p()) return i;

            // Pointer cannot be here.
            bool wrongLocation = false;
            assert(wrongLocation);
            if(!wrongLocation) wrongLocation = true;

            return -1;
        }

        // Modify corners of the mesh.
        int modifyCorners(std::vector<Vertex> &prevVertex, std::vector<Face> &prevFace)
        {
            const int numV = numVertices();
            std::vector<int> neighbours;

            assert(numV != 0);

            int failure = 0;

            // Go through all the mesh vertices.
            for(int i = 0; i < numV; ++i) {
                if(vertex[i].type == 1) { // Concave vertex

                    // Find 1-ring neighbourhood of the corner vertex.
                    getRing(i, neighbours);
                    const int nSize = neighbours.size();

                    // Find the internal angle for the corner and split it into equal parts.
                    const double intAngle = 2.0 * M_PI - vertex[i].alpha;
                    const double angle    = intAngle / (nSize - 1);

                    // Find the closest neighbour to the corner.
                    double minLen = std::numeric_limits<double>::max();
                    for(int j = 0; j < nSize; ++j) {
                        const double len = (vertex[i].p()  - vertex[neighbours[j]].p()).length();
                        minLen = std::min(minLen, len);
                    }

                    // Rotate the first neighbour of the corner vertex by the angle 'angle'
                    // to get equally-distant and distributed mesh vertices around this corner.
                    double tmpAngle = 0.0;
                    Vector2d &firstV = vertex[neighbours[0]].p();

                    Vector2d v;
                    for(int j = 0; j < nSize; ++j) {
                        tmpAngle = j * angle;

                        v = firstV.translated(-vertex[i].p());
                        v = v.rotated(tmpAngle);

                        const double len = v.length();
                        const double scaleFactor = minLen / len;

                        v = v.scaled(Vector2d(scaleFactor, scaleFactor));
                        v = v.translated(vertex[i].p());

                        vertex[neighbours[j]].p() = v;
                    }

                    // Scale mesh vertices around the corner to avoid interior fold-overs.
                    failure = scaleVertices(i, neighbours);
                    if(failure) return failure;

                    // Fix barycentric coordinates for all the neighbours of the corner vertex.
                    failure = fixCoordinates(neighbours, prevVertex, prevFace);
                    if(failure) return failure;
                }
            }

            return failure;
        }

        // Scale mesh vertices around the corner.
        inline int scaleVertices(const int centerInd, const std::vector<int> &neighbours)
        {
            const int nSize = neighbours.size();

            int count = 0;
            const int maxCount = 25;

            // Until all the interior fold-overs are gone, scale neighbours down by 1/2.
            while(foldover(neighbours)) {
                for(int i = 0; i < nSize; ++i) {
                    vertex[neighbours[i]] += vertex[centerInd];
                    vertex[neighbours[i]] *= 0.5;
                }
                ++count;
                if(count > maxCount) {
                    std::cerr << "ERROR: Max number of iterations is reached! Function: scaleVertices() failed!\n" << std::endl;
                    return 1;
                }
            }

            return 0;
        }

        // Check if no mesh vertex around the corner creates a foldover.
        bool foldover(const std::vector<int> &firstRing)
        {
            const int nSizeFirst = firstRing.size();

            assert(nSizeFirst != 0);

            // Go through all the neighbours of the corner vertex
            // and find 1-ring neighbourhood for each of them.
            std::vector<int> secondRing;
            for(int i = 0; i < nSizeFirst; ++i) {
                getRing(firstRing[i], secondRing);

                const int nSizeSecond = secondRing.size();
                int size              = nSizeSecond;

                assert(nSizeSecond != 0);

                if(i == 0 || i == nSizeFirst - 1) size--; // for boundary vertices the size is 1 element smaller
                const Vertex &v2 = vertex[firstRing[i]];

                assert(size != 0);

                // Find determinant for each triangle and check if it is negative/positive.
                for(int j = 0; j < size; ++j) {
                    const Vertex &v3 = vertex[secondRing[j]];
                    const Vertex &v1 = vertex[secondRing[(j + 1) % nSizeSecond]];

                    const double x1 = v1.p().x; const double y1 = v1.p().y;
                    const double x2 = v2.p().x; const double y2 = v2.p().y;
                    const double x3 = v3.p().x; const double y3 = v3.p().y;

                    // Compute the following determinant:
                    //
                    //       | x1 y1 1 |
                    // det = | x2 y2 1 |
                    //       | x3 y3 1 |
                    //
                    // if det > 0 or det == 0, no interior fold-over is detected;
                    // if det < 0, an interior fold-over happend;

                    const double det = x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2);

                    // std::cout << "det = " << det << "\n";

                    if(det < 0.0) return true; // fold-over found
                }
            }

            return false; // no fold-over is detected
        }

        // Fix barycentric coordinates around the corner.
        int fixCoordinates(const std::vector<int> &neighbours, const std::vector<Vertex> &prevVertex, const std::vector<Face> &prevFace)
        {
            const double eps = epsilon();
            const int nSize  = neighbours.size();
            const int numF   = prevFace.size();
            const int numC   = vertex[0].b().size();

            // Go through all the neighbours.
            for(int i = 0; i < nSize; ++i) {

                // Go through all the faces of the mesh before a ternary step.
                int count = 0;
                for(int j = 0; j < numF; ++j) {

                    // Compute triangle coordinates.
                    std::vector<double> lambda;
                    computeTriangleCoordinates(prevVertex[prevFace[j].v[0]], prevVertex[prevFace[j].v[1]], prevVertex[prevFace[j].v[2]], neighbours[i], lambda);

                    // If we found a face to which current neighbour belongs,
                    // linearly sample barycentric cordinates from the vertices of this face.
                    if(lambda[0] > -eps && lambda[1] > -eps && lambda[2] > -eps) {
                        for(int k = 0; k < numC; ++k) {
                            vertex[neighbours[i]].b()[k] = std::fabs(lambda[0]) * prevVertex[prevFace[j].v[0]].b()[k] +
                                                           std::fabs(lambda[1]) * prevVertex[prevFace[j].v[1]].b()[k] +
                                                           std::fabs(lambda[2]) * prevVertex[prevFace[j].v[2]].b()[k] ;
                        }

                        vertex[neighbours[i]].pm() = std::fabs(lambda[0]) * prevVertex[prevFace[j].v[0]].pm() +
                                                     std::fabs(lambda[1]) * prevVertex[prevFace[j].v[1]].pm() +
                                                     std::fabs(lambda[2]) * prevVertex[prevFace[j].v[2]].pm() ;

                        break;
                    } else ++count;
                }

                // If we did not find at least one face to which current neighbour belongs, return an error.
                if(count == numF) {
                    std::cerr << "ERROR: No appropriate face was found! fixCoordinates() function failed!\n" << std::endl;
                    return 1;
                }
            }

            return 0;
        }

        // Compute triangle coordinates.
        inline void computeTriangleCoordinates(const Vertex &v1, const Vertex &v2, const Vertex &v3, const int meshInd, std::vector<double> &coordinates) const
        {
            coordinates.resize(3);

            // Compute some related sub-areas.
            const double area_second = 0.5 * crossProduct(v2, v3, vertex[meshInd]);
            const double area_third  = 0.5 * crossProduct(v3, v1, vertex[meshInd]);

            // Compute the total inverted area of the triangle.
            const double total_area          = 0.5 * crossProduct(v1, v2, v3);
            const double inverted_total_area = 1.0 / total_area;

            // Compute all three coordinate functions.
            coordinates[0] = area_second * inverted_total_area;
            coordinates[1] = area_third  * inverted_total_area;
            coordinates[2] = 1.0 - coordinates[0] - coordinates[1];
        }

        // Compute cross product.
        inline double crossProduct(const Vertex &v1, const Vertex &v2, const Vertex &v3) const
        {
            const double x1 = v1.p().x; const double y1 = v1.p().y;
            const double x2 = v2.p().x; const double y2 = v2.p().y;
            const double x3 = v3.p().x; const double y3 = v3.p().y;

            return x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2);
        }

        // Subdivide mesh with the modified Loop scheme.
        void subdivideMesh()
        {
            ///////////////////////////////////
            // Do ODD step - add new points. //
            ///////////////////////////////////

            // Initiate some data.
            const int numV  = numVertices();
            const int numE  = numEdges();
            const int numHE = numHalfedges();

            vertex.resize(numV + numE);
            std::vector<int> neighbours(4);

            int idxV = numV;

            // Loop over all halfedges.
            for(int i = 0; i < numHE; ++i) {

                // Consider only the "big brother" of each halfedge-pair.
                const int i_ = halfedge[i].neigh;
                if(i > i_) {

                    // Local configuration:
                    /*
                             w1
                            /||\
                          /  ||  \
                        /   ^||    \
                      w2     []     w4
                        \   *||    /
                          \  ||  /
                            \||/
                             w3
                    */
                    // The current halfedge and its direction are indicated by "*" and "^".
                    // Create an odd vertex at the position "[]".
                    // NOTE: the right triangle and w4 do not exist if "*" is at the boundary.

                    neighbours[0] = halfedge[i].dest;
                    neighbours[1] = halfedge[halfedge[i].next].dest;
                    neighbours[2] = halfedge[halfedge[i].prev].dest;

                    // Mask for odd interior vertices.
                    if(i_ >= 0) {
                        neighbours[3] = halfedge[halfedge[i_].next].dest;
                        createOddIntVertex(neighbours, vertex[idxV]);
                    }
                    // Mask for odd boundary vertices.
                    else {
                        createOddBdyVertex(neighbours, vertex[idxV]);
                    }
                    idxV++;

                } // if a halfedge is the "big brother"

            } // for all halfedges

            // Update halfedge connectivity.
            updateHEConnectivity(numV);

            /////////////////////////////////////
            // Do EVEN step - move old points. //
            /////////////////////////////////////

            for(int i = 0; i < numV; ++i) {
                getRing(i, neighbours);

                // Mask for even interior vertices.
                if(vertex[i].type == 3) {
                    moveEvenIntVertex(neighbours, vertex[i]);
                }
                // Mask for even boundary vertices.
                else {
                    moveEvenBdyVertex(neighbours, vertex[i]);
                }
            }
        }

        // Update halfedge connectivity in the mesh after each subdivision step.
        void updateHEConnectivity(const int numV)
        {
            const int numHE = numHalfedges();
            const int numF  = numHE / 3;

            halfedge.resize(4 * numHE);

            // Update an outgoing halfedge of each even vertex.
            for(int i = 0; i < numV; ++i) {
                const int k = vertex[i].out;

                // An outgoing halfedge k is the local edge e: k % 3 in the triangle t: k / 3.
                // The new halfedges of that triangle have indices starting at...
                const int base = 3 * (numF + 3 * (k / 3));

                // In general, an offset of the new halfedge e: (i1,j1) is 3 * j1 + i1.
                // Here, we want e: (j,j+1) with j = k % 3.
                vertex[i].out = base + 3 * ((k + 1) % 3) + k % 3;
            }

            // Fix the connectivity between halfedges that replace the "old" edges,
            // set their destination, and set an outgoing halfedge for each odd vertex.

            int idxV = numV;

            // Loop over all halfedges.
            for(int i = 0; i < numHE; ++i) {

                // Consider only the "big brother" of each halfedge-pair.
                const int i_ = halfedge[i].neigh;
                if(i > i_) {

                    // Get indices of halfedges that will replace the current halfedge.
                    // Those are e: (k,k+1) and e: (k,k+2) with k = i % 3 (as above).
                    int base = 3 * (numF + 3 * (i / 3));
                    const int i0 = base + 3 * ((i + 1) % 3) + i % 3;
                    const int i1 = base + 3 * ((i + 2) % 3) + i % 3;

                    // Set their destination.
                    halfedge[i0].dest = idxV;
                    halfedge[i1].dest = halfedge[i].dest;

                    // Set an outgoing halfedge for odd vertex on the current halfedge.
                    vertex[idxV].out = i1;

                    // Is it an interior halfedge?
                    if(i_ >= 0) {
                        // Get indices of halfedges that will replace the neighbouring halfedge.
                        base = 3 * (numF + 3 * (i_ / 3));
                        const int i0_ = base + 3 * ((i_ + 1) % 3) + i_ % 3;
                        const int i1_ = base + 3 * ((i_ + 2) % 3) + i_ % 3;

                        // Mutually connect new halfedges.
                        halfedge[i0].neigh = i1_;
                        halfedge[i1].neigh = i0_;

                        halfedge[i0_].neigh = i1;
                        halfedge[i1_].neigh = i0;

                        // Also set destination of new neighbouring halfedges.
                        halfedge[i0_].dest = idxV;
                        halfedge[i1_].dest = halfedge[i_].dest;
                    }
                    else {
                        // If the current halfedge is at the boundary, then their children are, too.
                        halfedge[i0].neigh = -1;
                        halfedge[i1].neigh = -1;
                    }
                    idxV++;

                } // if a halfedge is the "big brother"

            } // for all halfedges

            // Fix connectivity and destination of new halfedges "inside" old triangles.
            int Ej = 0;
            int Ejj[3], Ejmj[3];
            for(int i = 0; i < numF; ++i) {

                // Explicitely set indices of new halfedges e: (0,1), e: (1,2), e: (2,0).
                const int base = 3 * (numF + 3 * i);
                Ejmj[1] = base + 3;
                Ejmj[2] = base + 7;
                Ejmj[0] = base + 2;

                // And e: (0,0), e: (1,1), e: (2,2).
                Ejj[0] = base;
                Ejj[1] = base + 4;
                Ejj[2] = base + 8;

                for(int j = 0; j < 3; ++j, Ej++) {
                    // Connect new halfedges with each other.
                    halfedge[Ej].neigh     = Ejj[j];
                    halfedge[Ejj[j]].neigh = Ej;

                    // Set their destination.
                    idxV = halfedge[Ejmj[j]].dest;
                    halfedge[Ej].dest               = idxV;
                    halfedge[Ejj[(j + 1) % 3]].dest = idxV;
                }
            }

            // Set next and previous halfedges for all new halfedges.
            int idxE = numHE;
            for(int i = numF; i < 4 * numF; ++i) {
                const int base = idxE;
                for(int j = 0; j < 3; ++j, idxE++) {
                    halfedge[idxE].prev = base + (j + 2) % 3;
                    halfedge[idxE].next = base + (j + 1) % 3;
                }
            }
        }

        // Create an interior odd vertex.
        inline void createOddIntVertex(const std::vector<int> &neighbours, Vertex &newV) const
        {
            assert(neighbours.size() == 4);

            double c[4] = { 0.375, 0.125, 0.375, 0.125 };

            // Zorin modification.
            if(vertex[neighbours[0]].type == 0 || vertex[neighbours[0]].type == 1) {

                const int k        = vertex[neighbours[0]].val - 1;
                const double alpha = vertex[neighbours[0]].alpha;
                const double theta = (2.0 * M_PI - alpha) / k;
                double gamma       = 0.5 - 0.25 * cos(theta);

                // If both vertices of an edge are convex or concave corners,
                // we return to the original weights 3/8.
                if(vertex[neighbours[2]].type == 0 || vertex[neighbours[2]].type == 1) gamma = 0.375;

                // Both weights must behave as 1/2.
                c[0] = 0.75 - gamma;
                c[2] = gamma;
            }

            // Zorin modification.
            if(vertex[neighbours[2]].type == 0 || vertex[neighbours[2]].type == 1) {

                const int k        = vertex[neighbours[2]].val - 1;
                const double alpha = vertex[neighbours[2]].alpha;
                const double theta = (2.0 * M_PI - alpha) / k;
                double gamma       = 0.5 - 0.25 * cos(theta);

                // If both vertices of an edge are convex or concave corners,
                // we return to the original weights 3/8.
                if(vertex[neighbours[0]].type == 0 || vertex[neighbours[0]].type == 1) gamma = 0.375;

                // Both weights must behave as 1/2.
                c[2] = 0.75 - gamma;
                c[0] = gamma;
            }

            // Create a new interior vertex.
            newV  = vertex[neighbours[0]];
            newV *= c[0];
            newV += c[1] * vertex[neighbours[1]];
            newV += c[2] * vertex[neighbours[2]];
            newV += c[3] * vertex[neighbours[3]];

            // Update its flags.
            newV.val   = 6;
            newV.out   = -1;
            newV.type  = 3;
            newV.alpha = -1.0;
        }

        // Create a boundary odd vertex.
        inline void createOddBdyVertex(const std::vector<int> &neighbours, Vertex &newV) const
        {
            assert(neighbours.size() > 0);

            // Create a new boundary vertex.
            newV  = vertex[neighbours[0]];
            newV += vertex[neighbours[2]];
            newV *= 0.5;

            // Update its flags.
            newV.val   = 4;
            newV.out   = -1;
            newV.type  = 2;
            newV.alpha = -1.0;
        }

        // Move an interior even vertex.
        inline void moveEvenIntVertex(const std::vector<int> &neighbours, Vertex &center) const
        {
            const int n = neighbours.size();

            assert(n > 0);
            assert(center.alpha == -1.0);
            assert(center.type == 3);

            // Compute weights for arbitrary valency.
            double c0  = 0.375 + 0.25 * cos(2.0 * M_PI / n);
                   c0 *= c0 * 1.6;

            // Compute the barycenter of a set of vertices.
            Vertex barycenter = vertex[neighbours[0]];
            for(int i = 1; i < n; ++i) barycenter += vertex[neighbours[i]];

            const double nInv = 1.0 / n;

            // Move an old interior vertex.
            center *= c0;
            center += (nInv * (1.0 - c0)) * barycenter;
        }

        // Move a boundary even vertex.
        inline void moveEvenBdyVertex(const std::vector<int> &neighbours, Vertex &center) const
        {
            assert(center.type != 3);

            if(center.type == 2) { // FLAT vertex

                assert(center.alpha == -1.0);
                assert(neighbours.size() > 1);

                const int n = center.val;

                assert(n > 1);
                assert(int(neighbours.size()) == n);

                center *= 2.0;
                center += vertex[neighbours[0]];
                center += vertex[neighbours[n-1]];
                center *= 0.25;
            }

            // CONVEX AND CONCAVE vertices are interpolated!
        }

        // Subdivide the morphed mesh.
        void subdivideMorphedMesh(const std::vector<Halfedge> &prevHalfedges,
                                  const std::vector<Halfedge> &nextHalfedges,
                                  const std::vector<Vertex>    &prevVertices,
                                        std::vector<Vertex>    &nextVertices)
        {
            ///////////////////////////////////
            // Do ODD step - add new points. //
            ///////////////////////////////////

            const int numV  = prevVertices.size();
                  int numHE = prevHalfedges.size();

            std::vector<int> neighbours(4);
            int idxV = numV;

            // Loop over all halfedges.
            for(int i = 0; i < numHE; ++i) {

                // Consider only the "big brother" of each halfedge-pair.
                const int i_ = prevHalfedges[i].neigh;
                if(i > i_) {

                    // Local configuration:
                    /*
                             w1
                            /||\
                          /  ||  \
                        /   ^||    \
                      w2     []     w4
                        \   *||    /
                          \  ||  /
                            \||/
                             w3
                    */
                    // The current halfedge and its direction are indicated by "*" and "^".
                    // Create an odd vertex at the position "[]".
                    // NOTE: the right triangle and w4 do not exist if "*" is at the boundary.

                    neighbours[0] = prevHalfedges[i].dest;
                    neighbours[1] = prevHalfedges[prevHalfedges[i].next].dest;
                    neighbours[2] = prevHalfedges[prevHalfedges[i].prev].dest;

                    // Mask for odd interior vertices.
                    if(i_ >= 0) {
                        neighbours[3] = prevHalfedges[prevHalfedges[i_].next].dest;
                        createOddIntMorphedVertex(neighbours, prevVertices, nextVertices[idxV]);
                    }
                    // Mask for odd boundary vertices.
                    else {
                        createOddBdyMorphedVertex(neighbours, prevVertices, nextVertices[idxV]);
                    }
                    idxV++;

                } // if a halfedge is the "big brother"

            } // for all halfedges

            /////////////////////////////////////
            // Do EVEN step - move old points. //
            /////////////////////////////////////

            numHE = nextHalfedges.size();

            for(int i = 0; i < numV; ++i) {

                int n = 0, prev, curr, next;

                // Find neighbours of the vertex.
                curr = next = nextVertices[i].out;

                assert(nextVertices[i].val > 0);
                neighbours.resize(nextVertices[i].val);
                do {
                    neighbours[n] = nextHalfedges[curr].dest;

                    prev = nextHalfedges[curr].prev;
                    curr = nextHalfedges[prev].neigh;
                    n++;
                } while((curr >= 0) && (curr != next));

                // Add one more neighbour for boundary vertices.
                if(curr < 0) {
                    curr          = nextHalfedges[prev].prev;
                    neighbours[n] = nextHalfedges[curr].dest;
                }

                // Mask for even interior vertices.
                if(nextVertices[i].type == 3) {
                    moveEvenIntMorphedVertex(neighbours, nextVertices, prevVertices[i], nextVertices[i]);
                }
                // Mask for even boundary vertices.
                else {
                    moveEvenBdyMorphedVertex(neighbours, nextVertices, prevVertices[i], nextVertices[i]);
                }
            }
        }

        // Create an interior odd morphed vertex.
        inline void createOddIntMorphedVertex(const std::vector<int> &neighbours, const std::vector<Vertex> &prevVertex, Vertex &newV) const
        {
            assert(neighbours.size() == 4);

            double c[4] = { 0.375, 0.125, 0.375, 0.125 };

            // Zorin modification.
            if(prevVertex[neighbours[0]].type == 0 || prevVertex[neighbours[0]].type == 1) {

                const int k        = prevVertex[neighbours[0]].val - 1;
                const double alpha = prevVertex[neighbours[0]].alpha;
                const double theta = (2.0 * M_PI - alpha) / k;
                double gamma       = 0.5 - 0.25 * cos(theta);

                // If both vertices of an edge are convex or concave corners,
                // we return to the original weights 3/8.
                if(prevVertex[neighbours[2]].type == 0 || prevVertex[neighbours[2]].type == 1) gamma = 0.375;

                // Both weights must behave as 1/2.
                c[0] = 0.75 - gamma;
                c[2] = gamma;
            }

            // Zorin modification.
            if(prevVertex[neighbours[2]].type == 0 || prevVertex[neighbours[2]].type == 1) {

                const int k        = prevVertex[neighbours[2]].val - 1;
                const double alpha = prevVertex[neighbours[2]].alpha;
                const double theta = (2.0 * M_PI - alpha) / k;
                double gamma       = 0.5 - 0.25 * cos(theta);

                // If both vertices of an edge are convex or concave corners,
                // we return to the original weights 3/8.
                if(prevVertex[neighbours[0]].type == 0 || prevVertex[neighbours[0]].type == 1) gamma = 0.375;

                // Both weights must behave as 1/2.
                c[2] = 0.75 - gamma;
                c[0] = gamma;
            }

            // Create a new interior vertex.
            newV.pm()  = prevVertex[neighbours[0]].pm();
            newV.pm() *= c[0];
            newV.pm() += c[1] * prevVertex[neighbours[1]].pm();
            newV.pm() += c[2] * prevVertex[neighbours[2]].pm();
            newV.pm() += c[3] * prevVertex[neighbours[3]].pm();
        }

        // Create a boundary odd morphed vertex.
        inline void createOddBdyMorphedVertex(const std::vector<int> &neighbours, const std::vector<Vertex> &prevVertex, Vertex &newV) const
        {
            assert(neighbours.size() > 0);

            // Create a new boundary vertex.
            newV.pm()  = prevVertex[neighbours[0]].pm();
            newV.pm() += prevVertex[neighbours[2]].pm();
            newV.pm() *= 0.5;
        }

        // Move an interior even morphed vertex.
        inline void moveEvenIntMorphedVertex(const std::vector<int> &neighbours, const std::vector<Vertex> &nextVertex, const Vertex &center, Vertex &newV) const
        {
            const int n = neighbours.size();

            assert(n > 0);
            assert(center.alpha == -1.0);
            assert(center.type == 3);

            // Compute weights for arbitrary valency.
            double c0  = 0.375 + 0.25 * cos(2.0 * M_PI / n);
                   c0 *= c0 * 1.6;

            // Compute the barycenter of a set of vertices.
            Vector2d barycenter;
            for(int i = 0; i < n; ++i) barycenter += nextVertex[neighbours[i]].pm();

            const double nInv = 1.0 / n;

            // Move an old interior vertex.
            newV.pm()  = center.pm();
            newV.pm() *= c0;
            newV.pm() += (nInv * (1.0 - c0)) * barycenter;
        }

        // Move a boundary even morphed vertex.
        inline void moveEvenBdyMorphedVertex(const std::vector<int> &neighbours, const std::vector<Vertex> &nextVertex, const Vertex &center, Vertex &newV) const
        {
            assert(center.type != 3);

            if(center.type == 2) { // FLAT vertex

                assert(center.alpha == -1.0);
                assert(neighbours.size() > 1);

                const int n = center.val;

                assert(n > 1);
                assert(int(neighbours.size()) == n);

                newV.pm()  = center.pm();
                newV.pm() *= 2.0;
                newV.pm() += nextVertex[neighbours[0]].pm();
                newV.pm() += nextVertex[neighbours[n-1]].pm();
                newV.pm() *= 0.25;
            } else { // CONVEX AND CONCAVE vertices are interpolated!
                newV.pm() = center.pm();
            }
        }
    };

} // namespace warpit

#endif // TRIMESH_HPP
