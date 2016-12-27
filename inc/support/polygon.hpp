//* Author: Dmitry Anisimov (c) 2016-2017
//* With any questions please contact: danston@ymail.com

// Todo list:
// 1. Add segment barycentric coordinates and improve the quality of some functions based on these coordinates, for example the function projectOnTheClosestEdge().
// 2. Improve the self-intersected function.

#ifndef POLYGON_HPP
#define POLYGON_HPP

// STL includes.
#include <fstream>

// Local includes.
#include "vertex.hpp"

// warpit namespace
namespace warpit {

    // The polygon class.
    class Polygon
    {
    public:
        // The default constructor.
        Polygon() : vertex(), eps(1.0e-14) { }

        // Construct a polygon from the set of vertices.
        Polygon(const std::vector<Vertex> &initialVertices) : vertex(initialVertices), eps(1.0e-14) { }

        // Construct a rectangular polygon from its width and height.
        Polygon(const double width, const double height) : vertex(4), eps(1.0e-14)
        {
            vertex[0] = Vertex(Vector2d(0.0, 0.0));
            vertex[1] = Vertex(Vector2d(width, 0.0));
            vertex[2] = Vertex(Vector2d(width, height));
            vertex[3] = Vertex(Vector2d(0.0, height));
        }

        // Check if the polygon contains any vertices.
        inline bool empty() const
        {
            if(numVertices() == 0) return true;
            return false;
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

        // Return vertices of the polygon.
        inline std::vector<Vertex>& vertices()
        {
            return vertex;
        }

        // Return const vertices of the polygon.
        inline const std::vector<Vertex>& vertices() const
        {
            return vertex;
        }

        // Return the number of vertices in the polygon.
        inline int numVertices() const
        {
            return vertex.size();
        }

        // Clear the polygon.
        inline void clear()
        {
            vertex.clear();
        }

        // Check if two vertices create an edge of the polygon.
        inline bool isEdge(const Vertex &v1, const Vertex &v2) const
        {
            const int numV = numVertices();

            for(int i = 0; i < numV; ++i)
                if(vertex[i].p() == v1.p() && vertex[(i + 1) % numV].p() == v2.p()) return true;

            return false;
        }

        // Check if the vertex is on some edge of the polygon.
        inline bool isOnEdge(const int i, const Vertex &v) const
        {
            const int numV = numVertices();

            const int    iPlus  = (i + 1) % numV;
            const double length = (vertex[i].p() - vertex[iPlus].p()).length();
            const double actualDistance = vertexLineDistance(vertex[i], vertex[iPlus], v);

            const double distance1 = (vertex[i].p() - v.p()).length();
            const double distance2 = (vertex[iPlus].p() - v.p()).length();

            return actualDistance < eps && distance1 <= length && distance2 <= length;
        }

        // Check if the vertex belongs to the polygon's boundary.
        // If it does, return index of the corresponding edge.
        inline int findEdgeIndex(const Vertex &v) const
        {
            const int numV = numVertices();

            for(int i = 0; i < numV; ++i)
                if(isOnEdge(i, v)) return i;

            return -1;
        }

        // Insert a new vertex in the polygon.
        inline void insert(const Vertex &v, const int index)
        {
            const int numV = numVertices();
            const int realIndex = (index + 1) % numV;
            vertex.insert(vertex.begin() + realIndex, v);
        }

        // Compute the barycenter of the polygon.
        inline Vector2d barycenter() const
        {
            const int numV = numVertices();

            Vector2d bary;
            for(int i = 0; i < numV; ++i) bary += vertex[i].p();
            bary *= 1.0 / numV;

            return bary;
        }

        // Move the polygon to the provided barycenter.
        inline void setBarycenter(const Vector2d &bary = Vector2d())
        {
            const int numV = numVertices();

            const Vector2d currentBary = barycenter();
            const Vector2d delta = bary - currentBary;

            for(int i = 0; i < numV; ++i) vertex[i].p().translate(delta);
        }

        // Fit the size of the polygon to a given size.
        inline void fitToSize(const double size = 1.0)
        {
            double minX = std::numeric_limits<double>::max(), maxX = -std::numeric_limits<double>::max();
            double minY = std::numeric_limits<double>::max(), maxY = -std::numeric_limits<double>::max();

            const int numV = numVertices();

            for(int i = 0; i < numV; ++i) {
                minX = std::min(vertex[i].p().x, minX);
                minY = std::min(vertex[i].p().y, minY);

                maxX = std::max(vertex[i].p().x, maxX);
                maxY = std::max(vertex[i].p().y, maxY);
            }

            const double diameter = (Vector2d(maxX, maxY) - Vector2d(minX, minY)).length();
            for(int i = 0; i < numV; ++i) vertex[i].p() *= (0.5 * size) / diameter;
        }

        // Check if the vertex belongs to the boundary of the polygon.
        bool belongsToBoundary(const Vertex &v) const
        {
            const int numV = numVertices();

            for(int i = 0; i < numV; ++i)
                if(isOnEdge(i, v)) return true;

            return false;
        }

        // Check if the vertex belongs to this polygon.
        bool contains(const Vertex &v) const
        {
            const double x = v.p().x;
            const double y = v.p().y;

            const int numV = numVertices();

            int j = numV - 1;
            bool contained = false;

            for(int i = 0; i < numV; i++) {
                const double yi = vertex[i].p().y;
                const double yj = vertex[j].p().y;
                const double xi = vertex[i].p().x;
                const double xj = vertex[j].p().x;

                if((vertex[i].p() - v.p()).length() < eps) return true;

                if(((yi < y && yj >= y) || (yj < y && yi >= y)) && (xi <= x || xj <= x)) {
                    const double xL = xi + (y - yi)  / (yj - yi) * (xj - xi);
                    contained ^= xL < x;

                    if(std::fabs(x - xL) < eps) return true;
                }

                j = i;
            }

            return contained;
        }

        // Find the closest to the given vertex edge of the polygon.
        int findClosestEdge(Vertex &v) const
        {
            const int numV = numVertices();

            double minDistance = std::numeric_limits<double>::max();
            int index = -1;

            for(int i = 0; i < numV; ++i) {
                const Vertex &v1 = vertex[i];
                const Vertex &v2 = vertex[(i + 1) % numV];

                const double distance = vertexLineDistance(v1, v2, v);
                if(distance < minDistance) {
                    index = i;
                    minDistance = distance;
                }
            }

            return index;
        }

        // Find the closest vertex.
        int findClosestVertex(Vertex &v) const
        {
            const int numV = numVertices();

            for(int i = 0; i < numV; ++i)
                if(v.p() == vertex[i].p()) return i;

            return -1;
        }

        // Project the vertex on the closest edge of the polygon.
        void projectVertexOnClosestEdge(Vertex &v) const
        {
            const int numV  = numVertices();
            const int index = findClosestEdge(v);

            const Vector2d &p1 = vertex[index].p();
            const Vector2d &p2 = vertex[(index + 1) % numV].p();

            Vector2d direction = p2 - p1;
            const double length = direction.length();
            direction *= 1.0 / length;
            const Vector2d tmp = v.p() - p1;

            v.p() = tmp.scalarProduct(direction) * direction + p1;

            if(!contains(v) && !belongsToBoundary(v)) v.p() = 0.5 * p1 + 0.5 * p2;
        }

        // Find the bounding box of the polygon.
        void boundingBox(Vector2d &blCorner, Vector2d &urCorner) const
        {
            const int numV = numVertices();

            double minX = std::numeric_limits<double>::max(), maxX = -std::numeric_limits<double>::max();
            double minY = std::numeric_limits<double>::max(), maxY = -std::numeric_limits<double>::max();

            for(int i = 0; i < numV; ++i) {
                minX = std::min(minX, vertex[i].p().x);
                maxX = std::max(maxX, vertex[i].p().x);
                minY = std::min(minY, vertex[i].p().y);
                maxY = std::max(maxY, vertex[i].p().y);
            }

            blCorner = Vector2d(minX, minY);
            urCorner = Vector2d(maxX, maxY);
        }

        // Check if the polygon is self-intersected. Slow O(n^2) algorithm.
        bool selfIntersected() const
        {
            const int numV = numVertices();
            Vector2d interPoint;

            for(int i = 0; i < numV; ++i) {
                for(int j = 0; j < i; ++j) {
                    if(areIntersected(vertex[i].p(), vertex[(i + 1) % numV].p(), vertex[j].p(), vertex[(j + 1) % numV].p(), interPoint)) {
                        const double d1 = (vertex[i].p() - interPoint).length();
                        const double d2 = (vertex[(i + 1) % numV].p() - interPoint).length();
                        const double d3 = (vertex[j].p() - interPoint).length();
                        const double d4 = (vertex[(j + 1) % numV].p() - interPoint).length();

                        if(d1 > eps && d2 > eps && d3 > eps && d4 > eps) return true;
                    }
                }

                for(int j = i + 1; j < numV; ++j) {
                    if(areIntersected(vertex[i].p(), vertex[(i + 1) % numV].p(), vertex[j].p(), vertex[(j + 1) % numV].p(), interPoint)) {
                        const double d1 = (vertex[i].p() - interPoint).length();
                        const double d2 = (vertex[(i + 1) % numV].p() - interPoint).length();
                        const double d3 = (vertex[j].p() - interPoint).length();
                        const double d4 = (vertex[(j + 1) % numV].p() - interPoint).length();

                        if(d1 > eps && d2 > eps && d3 > eps && d4 > eps) return true;
                    }
                }
            }

            return false;
        }

        // Load polygon from a file.
        void load(const std::string &fileName)
        {
            clear();

            if(fileName.empty()) return;

            std::ifstream readFile(fileName.c_str(), std::ios_base::in);

            if(!readFile) {
                std::cout << "\nError reading .poly file!" << std::endl;
                exit(EXIT_FAILURE);
            }

            std::string line;
            getline(readFile, line);
            if(line != "POLY") {
                std::cout << "\nThis is not a valid POLY file!" << std::endl;
                exit(EXIT_FAILURE);
            }

            int numV, dummyFirst;
            double dummyLast;

            readFile >> numV;

            Vertex newVertex;
            for(int i = 0; i < numV; ++i) {
                readFile >> dummyFirst >> newVertex.p().x >> newVertex.p().y >> dummyLast;
                vertex.push_back(newVertex);
            }

            readFile.close();

            fitToSize(1.0);
            setBarycenter();
        }

        // Save polygon as a .poly file.
        void save(const std::string &fileName)
        {
            if(fileName.empty()) return;

            std::ofstream saveFile(fileName.c_str(), std::ios_base::out);

            if(!saveFile) {
                std::cout << "\nError saving .poly file!" << std::endl;
                exit(EXIT_FAILURE);
            }

            const int numV = numVertices();

            saveFile << "POLY" << std::endl;
            saveFile <<  numV  << std::endl;

            for(int i = 0; i < numV; ++i)
                saveFile << i + 1 << " " << vertex[i].p().x << " "
                                         << vertex[i].p().y << " "
                                         << 0.0 << std::endl;

            saveFile.close();
        }

    private:
        // Vertices of the polygon.
        std::vector<Vertex> vertex;

        // Tolerance.
        double eps;

        // Check if two segments are intersected within each other.
        inline bool areIntersected(const Vector2d &v1, const Vector2d &v2, const Vector2d &v3, const Vector2d &v4, Vector2d &v) const
        {
            const double lambda1 = intersectionCoefficient(v1, v2, v3, v4);
            const double lambda2 = intersectionCoefficient(v3, v4, v1, v2);

            const Vector2d u = v2 - v1;
            v = v1 + (u * lambda1);

            return lambda1 >= -eps && lambda1 <= 1.0 + eps && lambda2 >= -eps && lambda2 <= 1.0 + eps;
        }

        // Intersection coefficient.
        inline double intersectionCoefficient(const Vector2d &v1, const Vector2d &v2, const Vector2d &v3, const Vector2d &v4) const
        {
            const Vector2d u = v2 - v1;
            const Vector2d v = v4 - v3;
            const Vector2d w = v4 - v1;

            return -(-w[0] * v[1] + w[1] * v[0]) / (u[0] * v[1] - u[1] * v[0]);
        }

        // The distance from the vertex to a line.
        static inline double vertexLineDistance(const Vertex &v1, const Vertex &v2, const Vertex &v)
        {
            Vector2d direction = v2.p() - v1.p();
            const double length = direction.length();
            direction *= 1.0 / length;

            return (v1.p() - v.p() - (v1.p() - v.p()).scalarProduct(direction) * direction).length();
        }
    };

} // namespace warpit

#endif // POLYGON_HPP
