//* Author: Dmitry Anisimov (c) 2016-2017
//* With any questions please contact: danston@ymail.com

#ifndef VERTEX_HPP
#define VERTEX_HPP

// STL includes.
#include <vector>

// Local includes.
#include "vector2d.hpp"

// warpit namespace
namespace warpit {

    // Vertex data structure.
    class Vertex
    {
    public:
        // Constructors.
        Vertex() : out(-1), val(0), type(3), alpha(-1.0), _p(), _pm(), _b() { }

        Vertex(const std::size_t bSize) : out(-1), val(0), type(3), alpha(-1.0), _p(), _pm(), _b(bSize, 0.0) { }

        Vertex(const Vertex &V) : out(V.out), val(V.val), type(V.type), alpha(V.alpha), _p(V.p()), _pm(V.pm()), _b(V.b()) { }

        Vertex(const Vector2d &pos) : out(-1), val(0), type(3), alpha(-1.0), _p(pos), _pm(), _b() { }

        // Flags.
        int    out;   // index of the outgoing edge
        int    val;   // valency of the vertex
        int    type;  // type of the vertex: 0 - CONVEX; 1 - CONCAVE; 2 - FLAT; 3 - INTERIOR
        double alpha; // external angle for each corner

        // Return position.
        inline Vector2d& p()
        {
            return _p;
        }

        // Return const position.
        inline const Vector2d& p() const
        {
            return _p;
        }

        // Return morphed position.
        inline Vector2d& pm()
        {
            return _pm;
        }

        // Return const morphed position.
        inline const Vector2d& pm() const
        {
            return _pm;
        }

        // Return barycentric coordinates attached to the vertex.
        inline std::vector<double>& b()
        {
            return _b;
        }

        // Return const barycentric coordinates attached to the vertex.
        inline const std::vector<double>& b() const
        {
            return _b;
        }

        // Overload some basic operators.

        // Multiplication by a constant from right.
        Vertex operator*(const double scalar) const
        {
            const std::size_t bSize = _b.size();
            Vertex newV(bSize);

            newV.p()  = scalar * _p;
            newV.pm() = scalar * _pm;

            assert(val >= 0);

            newV.out = out;
            newV.val = val;

            newV.alpha = alpha;

            assert(type >= 0 && type <= 3);

            newV.type = type;
            for(std::size_t i = 0; i < bSize; ++i) newV.b()[i] = scalar * _b[i];

            return newV;
        }

        // Multiplication by a constant from left.
        friend inline Vertex operator*(const double scalar, const Vertex &V)
        {
            return V * scalar;
        }

        // Addition of two vertices.
        Vertex operator+(const Vertex &V) const
        {
            const std::size_t bSize = _b.size();
            Vertex newV(bSize);

            newV.p()  = _p  + V.p();
            newV.pm() = _pm + V.pm();

            for(std::size_t i = 0; i < bSize; ++i) newV.b()[i] = _b[i] + V.b()[i];

            return newV;
        }

        // Addition of two vertices without creating a new vertex.
        void operator+=(const Vertex &V)
        {
            const std::size_t bSize = _b.size();

            _p  += V.p();
            _pm += V.pm();

            for(std::size_t i = 0; i < bSize; ++i) _b[i] += V.b()[i];
        }

        // Multiplication by a constant from right without creating a new vertex.
        void operator*=(const double scalar)
        {
            const std::size_t bSize = _b.size();

            _p  *= scalar;
            _pm *= scalar;

            for(std::size_t i = 0; i < bSize; ++i) _b[i] *= scalar;
        }

        // Overload the == operator.
        inline bool operator==(const Vertex &V) const
        {
            return V.p() == _p;
        }

    private:
        // Elements.
        Vector2d _p;            // position of the vertex
        Vector2d _pm;           // position of the morphed vertex
        std::vector<double> _b; // barycentric coordinates attached to each vertex
    };

} // namespace warpit

#endif // VERTEX_HPP
