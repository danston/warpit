//* Authors: Dmitry Anisimov
//* halfedge.hpp
//* With any questions please contact: danston@ymail.com

#ifndef HALFEDGE_HPP
#define HALFEDGE_HPP

// warpit namespace.
namespace warpit {

    // Halfedge data structure.
    class Halfedge
    {
    public:
        // Constructor.
        Halfedge() : prev(-1), next(-1), neigh(-1), dest(-1) { }

        // Flags.
        int prev;  // index of the previous halfedge
        int next;  // index of the next halfedge
        int neigh; // the neighboring halfedge
        int dest;  // index of the vertex at the end of the halfedge
    };

} // namespace warpit

#endif // HALFEDGE_HPP
