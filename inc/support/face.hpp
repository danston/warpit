//* Author: Dmitry Anisimov (c) 2016-2017
//* With any questions please contact: danston@ymail.com

#ifndef FACE_HPP
#define FACE_HPP

// warpit namespace
namespace warpit {

    // Face data structure.
    class Face
    {
    public:
        // Constructor.
        Face() { v[0] = v[1] = v[2] = v[3] = f[0] = f[1] = f[2] = f[3] = -1; }

        // Elements.
        int v[4]; // indices of the vertices
        int f[4]; // indices of the neighboring faces
    };

} // namespace warpit

#endif // FACE_HPP
