//* Authors: Teseo Schneider and Dmitry Anisimov (c) 2016-2017
//* With any questions please contact: danston@ymail.com

#ifndef UTILS_HPP
#define UTILS_HPP

#ifndef NDEBUG
#define CHECK_GL_VALID_STATE() warpit::utils::checkGLValidState(__FILE__, __LINE__);
#else
#define CHECK_GL_VALID_STATE() //
#endif

#ifdef WIN32
#include <GL/glu.h>
#else
#include <glu.h>
#endif

// STL includes.
#include <string>

// warpit namespace
namespace warpit {

    // The util class.
    class utils
    {
    public:
        // Check an OpenGL validity state.
        void static checkGLValidState(const std::string &file, const int line)
        {
            GLenum error = glGetError();
            if(error != GL_NO_ERROR) std::cout << "glx error: " << file << ":" << line << " > " << gluErrorString(error) << " - " << error << std::endl;
        }
    };

} // namespace warpit

#endif // UTILS_HPP
