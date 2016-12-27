//* Author: Dmitry Anisimov (c) 2016-2017
//* With any questions please contact: danston@ymail.com

#ifndef TEXTURE_HPP
#define TEXTURE_HPP

// STL includes.
#include <cassert>
#include <iostream>

// Local includes.
#include "utils.hpp"

// Qt includes.
#include <QGLFunctions>
#include <QImageReader>

// warpit namespace
namespace warpit {

    // The texture class.
    class Texture
    {
    public:
        // Constructors.
        Texture() : loaded(false), path(), qglFunctions() { qglFunctions.initializeGLFunctions(); }
        Texture(const std::string &_path) : loaded(false), path(_path), qglFunctions() { qglFunctions.initializeGLFunctions(); }

        // Destructor.
        ~Texture()
        {
            if(!loaded) return;
            glDeleteTextures(1, &name);

            CHECK_GL_VALID_STATE();
        }

        // Bind the texture.
        inline void bind()
        {
            assert(loaded);

            qglFunctions.glActiveTexture(GL_TEXTURE0);

            glEnable(target);
            glBindTexture(target, name);

            CHECK_GL_VALID_STATE();
        }

        // Unbind the texture.
        inline void unbind()
        {
            glDisable(target);

            CHECK_GL_VALID_STATE();
        }

        // Set a texture type.
        inline void setTextureType(const GLenum &type)
        {
            if(type == GL_TEXTURE_1D) {
                target = GL_TEXTURE_1D;
                finalize1D();
            }

            if(type == GL_TEXTURE_2D) {
                target = GL_TEXTURE_2D;
                finalize2D();
            }
        }

        // Return the texture width.
        inline int width() const
        {
            assert(loaded);
            return _width;
        }

        // Return the texture height.
        inline int height() const
        {
            assert(loaded);
            return _height;
        }

        // Return the original image.
        inline const QImage& getImage() const
        {
            assert(loaded);
            return _img;
        }

    private:
        // Global variables.
        bool loaded;
        std::string path;

        int _width;
        int _height;

        GLenum target;
        GLuint name;

        QImage _img;

        QGLFunctions qglFunctions;

        // Finilize the program with 1D texture.
        void finalize1D()
        {
            QImageReader reader(path.c_str());

            const bool read = reader.read(&_img);

            if(!read) {
                std::cout << "Failed to read: " << path.c_str() << " with message:" << reader.errorString().toStdString().c_str() << "; " << std::endl;
                assert(read);
                return;
            }

            assert(_img.width() > 0.0);

            QImage img = QGLWidget::convertToGLFormat(_img);

            CHECK_GL_VALID_STATE();

            glGenTextures(1, &name);
            glBindTexture(target, name);

            CHECK_GL_VALID_STATE();

            qglFunctions.glGenerateMipmap(GL_TEXTURE_1D);

            CHECK_GL_VALID_STATE();

            glTexParameteri(target, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
            glTexParameteri(target, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

            glTexParameteri(target, GL_TEXTURE_WRAP_S, GL_CLAMP);
            glTexParameteri(target, GL_TEXTURE_WRAP_T, GL_CLAMP);

            CHECK_GL_VALID_STATE();

            glTexImage1D(target, 0, GL_RGBA, img.height(), 0, GL_RGBA, GL_UNSIGNED_BYTE, img.bits());

            _width  = 1;
            _height = img.height();

            CHECK_GL_VALID_STATE();

            loaded = true;
        }

        // Finilize the program with 2D texture.
        void finalize2D()
        {
            QImageReader reader(path.c_str());

            const bool read = reader.read(&_img);

            if(!read) {
                std::cout << "Failed to read: " << path.c_str() << " with message:" << reader.errorString().toStdString().c_str() << "; " << std::endl;
                assert(read);
                return;
            }

            assert(_img.width() > 0.0);

            QImage img = QGLWidget::convertToGLFormat(_img);

            CHECK_GL_VALID_STATE();

            glGenTextures(1, &name);
            glBindTexture(target, name);

            CHECK_GL_VALID_STATE();

            qglFunctions.glGenerateMipmap(GL_TEXTURE_2D);

            CHECK_GL_VALID_STATE();

            glTexParameteri(target, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
            glTexParameteri(target, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

            glTexParameteri(target, GL_TEXTURE_WRAP_S, GL_REPEAT);
            glTexParameteri(target, GL_TEXTURE_WRAP_T, GL_REPEAT);

            CHECK_GL_VALID_STATE();

            glTexImage2D(target, 0, GL_RGBA, img.width(), img.height(), 0, GL_RGBA, GL_UNSIGNED_BYTE, img.bits());

            _width  = img.width();
            _height = img.height();

            CHECK_GL_VALID_STATE();

            loaded = true;
        }
    };

} // warpit namespace

#endif // TEXTURE_HPP
