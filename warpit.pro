#---------------------------------------#
#-Project created by QtCreator 04.02.15-#
#-Author: Dmitry Anisimov (c)-----------#
#---------------------------------------#

# Options.
QT += core gui opengl

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

# Choose between release and debug.

# CONFIG += qt release
CONFIG += warn_on debug

# QMAKE_CXXFLAGS_RELEASE += -O2 -O3 -DNDEBUG -fno-strict-overflow
QMAKE_CXXFLAGS_DEBUG += -O2 -O3 -Wall -Wextra

TARGET = warpit

TEMPLATE = app

# The global path to the project.

# MAC OS
DEFINES += PROJECT_PATH=\"\\\"$$/Users/danston/Documents/github/warpit\\\"\"

# WINDOWS
# DEFINES += PROJECT_PATH=\"\\\"C:\\\Users\\\your_username\\\path_to\\\warpit\\\"\"

# Paths to local includes.

# OS
INCLUDEPATH += /Users/danston/Documents/github/warpit/inc
INCLUDEPATH += /Users/danston/Documents/github/warpit/inc/support

# WINDOWS
# INCLUDEPATH += C:\Users\your_username\path_to\warpit\inc
# INCLUDEPATH += C:\Users\your_username\path_to\warpit\inc\support

# Sources.
SOURCES += src/main.cpp       \
           src/mainwindow.cpp \
           src/qtviewer.cpp   \
           src/program.cpp    \
           src/triangle.cpp

# Headers.
HEADERS += inc/mainwindow.hpp           \
           inc/qtviewer.hpp             \
           inc/program.hpp              \
           inc/support/vector2d.hpp     \
           inc/support/face.hpp         \
           inc/support/vertex.hpp       \
           inc/support/halfedge.hpp     \
           inc/support/trimesh.hpp      \
           inc/support/polygon.hpp      \
           inc/support/texture.hpp      \
           inc/support/triangle.hpp     \
           inc/support/triangulator.hpp \
           inc/support/utils.hpp

# Forms.
FORMS += ui/mainwindow.ui
