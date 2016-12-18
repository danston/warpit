//* Authors: Dmitry Anisimov
//* main.cpp
//* With any questions please contact: danston@ymail.com

// Local includes.
#include "mainwindow.hpp"

// Qt includes.
#include <QApplication>

// Main function.
int main(int argc, char *argv[])
{
    QApplication app(argc, argv);
    warpit::MainWindow window;
    window.show();

    return app.exec();
}
