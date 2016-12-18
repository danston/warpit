# warpit
Image warper based on generalized barycentric coordinates.

In order to compile the program do:

# PREPARATION:

# MAC OS:

1. Install and do (if not yet done):

a. Install Macports from here: https://www.macports.org, then open terminal and type:
b. sudo port install qt5-mac
c. sudo port install qt5-creator-mac

2. Open warpit.pro and edit these lines:

# MAC OS
DEFINES += PROJECT_PATH=\"\\\"$$/Users/your_username/path_to/warpit\\\"\"

# MAC OS
INCLUDEPATH += /Users/your_username/path_to/warpit/inc
INCLUDEPATH += /Users/your_username/path_to/warpit/inc/support

# WINDOWS:

1. Install and do (if not yet done):
a. Install Qt 5 and Qt Creator from here: http://www.qt.io/download-open-source/
b. When installing choose the version with Open GL 

2. Open warpit.pro and edit these lines:

# WINDOWS
DEFINES += PROJECT_PATH=\"\\\"C:\\\Users\\\your_username\\\path_to\\\warpit\\\"\"

# WINDOWS
INCLUDEPATH += C:\Users\your_username\path_to\warpit\inc
INCLUDEPATH += C:\Users\your_username\path_to\warpit\inc\support

# NOTE: 
In case of WINDOWS, uncomment the corresponding lines and comment the MAC OS related lines!

# COMPILATION:

# MAC OS:

1. First method:

a. Open Terminal
b. Type: cd path_to/warpit/bin
c. Type: qmake ..
d. Type: make -j7
e. Open the app warpit.app

2. Second method:

a. Configure warpit.pro with Qt Creator (just open it)
b. Choose Release vesion in the left bottom corner
c. Build the program (cmd + R)

# WINDOWS:

a. Configure warpit.pro with Qt Creator (just open it)
b. Choose Release version in the left bottom corner
c. Build the program (ctrl + R)

# NOTE:
The program can crash if you choose too many subdivision steps because of lack of the available memory. I am still debugging it.

# CONTACTS:
In case of any questions, please report them here: danston@ymail.com
