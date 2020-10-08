Compile and build:
------------------

Instructions for Ubuntu:
(1) To run this code you need Eigen, OpenGL, and Lapack. Lapack is usually preinstalled on your computer.
Eigen can be found at http://eigen.tuxfamily.org/index.php?title=Main_Page

(2) Create a file named "Makefile". The content of the "Makefile" should be the same "Makefile_Sample" except that you will need to change the path to eigen from "/usr/local/include/eigen3/" to the location of eigen in your system.

(3) Open a terminal, "cd" to this folder and run the command "make" (without the quotes).

(4) To start the simulation, run the command "./simDER option.txt" (without the quotes). More on option.txt later.

For Mac, you need the following modifications:

(1) Change the " #include <GLUT/glut.h> " in "main.cpp"
(2) Use "Makefile_Sample_OSX" instead of "Makefile_Sample".

Physical parameters:
------------------
(1) Angular velocity: a text file, e.g. omega.txt, should contain two columns: one column is the time and second column is the angular velocity (rpm). This is essentially the time series of the angular velocity.
(2) The name of this input file should be specified in "option.txt" file under the option name "input-file". This option.txt should be specified while running the simulation, i.e. ./simDER option.txt
(3) You can edit the parameters of the simulation by editing "option.txt" file. You can also specify an option using the following syntax:
./simDER option.txt -- option_name option_value
Eample: ./simDER option.txt -- RodLength 0.2
(4) Details on the options (we use SI units): 
    "RodLength" is the contour length of the helix.
    "helixradius" is the radius of the helix.
    "helixpitch" is the pitch of the helix.
    "headSize" is the radius of the bacterial head.
    "rodRadius" is the cross-sectional radius of the flagellum.
    "youngM" is the young's modulus.
    "Poisson" is the Poisson ratio.
    "deltaTime" is the time step size.
    "totalTime" is the time at which the simulation ends.
    "tol" and "stol" are small numbers used in solving the linear system. Fraction of a percent, e.g. 1.0e-3, is often a good choice.
    "maxIter" is the maximum number of iterations allowed before the solver quits.
    "density" is the mass per unit volume.
    "gVector" is the vector specifying acceleration due to gravity.
    "viscosity" is the viscosity of the fluid medium.
    "numVertices" is the number of nodes on the rod.
    "render" (0 or 1) indicates whether OpenGL visualization should be displayed.
    "saveData" (0 or 1) indicates whether the location of the head should be saved in "datafiles/" folder (this folder will be created by the program).
    "input-file" is the name of the file that provides the time series of angular velocity (as mentioned above).
    "epsilon" parameter is used Regularized Stokeslet Segment (RSS) model.
    "use-RSS" is a boolean variable. If true, RSS will be used. If false, RFT will be used. See "Propulsion of microorganisms by a helical flagellum" for expressions of eta_par and eta_per in Gray and Hancock model.
    "include-contact" is a boolean variable. If true, contact will be included. Otherwise, contact will be neglected.
