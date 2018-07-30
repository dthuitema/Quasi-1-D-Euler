This directory contains the files for a quasi-1D approach to solving the Euler equations in a converging-diverging nozzle. Contents:

- makefile
- main.cpp (execution file)
- euler_functions.cpp (functions)
- euler_functions.h (header file)
- plotting.py (python script for plotting)


Dependencies/Libraries:
The main code is contained in main.cpp. For the main code, the following libraries are required:
-Eigen

For the python plotting script:
-Python
-Numpy
-Matplotlib

Outputs:
- Data file containing geometry
- Data file containing state variables for each CV
- Data file containing density residual

How to Compile:
- Create new folder with the following programs in it:
	-main.cpp (main file)
	-euler_functions.cpp (functions)
	-euler_functions.h (function header file)
	-plotting.py (plotting python script)
- Open main.cpp
- Compile with: :compile
- then: make -k
- It will then use the make file to compile

How to Run:
- Execute the code from main.cpp
- Will ask for implicit or explicit solution (enter 1 for implicit, 2 for explicit)
- Enter pressure ratio
- Enter number of grid points/control volumes
- Enter epsilon for scalar dissipation for the Roe Flux Splitting scheme
- Program will then execute
- To plot, open and run the python script to plot
- Python script will ask for the name of the figures, and will display and save three figures:
	- Plot of the Mach number distribution
	- Plot of the pressure distribution (P/Pt)
	- Plot of the density residual

