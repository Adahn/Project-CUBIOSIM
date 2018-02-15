#########################################################################
# This README refers only to the directory rk4 in which the sequential
# method is implemented (not parallel!)
#########################################################################

# ------------------- File structure ---------------------
This directory consists of following files/directories:

    - Makefile: Configuration for compilation
                In line: all: ./bin/repr_rk45
                you can modify the name of the compiled program 'repr_rk45'

    - bin: The compiled executable will be stored in this directory
           If you want to change the output directory, you can modify it in the Makefile

    - inc: Includes following header files:
        - chCommandLine.h: TODO LUCAS
        - chTimer.hpp: TODO LUCAS
        - chTimer.h: TODO LUCAS
        - runge_kutta.hh: Solvers for Runge-Kutta methods RK4 and RK45 (normally no need to modify)
        - repressilator.hh: Class Repressilator describing the model: dYdt = S*R - d*Y
        - utility.hh: Helping functions

    - lib: Normally empty, just added for completeness

    - src:
        - repr_rngkutta.cc: Main function calling the solver (adapt parameters here!)

    - script.sh: Bash script to test program for several input sizes (good for testing!)




# ------------------- Getting started ----------------------

1) Compile the program: when you are in the directory rk4, tape the following
                        command:
                        > make

                        The program gets compiled automatically

                        Remark: If you add or delete some files, don't forget to
                                adapt the Makefile
                        Remark: If you want to delete compiled files
                                > make clean

2) Call program: TODO LUCAS


# -------------------- Modify programs -------------------------
If you want to modify the programs/parameters, consider the comments of this form:

// !! this is a sample comment !!

Probably you will essentially modify the main function in the file src/repr_rngkutta.cc,
check these type of comments and you will see where you can change parameters/
initial values etc.

If you want to modify the model of the Repressilator, you can do this in the file
inc/repressilator.hh 

Normally you don't need to modify the solver Runge-Kutta methods in the file
inc/runge_kutta.hh , except you optimise them
