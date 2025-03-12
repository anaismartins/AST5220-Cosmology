# Hans A. Winther (2020) (hans.a.winther@gmail.com)

SHELL := /bin/bash

# Set compiler (use =c++17 if you have this availiable)
CC = g++ -std=c++17

# Paths to GSL library
INC  = -I/${HOME}/local/include
LIBS = -L/${HOME}/local/lib -lgsl -lgslcblas

#=======================================================
# Options
#=======================================================
OPTIONS = 

# Add bounds checking
OPTIONS += -D_GLIBCXX_DEBUG

# Show warnings if atempting to evaluate a spline out of bounds
OPTIONS += -D_SPLINE_WARNINGS_ON

# Show info about the solution as we integrate
# OPTIONS += -D_FIDUCIAL_VERBOSE_ODE_SOLVER_TRUE

# Add OpenMP parallelization
OPTIONS += -D_USEOPEMP
CC += -fopenmp

#=======================================================

C = -O3 -g $(OPTIONS)

#=======================================================

SRC_DIR=src/cpp/
BUILD_DIR=build/
BIN_DIR=bin/
TARGETS := $(BIN_DIR)cmb
all: $(TARGETS)

# OBJECT FILES
OBJS = $(BUILD_DIR)Main.o $(BUILD_DIR)Utils.o $(BUILD_DIR)BackgroundCosmology.o $(BUILD_DIR)RecombinationHistory.o $(BUILD_DIR)Perturbations.o $(BUILD_DIR)PowerSpectrum.o $(BUILD_DIR)Spline.o $(BUILD_DIR)ODESolver.o

# DEPENDENCIES
Main.o                  : BackgroundCosmology.h RecombinationHistory.h Perturbations.h PowerSpectrum.h
Spline.o                : Spline.h
ODESolver.o             : ODESolver.h
Utils.o                 : Utils.h Spline.h ODESolver.h
BackgroundCosmology.o   : BackgroundCosmology.h Utils.h Spline.h ODESolver.h
RecombinationHistory.o  : RecombinationHistory.h BackgroundCosmology.h
Perturbations.o         : Perturbations.h BackgroundCosmology.h RecombinationHistory.h
PowerSpectrum.o         : PowerSpectrum.h BackgroundCosmology.h RecombinationHistory.h Perturbations.h
Examples.o              : Utils.h Spline.h ODESolver.h

examples: Examples.o Utils.o Spline.o ODESolver.o
	${CC} -o $@ $^ $C $(INC) $(LIBS)

$(BIN_DIR)%: $(OBJS)
	${CC} -o $@ $^ $C $(INC) $(LIBS)

$(BUILD_DIR)%.o: $(SRC_DIR)%.cpp
	${CC} -c -o $@ $< $C $(INC) 

clean:
	rm -rf $(TARGETS) $(BUILD_DIR)*.o

python-lib: $(OBJS)
	${CC} -shared -fPIC -o $(BIN_DIR)libcmb.so $(SRC_DIR)BackgroundCosmology.cpp $(SRC_DIR)ODESolver.cpp $(SRC_DIR)Utils.cpp $(SRC_DIR)Spline.cpp $(LIBS)

.PRECIOUS: $(BUILD_DIR)%.o