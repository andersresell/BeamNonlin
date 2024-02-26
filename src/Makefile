
#*******************************************************************************
# Default build settings. Override values in .configure
#*******************************************************************************
BUILD_TYPE ?= release
SINGLE_PRECISION ?= 0

BUILD_DIR_DEBUG = $(BeamNonlinHome)/build_debug
BUILD_DIR_RELEASE = $(BeamNonlinHome)/build_release


CXX_FLAGS_COMMON = -Wall -Wno-unused-function -fmax-errors=1  -fopenmp -lyaml-cpp
CXX_FLAGS = $(CXX_FLAGS_COMMON)

ifeq ($(BUILD_TYPE), debug)
	BUILD_DIR = $(BUILD_DIR_DEBUG)
	CXX_FLAGS += -g
	#CXX_FLAGS += -D_GLIBCXX_DEBUG #_GLIBCXX_DEBUG enables range checking for [] operator for std vector
else ifeq ($(BUILD_TYPE), release)
	BUILD_DIR = $(BUILD_DIR_RELEASE)
	CXX_FLAGS +=  -DNDEBUG -O3 -march=native 
else
    $(error Unsupported build type specified: $(BUILD_TYPE). Valid types are "debug" and "release")
endif



EXE = $(BUILD_DIR)/BeamNonlin
INC_DIR = $(BeamNonlinHome)/include


#Add source files here
SRC = BeamNonlin.cpp \
	  Utils.cpp \
	  Solver.cpp \
	  SolverUtils.cpp \
	  Geometry.cpp \
	  Config.cpp \
	  InputParser.cpp

OBJ = $(patsubst %.cpp, $(BUILD_DIR)/%.o, $(SRC)) 

# GLIBCXX_DEBUG doesn't work with yaml-cpp
#$(BUILD_DIR)/InputParser.o: CXX_FLAGS = $(CXX_FLAGS_COMMON) -g 



all: prep $(EXE)

$(EXE): $(OBJ) 
	g++ -o $(EXE) $^ $(LDLIBS) $(CXX_FLAGS)




#*******************************************************************************
# Build C++ object files
#*******************************************************************************
#The -MMD option creates a list of dependencies, and it's included later
$(BUILD_DIR)/%.o: %.cpp  #adding .configure as a dependency, so that compilation is redone if this file is changed
	g++ -MMD -c -o $@ $< $(CXX_FLAGS)

-include $(OBJ:.o=.d)



prep:
	mkdir -p $(BUILD_DIR_DEBUG)
	mkdir -p $(BUILD_DIR_RELEASE)

.PHONY: clean 

clean: 
	-rm -rf $(BUILD_DIR_DEBUG)
	-rm -rf $(BUILD_DIR_RELEASE)





