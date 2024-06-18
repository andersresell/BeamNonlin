
#*******************************************************************************
# Default build settings. Override values in .configure
#*******************************************************************************

#*******************************************************************************
# Some usefult makefile examples:
#*******************************************************************************
#override flag: (use with caution)
#$(BUILD_DIR)/InputParser.o: CXX_FLAGS = $(CXX_FLAGS_COMMON) -g 



BUILD_TYPE ?= release
SINGLE_PRECISION ?= 0

BUILD_DIR_DEBUG = $(BeamNonlinHome)/build-debug
BUILD_DIR_RELEASE = $(BeamNonlinHome)/build-release


CXX_FLAGS_COMMON = -Wall -Wno-unused-function -fmax-errors=1 -lyaml-cpp -fopenmp 
ifeq ($(SINGLE_PRECISION), 1)
	CXX_FLAGS_COMMON += -DSINGLE_PRECISION
endif
#CXX_FLAGS_COMMON += -fno-math-errno -fno-trapping-math #use these? have heard it's possibly unsafe to use -ffast-math
CXX_FLAGS = $(CXX_FLAGS_COMMON)

ifeq ($(BUILD_TYPE), debug)
	BUILD_DIR = $(BUILD_DIR_DEBUG)
	CXX_FLAGS += -g -D_GLIBCXX_DEBUG  #_GLIBCXX_DEBUG enables range checking for [] operator for std vector
	YAML_CPP_INSTALL_DIR = $(HOME)/yaml-cpp/debug_install
#specify includes and library paths for yaml-cpp if the non default debug build is used
	CXX_FLAGS += -I$(YAML_CPP_INSTALL_DIR)/include -L$(YAML_CPP_INSTALL_DIR)/lib 
else ifeq ($(BUILD_TYPE), release)
	BUILD_DIR = $(BUILD_DIR_RELEASE)
	CXX_FLAGS += -DNDEBUG  -O3 -march=native #-flto
else
    $(error Unsupported build type specified: $(BUILD_TYPE). Valid types are "debug" and "release")
endif


EXE = $(BUILD_DIR)/BeamNonlin
INC_DIR = $(BeamNonlinHome)/include


#Add source files here
SRC = BeamNonlin.cpp \
	  Utils.cpp \
	  Solver.cpp \
	  Numerics.cpp \
	  SolverUtils.cpp \
	  Geometry.cpp \
	  Config.cpp \
	  InputParser.cpp \
	  Borehole.cpp \
	  HoleContact.cpp \
	  CrisfieldBeam.cpp \
	  BattiniBeam.cpp \
	  Testing.cpp

OBJ = $(patsubst %.cpp, $(BUILD_DIR)/%.o, $(SRC)) 

ifdef USER_DEFINED_EXTERNAL_FORCE_FILE
SRC_USER = $(USER_DEFINED_EXTERNAL_FORCE_FILE)
OBJ_USER = $(BUILD_DIR)/$(notdir $(SRC_USER:.cpp=.o))
#FLAGS_USER = -DUSER_DEFINED_EXTERNAL_FORCE_FILE 
else
SRC_USER = UserFunction.cpp
OBJ_USER = $(patsubst %.cpp, $(BUILD_DIR)/%.o, $(SRC_USER))
#FLAGS_USER =
endif

OBJ += $(OBJ_USER) 

$(info obj= $(OBJ))


$(info src_user= $(SRC_USER))

$(info obj_user= $(OBJ_USER))

all: prep $(EXE)

$(EXE): $(OBJ) 
	g++ -o $(EXE) $^ $(CXX_FLAGS)


#*******************************************************************************
# Build C++ object files
#*******************************************************************************
#The -MMD option creates a list of dependencies, and it's included later
$(BUILD_DIR)/%.o: %.cpp  #adding .configure as a dependency, so that compilation is redone if this file is changed
	g++ -MMD -c -o $@ $< $(CXX_FLAGS)


$(BUILD_DIR)/$(SRC:.cpp=.o): $(SRC)  
	g++ -MMD -c -o "$@" "$<" $(CXX_FLAGS)  

$(BUILD_DIR)/$(notdir $(SRC_USER:.cpp=.o)): $(SRC_USER)  
	g++ -MMD -c -o "$@" "$<" $(CXX_FLAGS)  #$(FLAGS_USER)


-include $(OBJ:.o=.d)

prep:
	mkdir -p $(BUILD_DIR_DEBUG)
	mkdir -p $(BUILD_DIR_RELEASE)

.PHONY: clean 

clean: 
	-rm -rf $(BUILD_DIR_DEBUG)
	-rm -rf $(BUILD_DIR_RELEASE)





