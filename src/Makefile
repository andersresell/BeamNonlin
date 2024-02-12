
#*******************************************************************************
# Default build settings. Override values in .configure
#*******************************************************************************
BUILD_TYPE ?= release
SINGLE_PRECISION ?= 0

BUILD_DIR_DEBUG = $(DrillSimHome)/build_debug
BUILD_DIR_RELEASE = $(DrillSimHome)/build_release


CXX_FLAGS =  -Wall -Wno-unused-function -fmax-errors=1 
CXX_FLAGS_NO_GLIBCXX_DEBUG = $(CXX_FLAGS)

ifeq ($(BUILD_TYPE), debug)
	BUILD_DIR = $(BUILD_DIR_DEBUG)
	CXX_FLAGS_NO_GLIBCXX_DEBUG += -g
	CXX_FLAGS += -g -D_GLIBCXX_DEBUG #_GLIBCXX_DEBUG enables range checking for [] operator for std vector
else ifeq ($(BUILD_TYPE), release)
	BUILD_DIR = $(BUILD_DIR_RELEASE)
	CXX_FLAGS +=  -DNDEBUG -O3 -march=native 
	CUDA_FLAGS += -DNDEBUG -O3
else
    $(error Unsupported build type specified: $(BUILD_TYPE). Valid types are "debug" and "release")
endif



EXE = $(BUILD_DIR)/BeamNonlin
INC_DIR = $(DrillSimHome)/include

#Add source files here
SRC = BeamNonlin.cpp \


OBJ = $(patsubst %.cpp, $(BUILD_DIR)/%.o, $(SRC)) 

all: prep $(EXE)

$(EXE): $(OBJ) $(OBJ_CUDA)
	g++ -o $(EXE) $^ $(LDLIBS) $(CXX_FLAGS)

#*******************************************************************************
# Build C++ object files
#*******************************************************************************
#The -MMD option creates a list of dependencies, and it's included later
$(BUILD_DIR)/%.o: %.cpp $(DrillSimHome)/.configure #adding .configure as a dependency, so that compilation is redone if this file is changed
	g++ -MMD -c -o $@ $< $(CXX_FLAGS)



prep:
	mkdir -p $(BUILD_DIR_DEBUG)
	mkdir -p $(BUILD_DIR_RELEASE)

.PHONY: clean 

clean: $(DrillSimHome)/.configure
	-rm -rf $(BUILD_DIR_DEBUG)
	-rm -rf $(BUILD_DIR_RELEASE)





