########################################################################################################
# General parameters of compilations 
CXX=g++
CXX_FLAGS = -std=c++14 -Wall -O3 -I /ebio/abt6/yvoichek/smallproj/prefix/include/ -I ./include/ -pthread -msse4.2  # -march=native 
#-g -mno-tbm
LDFLAGS :=  -L/ebio/abt6/yvoichek/smallproj/prefix/lib -lstdc++ -lboost_program_options
########################################################################################################

BUILD_DIR = build
BIN_DIR = bin

########################################################################################################
# Define object to be created
SRC_H_OBJ = $(wildcard src/*.h)
SRC_CPP_OBJ = $(subst .h,.cpp,$(SRC_H_OBJ))
SRC_CPP_KMC_OBJ = $(wildcard include/KMC/kmc_api/*.cpp)
OBJ = $(SRC_CPP_OBJ:src/%.cpp=$(BUILD_DIR)/%.o) 
OBJ_KMC = $(SRC_CPP_KMC_OBJ:include/KMC/kmc_api/%.cpp=$(BUILD_DIR)/%.o)
########################################################################################################
# Define BIN to be created
SRC_CPP_BIN = $(filter-out $(SRC_CPP_OBJ),$(wildcard src/*.cpp))
BIN = $(filter-out associate_kmers_with_phenotypes,$(SRC_CPP_BIN:src/%.cpp=%))
########################################################################################################
# Define dependencies - Gcc/Clang will create these .d files containing dependencies.
DEP = $(OBJ:%.o=%.d)


########################################################################################################
# Define rules

all: $(BIN) associate_kmers_with_phenotypes
	@echo "Built everything"

associate_kmers_with_phenotypes : $(OBJ) $(OBJ_KMC)
	mkdir -p bin
	@echo $^
	$(CXX) $(CXX_FLAGS)  $^ src/$(notdir $@).cpp  -o bin/$@ $(LDFLAGS)

$(BIN) : $(OBJ_KMC) $(OBJ)
	@echo $^
	mkdir -p bin
	$(CXX) $(CXX_FLAGS)  $^ src/$(notdir $@).cpp  -o bin/$@

# Include all .d files
-include $(DEP)

.SECONDEXPANSION:
$(OBJ_KMC) : $$(patsubst %.o,%.cpp,include/KMC/kmc_api/$$(notdir $$@))
	mkdir -p $(@D)
	@echo $<
	@echo $@
	$(CXX) $(CXX_FLAGS) -MMD -c $< -o $@

# Build target for every single object file.
# The potential dependency on header files is covered
# by calling `-include $(DEP)`.
$(OBJ) : $$(patsubst %.o,%.cpp,src/$$(notdir $$@)) $(OBJ_KMC)
	mkdir -p $(@D)
	# The -MMD flags additionaly creates a .d file with
	# # the same name as the .o file
	$(CXX) $(CXX_FLAGS) -MMD -c $< -o $@
.PHONY : clean
clean :
	# This should remove all generated files.
	-rm $(BUILD_DIR) $(BIN_DIR) -Rf
########################################################################################################
