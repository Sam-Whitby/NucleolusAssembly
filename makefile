HEADERS  = ./src
SRC_DIR  := ./src
OBJ_DIR  := ./obj
SRC_FILES := $(wildcard $(SRC_DIR)/*.cpp)
OBJ_FILES := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRC_FILES))

FLAGS  = -std=c++11 -O3 -MMD
CXX   := g++

# Default builds both thesis simulation binaries
all: run_nucleolus run_condensate

# Compile all src/*.cpp → obj/*.o
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(FLAGS) -c -o $@ $< -I$(HEADERS)

run_nucleolus: run_nucleolus.o $(OBJ_FILES)
	$(CXX) $(FLAGS) -o run_nucleolus $^ -I$(HEADERS)

run_nucleolus.o: run_nucleolus.cpp
	$(CXX) -c $(FLAGS) run_nucleolus.cpp -I$(HEADERS)

run_condensate: run_condensate.o $(OBJ_FILES)
	$(CXX) $(FLAGS) -o run_condensate $^ -I$(HEADERS)

run_condensate.o: run_condensate.cpp
	$(CXX) -c $(FLAGS) run_condensate.cpp -I$(HEADERS)

run_hier: run_hier.o $(OBJ_FILES)
	$(CXX) $(FLAGS) -o run_hier $^ -I$(HEADERS)

run_hier.o: run_hier.cpp
	$(CXX) -c $(FLAGS) run_hier.cpp -I$(HEADERS)

run_custom: run_custom.o $(OBJ_FILES)
	$(CXX) $(FLAGS) -o run_custom $^ -I$(HEADERS)

run_custom.o: run_custom.cpp
	$(CXX) -c $(FLAGS) run_custom.cpp -I$(HEADERS)

run_polymer: run_polymer.o $(OBJ_FILES)
	$(CXX) $(FLAGS) -o run_polymer $^ -I$(HEADERS)

run_polymer.o: run_polymer.cpp
	$(CXX) -c $(FLAGS) run_polymer.cpp -I$(HEADERS)

clean:
	rm -f *.o $(OBJ_DIR)/*.o *.d $(OBJ_DIR)/*.d \
	      run_nucleolus run_condensate run_hier run_custom run_polymer

-include $(OBJ_FILES:.o=.d)

