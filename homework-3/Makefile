#define compiler
CC = g++

# directories 
INC_DIR = ./include
SRC_DIR = ./src
OBJ_DIR = ./obj

# *_FILES use "wildcards" to create an automatic list of all the header (.hpp), source (.cpp), and object (.o) files based on the directory.
INC_FILES := $(wildcard $(INC_DIR)/*.hpp)
SRC_FILES := $(wildcard $(SRC_DIR)/*.cpp)
OBJ_FILES := $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(SRC_FILES))

# -I$(INC_DIR) tells GCC to include files in ./include
CFLAGS = -I$(INC_DIR)

# this DEPS says that if the files in ./include change, things should get recompiled. 
DEPS = $(INC_FILES)

# Define targets - first target is always run
all: main

# @ and < are "automatic variables"
# @ stands for the target
# ^ stands for all prerequisites

main: $(OBJ_FILES) main.cpp
	$(CC) $^ $(CFLAGS) -o $@

# "@" refers to the target (left hand side)
# We use "<", which stands for the first prerequisite. 
# If we use "^", we build all prereqs, including $(DEPS).
# If we use "^", we don't build the $(DEPS) prereq - these are just the headers, and we include them so if they change, the makefile gets updated.
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(DEPS)
	$(CC) -c $< -o $@ $(CFLAGS)

clean:
	rm -f main a.out $(OBJ_DIR)/*.o *~ *.o 


