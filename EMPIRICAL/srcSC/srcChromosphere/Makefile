include ../../../../Makefile.def
include ../../../../Makefile.conf

OBJ_FILES = chromosphere.cpp.o 
SRC_FILES = chromosphere.cpp
INC_DIRS := . include/
# Add a prefix to INC_DIRS. So moduleA would become -ImoduleA. GCC understands this -I flag
INC_FLAGS := $(addprefix -I,$(INC_DIRS))

CPPFLAGS := $(INC_FLAGS)

%.cpp.o: %.cpp
	$(CXX) --std=c++11 -g $(CFLAGS) $(CPPFLAGS) -c $< -o $@

chromosphere: $(OBJ_FILES)
	echo finished

clean:
	rm -f chromosphere.cpp.o
