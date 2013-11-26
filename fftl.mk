CXX = clang++ -std=gnu++11 -D__STRICT_ANSI__
CXXFLAGS = -Wall -g -D FFTL_UNIT_TEST
INCLUDES = 
LIBS =
OBJS = fftl.o
PROGRAM = fftl.out

all:$(PROGRAM)

$(PROGRAM): $(OBJS)
	$(CXX) $(CXXFLAGS) $^ $(INCLUDES) $(LIBS) -o $(PROGRAM)

.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(LIBS) -c $<

.PHONY: clean
clean:
	rm -f *o $(PROGRAM)
