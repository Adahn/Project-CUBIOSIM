CXX=g++ 
CXXFLAGS=-std=c++11 -Wall 
CPPFLAGS=

EXEC=main
SRC=$(wildcard *.cc)
OBJ=$(SRC:.cc=.o)

all: $(EXEC)

main: $(OBJ)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(CPPFLAGS)

%.o: %.c
	$(CXX) -o $@ -c $< $(CXXFLAGS)

clean:
	rm -f *.o
	rm -f *~ 
	rm solution*

mrproper: clean
	rm -f $(EXEC)
