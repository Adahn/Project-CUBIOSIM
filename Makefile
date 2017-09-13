SRC=src/example_ode.cc 
OBJ=example_ode.o  

CC=g++ #nvcc
CFLAGS=-Wall #-O2 -std=c++11

INC=-I./inc 
LIB=-L./lib 
LDFLAGS=$(INC) $(LIB)

# all
.PHONY: all
all: ./bin/example_ode

%.o: src/%.cc
	$(CC) $(LDFLAGS) $(CFLAGS) -c -o $@ $<

./bin/example_ode: $(OBJ)
	$(CC) $(LDFLAGS) $(CFLAGS)  -o $@ $^


# clean
.PHONY: clean
clean:
	rm -f *.o


exec:
	time ./bin/example_ode

