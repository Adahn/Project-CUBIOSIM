SRC=src/repr_rngkutta.cc
OBJ=repr_rngkutta.o

CC=g++ #nvcc
CFLAGS=-Wall -std=c++11

INC=-I./inc 
LIB=-L./lib 
LDFLAGS=$(INC) $(LIB)

.PHONY: all
all: ./bin/repr_rk4

%.o: src/%.cc
	$(CC) $(LDFLAGS) $(CFLAGS) -c -o $@ $<

./bin/repr_rk4: $(OBJ)
	$(CC) $(LDFLAGS) $(CFLAGS)  -o $@ $^


.PHONY: clean
clean:
	rm -f *.o
	rm -f bin/*

exec:
	time ./bin/repr_rk4

