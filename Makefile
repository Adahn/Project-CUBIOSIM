SRC=src/repressilator.cc 
OBJ=repressilator.o  

CC=g++ #nvcc
CFLAGS=-Wall #-O2 -std=c++11

INC=-I./inc 
LIB=-L./lib 
LDFLAGS=$(INC) $(LIB)

# all
.PHONY: all
all: ./bin/repressilator

%.o: src/%.cc
	$(CC) $(LDFLAGS) $(CFLAGS) -c -o $@ $<

./bin/repressilator: $(OBJ)
	$(CC) $(LDFLAGS) $(CFLAGS)  -o $@ $^


# clean
.PHONY: clean
clean:
	rm -f *.o
	rm -f bin/*


exec:
	time ./bin/repressilator


