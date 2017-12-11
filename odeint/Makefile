SRC=src/repressilator_n.cc 
OBJ=repressilator_n.o  

CC=g++ #nvcc
CFLAGS=-Wall -std=c++11

INC=-I./inc 
LIB=-L./lib 
LDFLAGS=$(INC) $(LIB)

.PHONY: all
all: ./bin/repressilator_n

%.o: src/%.cc
	$(CC) $(LDFLAGS) $(CFLAGS) -c -o $@ $<

./bin/repressilator_n: $(OBJ)
	$(CC) $(LDFLAGS) $(CFLAGS)  -o $@ $^


.PHONY: clean
clean:
	rm -f *.o
	rm -f bin/*

exec:
	time ./bin/repressilator_n

