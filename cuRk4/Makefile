SRC=src/repr_rngkutta.cu
OBJ=repr_rngkutta.o

CC=nvcc
CFLAGS=-std=c++11 -g -G -lineinfo -arch=compute_35
#CFLAGS=-std=c++11 -g -G -lineinfo -arch=compute_20

INC=-I./inc 
LIB=-L./lib 
LDFLAGS=$(INC) $(LIB)

.PHONY: all
all: ./bin/repr_rk4

#%.o: src/%.cu
#	$(CC) $(LDFLAGS) $(CFLAGS) -c -o $@ $<

repr_rngkutta.o: src/repr_rngkutta.cu inc/runge_kutta.hh inc/utility.hh
	$(CC) $(LDFLAGS) $(CFLAGS) -c -o $@ $<

./bin/repr_rk4: $(OBJ)
	$(CC) $(LDFLAGS) $(CFLAGS)  -o $@ $^


.PHONY: clean
clean:
	rm -f *.o
	rm -f bin/*

exec:
	time ./bin/repr_rk4

