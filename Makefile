main: Read.file.o
	g++ -Wall -g Read.file.o -o main

Read.file.o: Read.file.cc Read.file.hh
	g++ -g -Wall -c Read.file.cc

clean:
	rm -rf *.o

mrproper: clean
	rm -rf *
