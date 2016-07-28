all:test

CC=g++
DEBUG=-g
OPT=-O3
WARNING=-Wall -Wunused-variable 
LDFLAGS=-pthread -lfftw3 -lm -lgsl -lgslcblas -llbfgs -std=c++11

src = main.cpp string.cpp node.cpp
inc = myString.h node.h

test : $(src) $(includes)
	$(CC) -o $@ $(src) $(LDFLAGS) #$(WARNING)

#disk1: disk1.o
	#$(CC)  $(LDFLAGS) -o $@ $^ 

.PHONY:
	clean run

clean:
	-rm test
run: test
	./test a

