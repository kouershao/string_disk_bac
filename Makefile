all : timestepfixed

CC=gcc
CFLAGS=-std=c++11 -O3 #-Wall -Wunused-variable 
LDFLAGS=-lstdc++ -lm -lfftw3 -lgsl -lgslcblas -llbfgs
FLAGS=$(CFLAGS) $(LDFLAGS)

src = main.cpp myString.cpp myNode.cpp
inc = myString.h myNode.h

timestepfixed : $(src) $(includes)
	$(CC) -o $@ $(src) $(FLAGS)

#disk1: disk1.o
	#$(CC)  $(LDFLAGS) -o $@ $^ 

.PHONY:
	clean run

clean:
	-rm timestepfixed
run: timestepfixed
	./timestepfixed t

