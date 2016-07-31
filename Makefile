all : test

CC=gcc
CFLAGS=-std=c++11 -O3#-Wall -Wunused-variable 
LDFLAGS=-lstdc++ -lfftw3 -lgsl -lgslcblas -llbfgs
FLAGS=$(CFLAGS) $(LDFLAGS)


src = main.cpp myString.cpp myNode.cpp
inc = myString.h myNode.h

test : $(src) $(includes) Makefile
	$(CC) -o $@ $(src) $(FLAGS)

#disk1: disk1.o
	#$(CC)  $(LDFLAGS) -o $@ $^ 

.PHONY:
	clean run ls

clean:
	-rm test
run: test
	[[ ! -f "log" ]] || mv log log.bak
	./test a >log 2>&1 &
	tail -f log
