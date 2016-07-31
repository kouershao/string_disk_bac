all : endfixed

CC=g++
WARNING=-Wall -Wunused-variable 
LDFLAGS=-pthread -lfftw3 -lm -lgsl -lgslcblas -llbfgs -std=c++11 -O3

src = main.cpp myString.cpp myNode.cpp
inc = myString.h myNode.h

endfixed : $(src) $(includes)
	$(CC) -o $@ $(src) $(LDFLAGS) #$(WARNING)

#disk1: disk1.o
	#$(CC)  $(LDFLAGS) -o $@ $^ 

.PHONY:
	clean run

clean:
	-rm endfixed
run: endfixed
	./endfixed e

