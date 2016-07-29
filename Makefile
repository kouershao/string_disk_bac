all : test

CC=g++
WARNING=-Wall -Wunused-variable 
LDFLAGS=-pthread -lfftw3 -lm-2.12 -lgsl -lgslcblas -llbfgs -std=c++11 -O3 \
#-Wl,--rpath=./lib \
#-Wl,--dynamic-linker=/lib64/ld-2.12.so


src = main.cpp myString.cpp myNode.cpp
inc = myString.h myNode.h

test : $(src) $(includes)
	$(CC) -o $@ $(src) $(LDFLAGS) #$(WARNING)

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
